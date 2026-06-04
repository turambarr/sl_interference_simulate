% batch_baseband_pss_vote.m
% 目标：
% 1) 从一个多帧基带 IQ 文件中批量检测所有 PSS 同步峰
% 2) 以每个峰作为一帧的起点，执行 SSS 解调，得到 4 个相位候选序列
% 3) 记录每帧候选、做序列相似性分组、剔除离群值
% 4) 对每个簇进行按位多数投票，输出最终结果
%
% 同步端说明：
% - 本脚本使用纯基带 PSS 模板
% - 匹配过程参考 PSSsync_baseband.m
% - 不再做 63 MHz 射频下变频；输入信号视为 zero-IF 基带信号

clear; clc; close all;

%% ===== 自动切换到脚本目录 =====
if ~isdeployed
    currentFilePath = fileparts(mfilename('fullpath'));
    cd(currentFilePath);
    fprintf('已自动切换到脚本所在目录: %s\n', pwd);
end

%% ===== 参数区 =====
% 输入文件：修改为你的多帧基带信号
inFile = 'multiFrame_clean_baseband.iq';
input_format = 'auto';
dat_header_bytes = 0;

% 同步扫描参数
sync_scan_start = 0;
sync_scan_len = 20e6;
min_peak_height = 50;
peak_select_mode = 'all';   % 'all'：保留所有峰；'first'：只用第一个；'strongest'：只用最强

% SSS 解调参数（沿用你现有流程）
read_length = 6992*3 + 1000;
sss_decode_start_idx = 1024 + 48;
target_offset = 10;

fs_source = 409.6e6;
fs_target = 60e6;
sro_ppm  = 0;
N_fft = 1024;
demod_all_subcarriers = true;
decision_rotate_deg = 45;
plot_constellation_enable = true;
plot_constellation_max_frames = 4;

% 基带 CFO 盲估参数（参考 PSSsync_baseband.m）
blind_cfo_enable = true;
blind_cfo_min_len = 4096;
blind_cfo_guard = 1000;
blind_cfo_m = 4;
blind_cfo_fft = 2^17;
blind_cfo_search_range = 200;

% 分组/剔除参数
max_groups = 4;
assign_thresh = 0.18;
outlier_thresh = 0.22;
unstable_ratio_thresh = 0.60;
full_wrong_ber_thresh = 0.40;

% 输出
out_all_csv = 'baseband_all_candidates.csv';
out_group_csv = 'baseband_group_summary.csv';
out_unstable_csv = 'baseband_unstable_positions.csv';
out_ber_csv = 'baseband_group_ber.csv';

%% ===== 1. 生成纯基带 PSS 模板 =====
fprintf('Step 1: 生成纯基带 PSS 模板...\n');
pss_local_bb = local_build_baseband_pss(fs_source);
L_pss = length(pss_local_bb);
fprintf('          PSS 模板长度: %d 点\n', L_pss);

%% ===== 2. 读取多帧基带信号 =====
fprintf('Step 2: 读取多帧基带信号...\n');
fInfo = dir(inFile);
if isempty(fInfo)
    error('找不到文件: %s', inFile);
end

if strcmpi(input_format, 'dat') || (strcmpi(input_format, 'auto') && endsWith(lower(inFile), '.dat'))
    dataBytes = max(0, fInfo.bytes - dat_header_bytes);
else
    dataBytes = fInfo.bytes;
end

totalSamples = floor(dataBytes / 4);
if sync_scan_start < 0 || sync_scan_start >= totalSamples
    error('sync_scan_start 越界: %d / total=%d', sync_scan_start, totalSamples);
end

sync_scan_len = min(round(sync_scan_len), totalSamples - sync_scan_start);
[x_sync, meta_sync] = read_iq_auto_local(inFile, sync_scan_start, sync_scan_len, input_format, dat_header_bytes);
if meta_sync.numSamplesRead <= 0
    error('同步扫描阶段未读取到有效数据。');
end

x_sync = double(x_sync(:));
x_sync = x_sync - mean(x_sync);
x_sync = x_sync / (mean(abs(x_sync)) + eps);

fprintf('          成功读取 %d 个采样点 (约 %.2f ms)\n', length(x_sync), length(x_sync)/fs_source*1000);

%% ===== 3. 能量起振点 + 基带 CFO 盲估 =====
fprintf('Step 3: 基带能量起振点 + CFO 盲估...\n');
pwr = abs(x_sync).^2;
smoothed_pwr = movmean(pwr, 50);
noise_floor = mean(smoothed_pwr(1:min(2000, length(x_sync))));
steady_signal = max(smoothed_pwr);
threshold_pwr = noise_floor + 0.1 * (steady_signal - noise_floor);
energy_start_idx = find(smoothed_pwr > threshold_pwr, 1, 'first');
if isempty(energy_start_idx)
    energy_start_idx = 1;
end

if blind_cfo_enable
    cfo_est_len = min(20000, length(x_sync) - energy_start_idx);
    if cfo_est_len > blind_cfo_min_len
        start_idx = min(length(x_sync), energy_start_idx + blind_cfo_guard);
        end_idx = min(length(x_sync), start_idx + cfo_est_len - 1);
        if end_idx - start_idx + 1 >= blind_cfo_min_len
            x_for_cfo = x_sync(start_idx:end_idx);
            delta_f = local_blind_cfo_estimate(x_for_cfo, fs_source, blind_cfo_m, blind_cfo_fft, blind_cfo_search_range, 0);
            fprintf('          估计到残余 CFO: %.2f kHz\n', delta_f / 1e3);
        else
            delta_f = 0;
            fprintf('          CFO 估计区间过短，跳过。\n');
        end
    else
        delta_f = 0;
        fprintf('          首个有效区域过短，跳过 CFO 盲估。\n');
    end
else
    delta_f = 0;
    fprintf('          已关闭 CFO 盲估。\n');
end

%% ===== 4. 全局 CFO 补偿，得到纯基带信号 =====
fprintf('Step 4: 全局 CFO 补偿...\n');
t_global = (0:length(x_sync)-1).' / fs_source;
baseband_rx = x_sync .* exp(-1j * 2 * pi * delta_f * t_global);

%% ===== 5. 基带滑动归一化互相关扫描 =====
fprintf('Step 5: 启动基带滑动归一化互相关扫描...\n');
h_matched_bb = flipud(conj(pss_local_bb));
corr_out_bb = filter(h_matched_bb, 1, baseband_rx);
corr_mag_sq_bb = abs(corr_out_bb).^2;

E_local_bb = sum(abs(pss_local_bb).^2);
win_ones_bb = ones(length(pss_local_bb), 1);
E_rx_moving_bb = filter(win_ones_bb, 1, abs(baseband_rx).^2);
E_rx_moving_bb(E_rx_moving_bb < 1e-10) = 1e-10;

corr_norm_pct_bb = (corr_mag_sq_bb ./ (E_local_bb .* E_rx_moving_bb)) * 100;

%% ===== 6. 自动识别所有 PSS 峰值 =====
fprintf('Step 6: 提取所有 PSS 同步峰值...\n');
min_peak_dist_bb = length(pss_local_bb);
[pks_bb, locs_bb] = findpeaks(corr_norm_pct_bb, 'MinPeakHeight', min_peak_height, 'MinPeakDistance', min_peak_dist_bb);

if isempty(locs_bb)
    fprintf('          未检测到高于门限的峰，改用全局最大值兜底。\n');
    [pks_bb, locs_bb] = max(corr_norm_pct_bb);
end

end_pos_global = sync_scan_start + (locs_bb - 1);
start_pos_global = end_pos_global - (L_pss - 1);

fprintf('          检测到 %d 个候选峰。\n', length(pks_bb));
for i = 1:length(pks_bb)
    fprintf('            - Peak %02d: 匹配度 %.2f%%, 末端=%d, 起始=%d\n', ...
        i, pks_bb(i), end_pos_global(i), start_pos_global(i));
end

switch lower(peak_select_mode)
    case 'first'
        peak_indices = 1:min(1, numel(start_pos_global));
    case 'strongest'
        [~, best_idx] = max(pks_bb);
        peak_indices = best_idx;
    case 'all'
        peak_indices = 1:numel(start_pos_global);
    otherwise
        error('未知 peak_select_mode: %s', peak_select_mode);
end

%% ===== 7. 对每个峰执行 SSS 解调，记录四个相位候选 =====
fprintf('Step 7: 对每个峰做 SSS 解调与记录...\n');
all_rows = {};
row = 0;

for k = 1:numel(peak_indices)
    pi_k = peak_indices(k);
    read_start_sample = start_pos_global(pi_k);
    peak_x = start_pos_global(pi_k);
    peak_y = pks_bb(pi_k);

    if read_start_sample < 0
        fprintf('          [跳过] 峰 #%d 起点越界: %d\n', pi_k, read_start_sample);
        continue;
    end

    try
        [hex_candidates, payload_label] = local_demod_four_candidates_baseband( ...
            inFile, read_start_sample, read_length, sss_decode_start_idx, target_offset, ...
            fs_source, fs_target, sro_ppm, delta_f, N_fft, demod_all_subcarriers, ...
            decision_rotate_deg, input_format, dat_header_bytes, plot_constellation_enable, plot_constellation_max_frames, k, peak_y);

        for r = 1:4
            row = row + 1;
            all_rows{row,1} = k;
            all_rows{row,2} = inFile;
            all_rows{row,3} = peak_x;
            all_rows{row,4} = peak_y;
            all_rows{row,5} = read_start_sample;
            all_rows{row,6} = decision_rotate_deg + (r-1)*90;
            all_rows{row,7} = payload_label;
            all_rows{row,8} = hex_candidates{r};
            all_rows{row,9} = length(hex_candidates{r});
        end

        fprintf('          [峰 %02d] peak_x=%8d | read_start=%8d | done\n', k, peak_x, read_start_sample);
    catch ME
        fprintf('          [峰 %02d] error: %s\n', k, ME.message);
    end
end

if row == 0
    error('没有成功解调的样本。');
end

Tall = cell2table(all_rows, 'VariableNames', ...
    {'frame_index','file_name','peak_x','peak_y','read_start_sample','candidate_angle_deg','payload_label','hex_seq','hex_len'});

% 先剔除长度异常（只保留众数长度）
if iscell(Tall.hex_len)
    len_vec = cellfun(@double, Tall.hex_len);
else
    len_vec = double(Tall.hex_len);
end
len_mode = mode(len_vec);
Tall = Tall(len_vec == len_mode, :);

writetable(Tall, out_all_csv);
fprintf('Saved all candidates to %s\n', out_all_csv);

%% ===== 8. 按序列值分组 =====
fprintf('Step 8: 序列聚类分组...\n');
seqs = Tall.hex_seq;
N = height(Tall);

cluster_ids = zeros(N,1);
cluster_centers = {};
cluster_count = 0;

for i = 1:N
    s = seqs{i};
    if cluster_count == 0
        cluster_count = 1;
        cluster_ids(i) = 1;
        cluster_centers{1} = s;
        continue;
    end

    d_best = inf;
    c_best = 0;
    for c = 1:cluster_count
        d = local_hex_dist(s, cluster_centers{c});
        if d < d_best
            d_best = d;
            c_best = c;
        end
    end

    if d_best <= assign_thresh
        cluster_ids(i) = c_best;
    else
        cluster_count = cluster_count + 1;
        cluster_ids(i) = cluster_count;
        cluster_centers{cluster_count} = s;
    end
end

Tall.cluster_id = cluster_ids;

% 取样本数最多的前 4 组
u = unique(cluster_ids);
counts = zeros(size(u));
for k = 1:numel(u)
    counts(k) = sum(cluster_ids == u(k));
end
[~, ord] = sort(counts, 'descend');
keep_clusters = u(ord(1:min(max_groups, numel(ord))));
Tall = Tall(ismember(Tall.cluster_id, keep_clusters), :);

%% ===== 9. 每组离群剔除 + 多数投票 =====
fprintf('Step 9: 离群剔除与多数投票...\n');
group_rows = {};
grow = 0;
unstable_rows = {};
urow = 0;
ber_rows = {};
brow = 0;

final_hex = cell(1, numel(keep_clusters));
final_cnt = zeros(1, numel(keep_clusters));

for gi = 1:numel(keep_clusters)
    cid = keep_clusters(gi);
    idxs = find(Tall.cluster_id == cid);
    seq_group = Tall.hex_seq(idxs);

    proto = local_majority_vote_hex(seq_group);

    keep = false(numel(seq_group),1);
    for j = 1:numel(seq_group)
        keep(j) = local_hex_dist(seq_group{j}, proto) <= outlier_thresh;
    end
    seq_clean = seq_group(keep);
    idxs_clean = idxs(keep);

    if isempty(seq_clean)
        seq_clean = seq_group;
    end

    [voted, vote_stats] = local_majority_vote_with_stats(seq_clean, unstable_ratio_thresh);
    final_hex{gi} = voted;
    final_cnt(gi) = numel(seq_clean);

    unstable_pos = find(vote_stats.winner_ratio < unstable_ratio_thresh | vote_stats.top_tie);
    unstable_count = numel(unstable_pos);

    for uu = 1:unstable_count
        p = unstable_pos(uu);
        urow = urow + 1;
        unstable_rows{urow,1} = cid;
        unstable_rows{urow,2} = p;
        unstable_rows{urow,3} = vote_stats.voted_char(p);
        unstable_rows{urow,4} = vote_stats.winner_ratio(p);
        unstable_rows{urow,5} = vote_stats.top_count(p);
        unstable_rows{urow,6} = vote_stats.second_count(p);
        unstable_rows{urow,7} = vote_stats.total_count(p);
    end

    voted_marked = local_mark_unstable_hex(voted, unstable_pos);

    ber_list = [];
    bits_comp_list = [];
    dropped_full_wrong = 0;
    dropped_no_bits = 0;
    for jj = 1:numel(seq_clean)
        [ber_j, bits_comp_j] = local_hex_ber_excluding_positions(seq_clean{jj}, voted, unstable_pos);
        if bits_comp_j <= 0
            dropped_no_bits = dropped_no_bits + 1;
            continue;
        end
        if ber_j >= full_wrong_ber_thresh
            dropped_full_wrong = dropped_full_wrong + 1;
            continue;
        end
        ber_list(end+1,1) = ber_j; %#ok<AGROW>
        bits_comp_list(end+1,1) = bits_comp_j; %#ok<AGROW>
    end

    if isempty(ber_list)
        avg_ber = NaN;
        used_n = 0;
        avg_bits_comp = 0;
    else
        avg_ber = mean(ber_list);
        used_n = numel(ber_list);
        avg_bits_comp = mean(bits_comp_list);
    end

    brow = brow + 1;
    ber_rows{brow,1} = cid;
    ber_rows{brow,2} = numel(idxs_clean);
    ber_rows{brow,3} = used_n;
    ber_rows{brow,4} = dropped_full_wrong;
    ber_rows{brow,5} = dropped_no_bits;
    ber_rows{brow,6} = avg_ber;
    ber_rows{brow,7} = avg_bits_comp;

    grow = grow + 1;
    group_rows{grow,1} = cid;
    group_rows{grow,2} = numel(seq_group);
    group_rows{grow,3} = numel(seq_clean);
    group_rows{grow,4} = unstable_count;
    group_rows{grow,5} = voted;
    group_rows{grow,6} = voted_marked;
    group_rows{grow,7} = avg_ber;
    group_rows{grow,8} = used_n;
    group_rows{grow,9} = dropped_full_wrong;
end

Tgroup = cell2table(group_rows, 'VariableNames', ...
    {'cluster_id','raw_count','clean_count','unstable_count','voted_hex','voted_hex_marked','avg_ber_excluding_unstable','ber_used_count','ber_dropped_full_wrong'});
writetable(Tgroup, out_group_csv);
fprintf('Saved grouped summary to %s\n', out_group_csv);

if urow > 0
    Tunstable = cell2table(unstable_rows, 'VariableNames', ...
        {'cluster_id','hex_position_1based','voted_char','winner_ratio','top_count','second_count','total_count'});
else
    Tunstable = cell2table(cell(0,7), 'VariableNames', ...
        {'cluster_id','hex_position_1based','voted_char','winner_ratio','top_count','second_count','total_count'});
end
writetable(Tunstable, out_unstable_csv);
fprintf('Saved unstable-position report to %s\n', out_unstable_csv);

Tber = cell2table(ber_rows, 'VariableNames', ...
    {'cluster_id','clean_count','ber_used_count','dropped_full_wrong','dropped_no_bits','avg_ber','avg_compared_bits'});
writetable(Tber, out_ber_csv);
fprintf('Saved BER summary to %s\n', out_ber_csv);

%% ===== 10. 输出最终投票结果 =====
fprintf('\n======= FINAL 4 SEQUENCES (majority vote by cluster) =======\n');
for gi = 1:numel(keep_clusters)
    snip = final_hex{gi}(1:min(48, length(final_hex{gi})));
    ucnt = Tgroup.unstable_count(gi);
    fprintf('Group-%d (cluster=%d, clean_n=%d, unstable=%d, avgBER=%.6f): %s...\n', ...
        gi, keep_clusters(gi), final_cnt(gi), ucnt, Tgroup.avg_ber_excluding_unstable(gi), snip);
end
fprintf('============================================================\n');

%% ===== 本地函数 =====
function pss_local_bb = local_build_baseband_pss(fs_target)
    fs_symbol = 60e6;
    hex_str = 'BD3CD0148871751F84CED8C1BE32AC96';
    hex_chars = char(hex_str);
    bits = zeros(1, 128);
    for i = 1:length(hex_chars)
        bin_str = dec2bin(hex2dec(hex_chars(i)), 4);
        bits((i-1)*4 + 1 : i*4) = bin_str - '0';
    end

    d_phi = (bits == 1) * (-pi/2) + (bits == 0) * (pi/2);
    phase_accum = cumsum(d_phi);
    m_syms = exp(1j * phase_accum);

    pss_base = [-m_syms, m_syms, m_syms, m_syms, m_syms, m_syms, m_syms, m_syms];
    N_cp_syms = 48;
    cp_syms = -pss_base(end - N_cp_syms + 1 : end);
    pss_base_with_cp = [cp_syms, pss_base];

    [P, Q] = rat(fs_target / fs_symbol);
    pss_up = resample(pss_base_with_cp, P, Q);
    pss_local_bb = pss_up / max(abs(pss_up));
    pss_local_bb = pss_local_bb(:);
end

function delta_f = local_blind_cfo_estimate(x_for_cfo, fs_source, M, Nfft_cfo, search_range, target_f)
    x_for_cfo = x_for_cfo(:);
    if length(x_for_cfo) < 4096
        delta_f = 0;
        return;
    end

    Z_spec = fftshift(fft(x_for_cfo.^M, Nfft_cfo));
    f_axis_cfo = (-Nfft_cfo/2 : Nfft_cfo/2 - 1) / Nfft_cfo * fs_source;
    target_f_mapped = mod(target_f + fs_source/2, fs_source) - fs_source/2;
    [~, center_idx] = min(abs(f_axis_cfo - target_f_mapped));
    local_lo = max(1, center_idx - search_range);
    local_hi = min(Nfft_cfo, center_idx + search_range);
    local_spec = abs(Z_spec(local_lo:local_hi));
    [~, local_max_idx] = max(local_spec);
    actual_max_idx = local_lo + local_max_idx - 1;
    f_measured = f_axis_cfo(actual_max_idx);
    delta_f = (f_measured - target_f_mapped) / M;
end

function [hex_candidates, payload_label] = local_demod_four_candidates_baseband( ...
    inFile, read_start_sample, read_length, sss_decode_start_idx, target_offset, ...
    fs_source, fs_target, sro_ppm, cfo_hz, N_fft, demod_all_subcarriers, ...
    decision_rotate_deg, input_format, dat_header_bytes, plot_constellation_enable, plot_constellation_max_frames, frame_index, peak_y)

    [x_raw, ~] = read_iq_auto_local(inFile, read_start_sample, read_length, input_format, dat_header_bytes);
    x_raw = double(x_raw);
    x_raw = x_raw - mean(x_raw);
    x_raw = x_raw / (mean(abs(x_raw)) + eps);

    t_vec = (0:length(x_raw)-1).' / fs_source;
    x_cfo = x_raw .* exp(-1j * 2 * pi * cfo_hz * t_vec);

    fs_eff = fs_source * (1 + sro_ppm/1e6);
    Wn = 35e6 / (fs_source / 2);
    b_lpf = fir1(30, Wn);
    x_cfo_filtered = filtfilt(b_lpf, 1, x_cfo);

    T_in = 1 / fs_eff;
    T_out = 1 / fs_target;
    t_out = 0 : T_out : (length(x_cfo_filtered)-3) * T_in;

    idx_frac = t_out / T_in + 1;
    idx_base = floor(idx_frac);
    mu = idx_frac - idx_base;

    valid_mask = (idx_base >= 2) & (idx_base <= length(x_cfo_filtered)-2);
    idx_base = idx_base(valid_mask);
    mu = mu(valid_mask);

    h0 = -(mu - 1) .* (mu - 2) .* mu / 6;
    h1 =  (mu - 1) .* (mu + 1) .* (mu - 2) / 2;
    h2 = -(mu + 1) .* mu .* (mu - 2) / 2;
    h3 =  (mu + 1) .* (mu - 1) .* mu / 6;

    x_sro = h0 .* x_cfo_filtered(idx_base - 1) + ...
            h1 .* x_cfo_filtered(idx_base) + ...
            h2 .* x_cfo_filtered(idx_base + 1) + ...
            h3 .* x_cfo_filtered(idx_base + 2);
    x_sro = x_sro(:);

    sss_start_idx = sss_decode_start_idx + target_offset;
    if sss_start_idx < 1 || sss_start_idx + N_fft > length(x_sro)
        error('SSS 解调窗口越界。');
    end

    x_sss_time = x_sro(sss_start_idx : sss_start_idx + N_fft - 1);
    x_sss_freq = fft(x_sss_time, N_fft) / sqrt(N_fft);

    pwr_bins = abs(x_sss_freq).^2;
    valid_mask_power = pwr_bins > mean(pwr_bins) * 0.1;
    dc_guard = 3;
    valid_mask_power(1:dc_guard) = false;
    valid_mask_power(N_fft-dc_guard+1:N_fft) = false;
    valid_idxs = find(valid_mask_power);

    rel_freq = zeros(N_fft,1);
    rel_freq(1:N_fft/2) = (0:N_fft/2 - 1).';
    rel_freq(N_fft/2+1:N_fft) = (-N_fft/2:-1).';

    freq_indices = rel_freq(valid_idxs);
    syms_valid = x_sss_freq(valid_idxs);
    [freq_indices, sort_idx] = sort(freq_indices);
    syms_valid = syms_valid(sort_idx);
    valid_idxs = valid_idxs(sort_idx);

    syms_pow4 = syms_valid.^4;
    df = diff(freq_indices);
    seg_start = [1; find(df > 1) + 1];
    seg_end   = [find(df > 1); length(freq_indices)];

    slope_list = [];
    w_list = [];
    for k = 1:length(seg_start)
        idx_seg = seg_start(k):seg_end(k);
        if numel(idx_seg) < 8
            continue;
        end
        f_seg = freq_indices(idx_seg);
        ph_seg = unwrap(angle(syms_pow4(idx_seg)));
        p_seg = polyfit(f_seg, ph_seg, 1);
        slope_list(end+1,1) = p_seg(1); %#ok<AGROW>
        w_list(end+1,1) = numel(idx_seg); %#ok<AGROW>
    end

    if isempty(slope_list)
        error('相位拟合失败。');
    end

    slope4 = sum(slope_list .* w_list) / sum(w_list);
    phase0_4 = angle(mean(syms_pow4 .* exp(-1j * slope4 * freq_indices)));

    full_freq_indices = [(0:N_fft/2-1), (-N_fft/2:-1)].';
    phase_correction = (slope4 * full_freq_indices + phase0_4) / 4;
    x_sss_freq_corr = x_sss_freq .* exp(-1j * phase_correction);

    if plot_constellation_enable && frame_index <= plot_constellation_max_frames
        if demod_all_subcarriers
            syms_for_plot = x_sss_freq_corr(:);
            plot_label = sprintf('All carriers: %d', N_fft);
        else
            syms_for_plot = x_sss_freq_corr(valid_idxs);
            plot_label = sprintf('Valid carriers: %d', length(valid_idxs));
        end

        figure('Name', sprintf('Frame %d SSS Constellation', frame_index), 'Position', [120, 120, 1100, 520]);
        subplot(1, 2, 1);
        plot(real(syms_for_plot), imag(syms_for_plot), 'b.', 'MarkerSize', 8);
        grid on; axis square;
        xlim([-2 2]); ylim([-2 2]);
        xlabel('I'); ylabel('Q');
        title(sprintf('Frame %d: 相位补偿后星座\n峰值 %.2f%% | %s', frame_index, peak_y, plot_label));

        subplot(1, 2, 2);
        syms_rot_demo = syms_for_plot .* exp(1j * decision_rotate_deg * pi/180);
        plot(real(syms_rot_demo), imag(syms_rot_demo), 'r.', 'MarkerSize', 8);
        grid on; axis square;
        xlim([-2 2]); ylim([-2 2]);
        xlabel('I'); ylabel('Q');
        title(sprintf('Frame %d: 再旋转 %d° 后判决视图', frame_index, decision_rotate_deg));
    end

    if demod_all_subcarriers
        syms_payload = x_sss_freq_corr(:);
        payload_label = sprintf('All carriers: %d', N_fft);
    else
        syms_payload = x_sss_freq_corr(valid_idxs);
        payload_label = sprintf('Valid carriers: %d', length(valid_idxs));
    end

    hex_candidates = cell(1,4);
    for r = 1:4
        syms_payload_rot = syms_payload .* exp(1j * decision_rotate_deg * pi/180) .* exp(-1j * (r-1) * pi/2);
        bits_I = real(syms_payload_rot) < 0;
        bits_Q = imag(syms_payload_rot) < 0;

        demod_bits = zeros(length(syms_payload_rot)*2, 1);
        demod_bits(1:2:end) = bits_I;
        demod_bits(2:2:end) = bits_Q;

        hex_candidates{r} = local_bits_to_hex(demod_bits);
    end
end

function hex_str = local_bits_to_hex(bits)
    hex_str = '';
    for i = 1:4:length(bits)
        if i+3 > length(bits)
            chunk = bits(i:end);
            val = 0;
            for b = 1:length(chunk)
                val = val + chunk(b) * 2^(length(chunk)-b);
            end
            hex_str = [hex_str, dec2hex(val)]; %#ok<AGROW>
            break;
        end
        c = bits(i:i+3);
        val = c(1)*8 + c(2)*4 + c(3)*2 + c(4);
        hex_str = [hex_str, dec2hex(val)]; %#ok<AGROW>
    end
end

function [x_out, meta_out] = read_iq_auto_local(inFile, st_sample, num_samples, format_mode, dat_hd_bytes)
    if strcmpi(format_mode, 'auto')
        if endsWith(lower(inFile), '.dat')
            fmt = 'dat';
        else
            fmt = 'iq';
        end
    else
        fmt = format_mode;
    end

    meta_out = struct('numSamplesRead', 0);
    x_out = [];

    fid = fopen(inFile, 'rb');
    if fid == -1
        return;
    end

    if strcmpi(fmt, 'dat')
        if dat_hd_bytes > 0
            fseek(fid, dat_hd_bytes, 'bof');
        end
        bytes_per_cplx = 4;
        offset_bytes = st_sample * bytes_per_cplx;
        fseek(fid, offset_bytes, 'cof');

        N_read = num_samples * 2;
        data = fread(fid, N_read, 'int16=>double');
        fclose(fid);

        if isempty(data)
            return;
        end
        if mod(length(data), 2) == 1
            data = [data; 0];
        end
        I = data(1:2:end);
        Q = data(2:2:end);
        x_out = I + 1j * Q;
        meta_out.numSamplesRead = length(x_out);
    else
        fclose(fid);
        [x_cplx, meta] = iq_read_int16_le(inFile, st_sample, num_samples);
        x_out = double(x_cplx);
        meta_out = meta;
    end
end

function d = local_hex_dist(a, b)
    L = min(length(a), length(b));
    if L == 0
        d = 1;
        return;
    end
    d = sum(a(1:L) ~= b(1:L)) / L;
end

function voted = local_majority_vote_hex(seq_cell)
    if isempty(seq_cell)
        voted = '';
        return;
    end
    L = min(cellfun(@length, seq_cell));
    voted_chars = repmat('0', 1, L);
    hexset = ['0':'9','A':'F'];
    for p = 1:L
        col = zeros(1, numel(seq_cell));
        for k = 1:numel(seq_cell)
            col(k) = seq_cell{k}(p);
        end
        cnt = zeros(size(hexset));
        for h = 1:numel(hexset)
            cnt(h) = sum(col == hexset(h));
        end
        [~, m] = max(cnt);
        voted_chars(p) = hexset(m);
    end
    voted = voted_chars;
end

function [voted, stats] = local_majority_vote_with_stats(seq_cell, unstable_ratio_thresh)
    if isempty(seq_cell)
        voted = '';
        stats = struct('winner_ratio', [], 'top_tie', [], 'top_count', [], 'second_count', [], 'total_count', [], 'voted_char', '');
        return;
    end

    L = min(cellfun(@length, seq_cell));
    voted_chars = repmat('0', 1, L);
    winner_ratio = zeros(1, L);
    top_tie = false(1, L);
    top_count = zeros(1, L);
    second_count = zeros(1, L);
    total_count = zeros(1, L);

    hexset = ['0':'9','A':'F'];
    for p = 1:L
        col = zeros(1, numel(seq_cell));
        for k = 1:numel(seq_cell)
            col(k) = seq_cell{k}(p);
        end

        cnt = zeros(size(hexset));
        for h = 1:numel(hexset)
            cnt(h) = sum(col == hexset(h));
        end

        [cnt_sorted, idx_sorted] = sort(cnt, 'descend');
        voted_chars(p) = hexset(idx_sorted(1));
        top_count(p) = cnt_sorted(1);
        if numel(cnt_sorted) >= 2
            second_count(p) = cnt_sorted(2);
        end
        total_count(p) = numel(seq_cell);
        winner_ratio(p) = cnt_sorted(1) / max(1, total_count(p));
        top_tie(p) = (cnt_sorted(1) == cnt_sorted(2));
    end

    voted = voted_chars;
    stats = struct();
    stats.winner_ratio = winner_ratio;
    stats.top_tie = top_tie;
    stats.top_count = top_count;
    stats.second_count = second_count;
    stats.total_count = total_count;
    stats.voted_char = voted_chars;
    stats.unstable_ratio_thresh = unstable_ratio_thresh;
end

function out = local_mark_unstable_hex(hex_str, unstable_pos)
    if isempty(hex_str)
        out = '';
        return;
    end

    L = length(hex_str);
    mask = false(1, L);
    unstable_pos = unique(unstable_pos(:).');
    unstable_pos = unstable_pos(unstable_pos >= 1 & unstable_pos <= L);
    mask(unstable_pos) = true;

    out = '';
    for i = 1:L
        if mask(i)
            out = [out '"' hex_str(i) '"']; %#ok<AGROW>
        else
            out = [out hex_str(i)]; %#ok<AGROW>
        end
    end
end

function [ber, bits_comp] = local_hex_ber_excluding_positions(seq_hex, ref_hex, exclude_pos)
    L = min(length(seq_hex), length(ref_hex));
    if L <= 0
        ber = NaN;
        bits_comp = 0;
        return;
    end

    mask = true(1, L);
    exclude_pos = unique(exclude_pos(:).');
    exclude_pos = exclude_pos(exclude_pos >= 1 & exclude_pos <= L);
    mask(exclude_pos) = false;

    idxs = find(mask);
    if isempty(idxs)
        ber = NaN;
        bits_comp = 0;
        return;
    end

    bit_errors = 0;
    for ii = 1:numel(idxs)
        p = idxs(ii);
        va = local_hex_char_to_val(seq_hex(p));
        vb = local_hex_char_to_val(ref_hex(p));
        bit_errors = bit_errors + local_popcount4(bitxor(va, vb));
    end

    bits_comp = 4 * numel(idxs);
    ber = bit_errors / bits_comp;
end

function v = local_hex_char_to_val(ch)
    if ch >= '0' && ch <= '9'
        v = double(ch) - double('0');
    elseif ch >= 'A' && ch <= 'F'
        v = 10 + double(ch) - double('A');
    elseif ch >= 'a' && ch <= 'f'
        v = 10 + double(ch) - double('a');
    else
        v = 0;
    end
end

function c = local_popcount4(v)
    lut = [0 1 1 2 1 2 2 3 1 2 2 3 2 3 3 4];
    c = lut(v + 1);
end
