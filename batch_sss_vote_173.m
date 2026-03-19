% batch_sss_vote_173.m
% 目标：
% 1) 遍历 sigtest1~sigtest173
% 2) 先检测滑动自相关峰值 x 坐标
% 3) read_start_sample = peak_x - 874
% 4) 对每个文件做 SSS 解调，得到 4 个旋转候选序列
% 5) 丢弃明显错误样本（解调失败/长度异常/离群）
% 6) 按序列相似性分成 4 组
% 7) 每组按位多数投票，输出 4 个最终序列

clear; clc;

%% ===== 参数区（按需改） =====
file_prefix = 'sigtest';
file_suffix = '.iq';
file_range = 1:173;

% 峰值检测参数（与 autocorr_test2_slide128.m 保持一致）
peak_start_sample0 = 0;
peak_num_samples = 50000;
peak_W = 874;
peak_D = 874;
peak_min = 0.5;

% 解调参数（与 sss_demodulation.m 对齐）
read_length = 6992*3 + 1000;
sss_decode_start_idx = 1048;
target_offset = 10;

fs_source = 409.6e6;
fs_target = 60e6;
sro_ppm  = 0;
cfo_hz   = 17556;
N_fft = 1024;
freq_shift_hz = 63e6;
demod_all_subcarriers = true;
decision_rotate_deg = 45;

% 分组/剔除参数
max_groups = 4;
assign_thresh = 0.18;   % 分配到簇的最大归一化汉明距离
outlier_thresh = 0.22;  % 相对簇中心超过该距离视为离群
unstable_ratio_thresh = 0.60; % 多数票占比低于该阈值，判为“不稳定位”
full_wrong_ber_thresh = 0.40; % BER 高于该阈值视为“全错样本”并剔除

% 输出
out_all_csv = 'sss_all_candidates_173.csv';
out_group_csv = 'sss_group_summary_173.csv';
out_unstable_csv = 'sss_unstable_positions_173.csv';
out_ber_csv = 'sss_group_ber_173.csv';

%% ===== 主流程 =====
all_rows = {};
row = 0;

fprintf('Start batch SSS voting over 173 files...\n');
for idx = file_range
    inFile = sprintf('%s%d%s', file_prefix, idx, file_suffix);
    if ~exist(inFile, 'file')
        fprintf('[%3d] %-16s | missing\n', idx, inFile);
        continue;
    end

    try
        % 1) 检测峰值 x 坐标
        [peak_x, peak_y] = local_detect_peak_x(inFile, peak_start_sample0, peak_num_samples, peak_W, peak_D, peak_min);
        read_start_sample = peak_x - peak_D;

        if read_start_sample < 0
            fprintf('[%3d] %-16s | invalid read_start=%d (skip)\n', idx, inFile, read_start_sample);
            continue;
        end

        % 2) 用该 read_start_sample 解调，取四个候选序列
        [hex_candidates, payload_label] = local_demod_four_candidates( ...
            inFile, read_start_sample, read_length, sss_decode_start_idx, target_offset, ...
            fs_source, fs_target, sro_ppm, cfo_hz, N_fft, freq_shift_hz, ...
            demod_all_subcarriers, decision_rotate_deg);

        for r = 1:4
            row = row + 1;
            all_rows{row,1} = idx;
            all_rows{row,2} = inFile;
            all_rows{row,3} = peak_x;
            all_rows{row,4} = peak_y;
            all_rows{row,5} = read_start_sample;
            all_rows{row,6} = decision_rotate_deg + (r-1)*90; % 候选角度
            all_rows{row,7} = payload_label;
            all_rows{row,8} = hex_candidates{r};
            all_rows{row,9} = length(hex_candidates{r});
        end

        fprintf('[%3d] %-16s | peak_x=%8d | read_start=%8d | done\n', idx, inFile, peak_x, read_start_sample);

    catch ME
        fprintf('[%3d] %-16s | error: %s\n', idx, inFile, ME.message);
    end
end

if row == 0
    error('没有成功解调的样本。');
end

Tall = cell2table(all_rows, 'VariableNames', ...
    {'file_index','file_name','peak_x','peak_y','read_start_sample','candidate_angle_deg','payload_label','hex_seq','hex_len'});

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

%% 3) 按序列值分组（目标四组）
seqs = Tall.hex_seq;
N = height(Tall);

cluster_ids = zeros(N,1);
cluster_centers = {}; % 每簇临时中心（字符串）
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

% 取样本数最多的前4组
u = unique(cluster_ids);
counts = zeros(size(u));
for k = 1:numel(u)
    counts(k) = sum(cluster_ids == u(k));
end
[~, ord] = sort(counts, 'descend');
keep_clusters = u(ord(1:min(max_groups, numel(ord))));

Tall = Tall(ismember(Tall.cluster_id, keep_clusters), :);

%% 4) 每组离群剔除 + 多数投票
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

    % 先做一次初始投票得到中心
    proto = local_majority_vote_hex(seq_group);

    % 按与 proto 的距离剔除明显错误结果
    keep = false(numel(seq_group),1);
    for j = 1:numel(seq_group)
        keep(j) = local_hex_dist(seq_group{j}, proto) <= outlier_thresh;
    end
    seq_clean = seq_group(keep);
    idxs_clean = idxs(keep);

    if isempty(seq_clean)
        seq_clean = seq_group; % 防止全被剔除
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

    % === 以投票序列为基准计算平均 BER（剔除不稳定位 + 全错样本） ===
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

%% 5) 输出四个最终序列值
fprintf('\n======= FINAL 4 SEQUENCES (majority vote by cluster) =======\n');
for gi = 1:numel(keep_clusters)
    snip = final_hex{gi}(1:min(48, length(final_hex{gi})));
    ucnt = Tgroup.unstable_count(gi);
    fprintf('Group-%d (cluster=%d, clean_n=%d, unstable=%d, avgBER=%.6f): %s...\n', ...
        gi, keep_clusters(gi), final_cnt(gi), ucnt, Tgroup.avg_ber_excluding_unstable(gi), snip);
end
fprintf('============================================================\n');

%% ===== 本地函数 =====
function [peak_x, peak_y] = local_detect_peak_x(inFile, startSample0, numSamples, W, D, peak_min)
    [x, meta] = iq_read_int16_le(inFile, startSample0, local_num_samples(inFile, numSamples, startSample0));
    L = meta.numSamplesRead;
    if L <= D || L < W
        error('peak detect: data too short (L=%d).', L);
    end

    x = double(x) / 32768;
    x = x - mean(x);

    rx_delayed = x(1+D:end);
    rx_base    = x(1:end-D);

    P_metric = filter(ones(1, W), 1, conj(rx_base) .* rx_delayed);
    R_base = filter(ones(1, W), 1, abs(rx_base).^2);
    R_delayed = filter(ones(1, W), 1, abs(rx_delayed).^2);

    M_real = real(P_metric ./ (sqrt(R_base .* R_delayed) + 1e-10));
    M_real_full = [M_real; zeros(D,1)];
    t_axis = (0:length(x)-1).' + startSample0;

    [peak_x, peak_y] = local_pick_peak(t_axis, M_real_full, peak_min);
end

function [hex_candidates, payload_label] = local_demod_four_candidates( ...
    inFile, read_start_sample, read_length, sss_decode_start_idx, target_offset, ...
    fs_source, fs_target, sro_ppm, cfo_hz, N_fft, freq_shift_hz, ...
    demod_all_subcarriers, decision_rotate_deg)

    [x_raw, ~] = iq_read_int16_le(inFile, read_start_sample, read_length);
    x_raw = double(x_raw);
    x_raw = x_raw - mean(x_raw);
    x_raw = x_raw / mean(abs(x_raw));

    t_vec = (0:length(x_raw)-1) / fs_source;
    x_shifted = x_raw .* exp(-1j * 2 * pi * freq_shift_hz * t_vec).';
    x_cfo = x_shifted .* exp(-1j * 2 * pi * cfo_hz * t_vec).';

    fs_eff = fs_source * (1 + sro_ppm/1e6);
    Wn = 35e6 / (fs_source / 2);
    b_lpf = fir1(30, Wn);
    x_cfo_filtered = filtfilt(b_lpf, 1, x_cfo);

    T_in = 1 / fs_eff;
    T_out = 1 / fs_target;
    t_out = 0 : T_out : (length(x_cfo_filtered)-3)*T_in;

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

    sss_start_idx_60 = sss_decode_start_idx + target_offset;
    if sss_start_idx_60 < 1 || sss_start_idx_60 + N_fft > length(x_sro)
        error('sss window out of range.');
    end

    x_sss_time = x_sro(sss_start_idx_60 : sss_start_idx_60 + N_fft - 1);
    x_sss_freq = fft(x_sss_time, N_fft) / sqrt(N_fft);

    pwr_bins = abs(x_sss_freq).^2;
    valid_mask_power = pwr_bins > mean(pwr_bins)*0.1;
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

    slope_list = []; w_list = [];
    for k = 1:length(seg_start)
        idx_seg = seg_start(k):seg_end(k);
        if numel(idx_seg) < 8, continue; end
        f_seg = freq_indices(idx_seg);
        ph_seg = unwrap(angle(syms_pow4(idx_seg)));
        p_seg = polyfit(f_seg, ph_seg, 1);
        slope_list(end+1,1) = p_seg(1); %#ok<AGROW>
        w_list(end+1,1) = numel(idx_seg); %#ok<AGROW>
    end
    if isempty(slope_list)
        error('phase fit failed.');
    end

    slope4 = sum(slope_list .* w_list) / sum(w_list);
    phase0_4 = angle(mean(syms_pow4 .* exp(-1j * slope4 * freq_indices)));

    full_freq_indices = [(0:N_fft/2-1), (-N_fft/2:-1)].';
    phase_correction = (slope4 * full_freq_indices + phase0_4) / 4;
    x_sss_freq_corr = x_sss_freq .* exp(-1j * phase_correction);

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

        demod_bits = zeros(length(syms_payload_rot)*2,1);
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

function [peak_x, peak_y] = local_pick_peak(t_axis, y, peak_min)
    n = length(y);
    if n < 3
        [peak_y, idx] = max(y);
        peak_x = t_axis(idx);
        return;
    end
    cand = find(y(2:end-1) > y(1:end-2) & y(2:end-1) >= y(3:end)) + 1;
    cand = cand(y(cand) >= peak_min);
    if isempty(cand)
        [peak_y, idx] = max(y);
    else
        [peak_y, iBest] = max(y(cand));
        idx = cand(iBest);
    end
    peak_x = t_axis(idx);
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
    for p = 1:L
        col = zeros(1, numel(seq_cell));
        for k = 1:numel(seq_cell)
            col(k) = seq_cell{k}(p);
        end
        % 在十六进制字符集合内投票
        hexset = ['0':'9','A':'F'];
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

    % 保留阈值字段，便于外部调试
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

function n = local_num_samples(inFile, numSamples, startSample0)
    if ~isempty(numSamples)
        n = numSamples;
        return;
    end
    info = dir(inFile);
    if isempty(info)
        error('找不到文件: %s', inFile);
    end
    Ntotal = floor(info.bytes / 4);
    if startSample0 >= Ntotal
        n = 0;
    else
        n = Ntotal - startSample0;
    end
end
