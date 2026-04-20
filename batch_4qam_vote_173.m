% batch_4qam_vote_173.m
% 针对173个文件，仅提取4QAM导频，使用固定 offset=0 进行盲相差补偿，
% 提取出四种旋转相位下的比特序列，并通过一致性投票（2-Pass）找出固定不变的比特位。

clear; clc; close all;

%% 1. 参数设置
num_files = 173;
file_prefix = 'sigtest';

% 固定同步与抽取参数
sss_decode_start_idx = 1024+48+1024+48;   
fs_source = 409.6e6;
fs_target = 60e6;
N_fft = 1024;
sro_ppm = 0;

% 我们仅关注预设的4QAM导频子载波
pilot_indices = [488:495, 528:535]; 
num_pilot = length(pilot_indices);
bits_per_sym = 2;
total_bits = num_pilot * bits_per_sym; % 16 * 2 = 32 bits

% 存储所有文件的 4 个相位假设数据
% 维度：[173个文件, 4种假设, 32个比特]
all_files_bits = NaN(num_files, 4, total_bits);
valid_file_mask = false(num_files, 1);

%% 2. 遍历所有文件处理
for file_idx = 1:num_files
    inFile = sprintf('%s%d.iq', file_prefix, file_idx);
    if ~exist(inFile, 'file')
        fprintf('文件 %s 不存在，跳过。\n', inFile);
        continue;
    end
    
    fprintf('\n=== 正在处理文件 %d/173: %s ===\n', file_idx, inFile);
    
    %% --- A. 读取与盲同步 ---
    sync_scan_start = 0;
    sync_scan_len   = 5e6;
    [x_scan, meta_scan] = iq_read_int16_le(inFile, sync_scan_start, sync_scan_len);
    if meta_scan.numSamplesRead <= 0
        warning('扫描数据读取失败'); continue;
    end
    x_scan = double(x_scan(:));
    x_scan = x_scan - mean(x_scan);
    x_scan = x_scan / (mean(abs(x_scan)) + eps);

    [est_f_center, ~] = estimate_center_pilot_blind_from_signal( ...
        x_scan, fs_source, 60e6, 0.35, 2.0e6, 1.2e6, 12.0, 20e6, ...
        6e6, 0.9e6, 4.5, 0.25e6, false, 0.6e6, 2, 256, 120e3);
    freq_shift_hz = est_f_center; 
    cfo_hz = 0;                   

    fs_symbol = fs_target;
    fc_sync = 63e6; 
    hex_str = 'BD3CD0148871751F84CED8C1BE32AC96';
    hex_chars = char(hex_str);
    bits_dummy = zeros(1, 128);
    for i = 1:length(hex_chars)
        bits_dummy((i-1)*4 + 1 : i*4) = dec2bin(hex2dec(hex_chars(i)), 4) - '0';
    end
    d_phi = (bits_dummy == 1) * (-pi/2) + (bits_dummy == 0) * (pi/2);
    m_syms = exp(1j * cumsum(d_phi));
    pss_base = [-m_syms, m_syms, m_syms, m_syms, m_syms, m_syms, m_syms, m_syms];
    pss_base_with_cp = [-pss_base(end - 47 : end), pss_base];

    [P_res, Q_res] = rat(fs_source / fs_symbol);
    pss_up = resample(pss_base_with_cp, P_res, Q_res);
    t_local = (0:length(pss_up)-1) / fs_source;
    pss_local = pss_up .* exp(1j * 2 * pi * fc_sync * t_local);
    pss_local = pss_local / max(abs(pss_local)); 

    h_matched = fliplr(conj(pss_local));
    corr_out = filter(h_matched, 1, x_scan);
    corr_mag_sq = abs(corr_out).^2;

    E_local = sum(abs(pss_local).^2);
    E_rx_moving = filter(ones(1, length(pss_local)), 1, abs(x_scan).^2);
    E_rx_moving(E_rx_moving < 1e-10) = 1e-10;

    corr_norm_pct = (corr_mag_sq ./ (E_local .* E_rx_moving)) * 100;
    [~, locs_all] = findpeaks(corr_norm_pct, 'MinPeakHeight', 50, 'MinPeakDistance', length(pss_local));
    if isempty(locs_all)
        [~, max_idx] = max(corr_norm_pct);
        locs_all = max_idx;
    end
    [~, max_pk_idx] = max(corr_norm_pct(locs_all));
    read_start_sample = locs_all(max_pk_idx) - length(pss_local);
    
    %% --- B. 信号读取与 DDC/SRO 预处理 ---
    read_length = 6992*3+1000;
    [x_raw, ~] = iq_read_int16_le(inFile, read_start_sample, read_length);
    x_raw = double(x_raw - mean(x_raw));
    x_raw = x_raw / mean(abs(x_raw));

    t_vec = (0:length(x_raw)-1) / fs_source;
    x_shifted = x_raw .* exp(-1j * 2 * pi * freq_shift_hz * t_vec).';
    x_cfo = x_shifted .* exp(-1j * 2 * pi * cfo_hz * t_vec).';

    fs_eff = fs_source * (1 + sro_ppm/1e6);
    b_lpf = fir1(30, 35e6 / (fs_source / 2));
    x_cfo_filtered = filtfilt(b_lpf, 1, x_cfo);

    T_in = 1 / fs_eff; T_out = 1 / fs_target;
    t_out = 0 : T_out : (length(x_cfo_filtered)-3)*T_in;
    idx_frac = t_out / T_in + 1; idx_base = floor(idx_frac); mu = idx_frac - idx_base;
    valid_mask = (idx_base >= 2) & (idx_base <= length(x_cfo_filtered)-2);
    idx_base = idx_base(valid_mask); mu = mu(valid_mask);

    h0 = -(mu - 1) .* (mu - 2) .* mu / 6;
    h1 =  (mu - 1) .* (mu + 1) .* (mu - 2) / 2;
    h2 = -(mu + 1) .* mu .* (mu - 2) / 2;
    h3 =  (mu + 1) .* (mu - 1) .* mu / 6;
    x_sro = h0.*x_cfo_filtered(idx_base-1) + h1.*x_cfo_filtered(idx_base) + ...
            h2.*x_cfo_filtered(idx_base+1) + h3.*x_cfo_filtered(idx_base+2);
    x_sro = x_sro(:);

    %% --- C. 提取 4QAM 数据进行解析 (Offset=0) ---
    offset = 0; % 固定使用offset=0
    sss_start_idx_60 = sss_decode_start_idx + offset;
    if sss_start_idx_60 + N_fft > length(x_sro)
        warning('下标越界，跳过'); continue;
    end

    x_sss_time = x_sro(sss_start_idx_60 : sss_start_idx_60 + N_fft - 1);
    x_sss_freq = fft(x_sss_time, N_fft) / sqrt(N_fft);

    rel_freq = zeros(N_fft, 1);
    rel_freq(1:N_fft/2) = (0:N_fft/2 - 1).';
    rel_freq(N_fft/2+1:N_fft) = (-N_fft/2:-1).';

    freq_indices_pilot = rel_freq(pilot_indices);
    syms_valid_pilot = x_sss_freq(pilot_indices);

    [freq_indices_pilot, sort_idx] = sort(freq_indices_pilot);
    syms_valid_pilot = syms_valid_pilot(sort_idx);

    full_freq_indices = [(0:N_fft/2-1), (-N_fft/2:-1)].';

    % 基于4QAM导频计算最优相偏
    [x_sss_freq_corr, slope_hat, phi_hat, ~] = blind_phase_slope_comp_4qam( ...
        x_sss_freq, syms_valid_pilot, freq_indices_pilot, full_freq_indices);

    % 取出最终校正后的重点4QAM星座点 (已排序过，所以直接重取对应位置即可)
    syms_4qam_corrected = x_sss_freq_corr(pilot_indices);
    % 为了输出比特顺序稳定，我们也按频点从小到大排一下
    syms_4qam_corrected = syms_4qam_corrected(sort_idx); 

    %% --- D. 4种相位假设并映射为比特 ---
    phase_shifts = [0, pi/2, pi, 3*pi/2];
    for p_idx = 1:4
        rotated_syms = syms_4qam_corrected * exp(1j * phase_shifts(p_idx));
        
        % QPSK/4QAM 硬判决到 bit
        % 映射规则: 实部 I>0 -> 0, I<0 -> 1; 虚部 Q>0 -> 0, Q<0 -> 1 (格雷码常规)
        bits_hyp = zeros(1, total_bits);
        bits_hyp(1:2:end) = real(rotated_syms) < 0;
        bits_hyp(2:2:end) = imag(rotated_syms) < 0;
        
        all_files_bits(file_idx, p_idx, :) = bits_hyp;
    end
    
    valid_file_mask(file_idx) = true;
end

%% 3. 比特一致性投票与分析 (2-Pass 算法复原与投票 4 序列)
fprintf('\n=======================================\n');
fprintf('所有文件解码完毕，开始进行跨文件 2-Pass 相位对齐与投票。\n');

valid_indices = find(valid_file_mask);
if isempty(valid_indices)
    error('没有有效的文件数据可供分析。');
end
num_valid = length(valid_indices);

% --- Pass 1: 初步对齐 ---
% 选定第一个文件的 Phase 1 为绝对基准(Anchor)
anchor_bits = squeeze(all_files_bits(valid_indices(1), 1, :)).'; 

pass1_aligned_bits = zeros(num_valid, total_bits);
for i = 1:num_valid
    f_idx = valid_indices(i);
    min_dist = inf;
    best_h = 1;
    for hyp = 1:4
        test_bits = squeeze(all_files_bits(f_idx, hyp, :)).';
        dist = sum(test_bits ~= anchor_bits);
        if dist < min_dist
            min_dist = dist;
            best_h = hyp;
        end
    end
    pass1_aligned_bits(i, :) = squeeze(all_files_bits(f_idx, best_h, :)).';
end

% 寻找初步高一致性位 (> 80%) 作 Mask
bit_agreement_ratio = sum(pass1_aligned_bits == mode(pass1_aligned_bits, 1), 1) / num_valid;
stable_mask = bit_agreement_ratio > 0.80;
if sum(stable_mask) == 0
    stable_mask = true(1, total_bits); 
end

% --- Pass 2: 利用 stable_mask 精确对齐 ---
anchor_stable_bits = anchor_bits(stable_mask);
final_aligned_bits = zeros(num_valid, total_bits);

for i = 1:num_valid
    f_idx = valid_indices(i);
    min_dist = inf;
    best_h = 1;
    for hyp = 1:4
        test_bits = squeeze(all_files_bits(f_idx, hyp, :)).';
        dist = sum(test_bits(stable_mask) ~= anchor_stable_bits);
        if dist < min_dist
            min_dist = dist;
            best_h = hyp;
        end
    end
    final_aligned_bits(i, :) = squeeze(all_files_bits(f_idx, best_h, :)).';
end

% 多数投票 (不管一不一样，暴力对齐后取众数)，得到唯一合成序列
final_mode_bits = mode(final_aligned_bits, 1);

% --- 提取最终序列的另外 3 个相位模糊变体 ---
% 1. 逆映射回 QPSK 星座符号 (映射规则: bits=0 -> I/Q>0即+1; bits=1 -> I/Q<0即-1)
syms_I = 1 - 2 * final_mode_bits(1:2:end);
syms_Q = 1 - 2 * final_mode_bits(2:2:end);
base_syms = syms_I + 1j * syms_Q;

% 2. 旋转 0, 90, 180, 270 度并重新判决
ph_shifts = [0, pi/2, pi, 3*pi/2];
voted_4phases_bits = zeros(4, total_bits);
for p = 1:4
    rot_s = base_syms * exp(1j * ph_shifts(p));
    voted_4phases_bits(p, 1:2:end) = real(rot_s) < 0;
    voted_4phases_bits(p, 2:2:end) = imag(rot_s) < 0;
end

% 输出到文件
out_file = 'voted_consensus_4phases.txt';
fid = fopen(out_file, 'w');
fprintf(fid, '=== 2-Pass 强制对齐投票得出的核心序列与 4 种相位拓展 ===\n');
fprintf(fid, '基于全部 %d 个文件混合对抗提取的 "众数" 特征。\n\n', num_valid);

for p = 1:4
    fprintf(fid, 'Phase %d (%.0f deg): ', p, (p-1)*90);
    fprintf(fid, '%d', voted_4phases_bits(p, :));
    fprintf(fid, '\n');
end

% 顺便附上一张全员对齐后的总矩阵 (可选，便于人工验证)
fprintf(fid, '\n--- 附加信息：全部 %d 个文件对齐后投票前的原始数据 ---\n', num_valid);
for i = 1:num_valid
    fprintf(fid, 'File %03d: ', valid_indices(i));
    fprintf(fid, '%d', final_aligned_bits(i, :));
    fprintf(fid, '\n');
end
fclose(fid);

fprintf('强制投票 2-Pass 已完成！合并众数后的 4 种相位变体序列已保存至 %s\n', out_file);

%% --- 本地函数区 ---

function [f_center_pilot_hz, info] = estimate_center_pilot_blind_from_signal( ...
    x, Fs, target_bw_hz, bw_tol_ratio, occ_bg_win_hz, occ_smooth_hz, occ_thresh_db, min_component_bw_hz, ...
    max_expected_cfo_hz, pilot_bg_win_hz, pilot_min_prom_db, pilot_width_ref_hz, ...
    refine_enable, refine_lpf_bw_hz, refine_iters, refine_fir_order, refine_max_delta_hz)

segLen = 32768; overlapRatio = 0.75; nfft = 65536;
N = length(x);
if N < 2048, error('样本过短。'); end
if N < segLen * 2
    segLen = 2^floor(log2(max(1024, floor(N/2))));
    segLen = min(segLen, N); nfft = max(nfft, 2^nextpow2(segLen));
end
[P_med, ~, f] = robust_welch_psd_local(x, Fs, segLen, overlapRatio, nfft);
Pdb = 20*log10(P_med + eps);
bin_hz = Fs / length(f);
occ_bg_bins = make_odd_local(max(9, round(occ_bg_win_hz / bin_hz)));
occ_sm_bins = make_odd_local(max(5, round(occ_smooth_hz / bin_hz)));
P_occ_base = movmedian(Pdb, occ_bg_bins); P_occ_env = movmean(P_occ_base, occ_sm_bins);
noise_floor_db = prctile(P_occ_env, 20); occ_thr_db = noise_floor_db + occ_thresh_db;
occ_mask = P_occ_env > occ_thr_db;
components = mask_to_components_local(occ_mask);
cand = [];
for i = 1:size(components,1)
    iL = components(i,1); iR = components(i,2);
    bw = f(iR) - f(iL);
    if bw < min_component_bw_hz, continue; end
    bw_err = abs(bw - target_bw_hz) / target_bw_hz;
    in_mean = mean(P_occ_env(iL:iR));
    guard = round(1.5e6 / bin_hz);
    oL = max(1, iL-guard):max(1, iL-1); oR = min(length(P_occ_env), iR+1):min(length(P_occ_env), iR+guard);
    if isempty(oL), outL = noise_floor_db; else, outL = mean(P_occ_env(oL)); end
    if isempty(oR), outR = noise_floor_db; else, outR = mean(P_occ_env(oR)); end
    edge_drop = in_mean - 0.5*(outL + outR);
    score = 3.0*(1 - min(bw_err,1)) + 0.12*(in_mean - noise_floor_db) + 0.18*max(edge_drop,0);
    if bw_err <= bw_tol_ratio
        cand = [cand; iL, iR, score]; 
    end
end
if isempty(cand)
    Nf = length(f); bw_bins = max(9, round(target_bw_hz / bin_hz));
    half_bw = floor(bw_bins/2); out_half = min(floor(1.5*bw_bins), floor((Nf-1)/2));
    score_sw = -inf(Nf,1);
    for k = (out_half+2):(Nf-out_half-1)
        i1 = k-half_bw; i2 = k+half_bw;
        if i1 < 1 || i2 > Nf, continue; end
        in_band = P_occ_env(i1:i2);
        l1 = max(1, k-out_half); l2 = max(1, i1-1); r1 = min(Nf, i2+1); r2 = min(Nf, k+out_half);
        if l2 < l1 || r2 < r1, continue; end
        out_vals = [P_occ_env(l1:l2); P_occ_env(r1:r2)];
        contrast = mean(in_band) - mean(out_vals);
        score_sw(k) = contrast - 0.15*std(in_band);
    end
    [~, k_best] = max(score_sw); f_band_center = f(k_best);
else
    [~, id] = max(cand(:,3)); iL = cand(id,1); iR = cand(id,2);
    f_band_center = (f(iL) + f(iR)) / 2;
end
search_half_hz = max(1.5e6, max_expected_cfo_hz); search_half_hz = min(search_half_hz, target_bw_hz/2 - 0.5e6);
[~, f_subbin, pinfo] = detect_center_pilot_once_local(Pdb, f, f_band_center, search_half_hz, pilot_bg_win_hz, pilot_min_prom_db, pilot_width_ref_hz);
f_center_pilot_hz = f_subbin;
if refine_enable
    f_refined = refine_tone_freq_time_local(x, Fs, f_subbin, refine_lpf_bw_hz, refine_iters, refine_fir_order);
else
    f_refined = f_subbin;
end
refine_delta = f_refined - f_subbin;
if abs(refine_delta) > refine_max_delta_hz
    f_center_pilot_hz = f_subbin; refine_used = 0;
else
    f_center_pilot_hz = f_refined; refine_used = 1;
end
info = pinfo;
info.f_fft_subbin_hz = f_subbin; info.f_refined_hz = f_refined;
end

function [P_med, P_mean, f] = robust_welch_psd_local(x, Fs, segLen, overlapRatio, nfft)
x = x(:); N = length(x); nfftEff = 2^nextpow2(N); w = ones(N, 1);
X = fftshift(fft(x .* w, nfftEff)); P_out = abs(X); P_out = (P_out / max(P_out));
P_med = P_out; P_mean = P_out; f = ((-nfftEff/2):(nfftEff/2-1)).' * (Fs / nfftEff);
end

function [f_coarse, f_subbin, info] = detect_center_pilot_once_local(Pdb, f, f_center, search_half_hz, bg_win_hz, prom_min_db, width_ref_hz)
df = abs(f(2)-f(1)); idx = abs(f - f_center) <= search_half_hz; f_roi = f(idx); y = Pdb(idx);
if isempty(y), error('搜索区间为空。'); end
[max_val, ksel] = max(y); f_coarse = f_roi(ksel);
if ksel <= 1 || ksel >= length(y), f_subbin = f_coarse;
else
    y1 = y(ksel-1); y2 = y(ksel); y3 = y(ksel+1); den = (y1 - 2*y2 + y3);
    delta = (abs(den)<1e-12)*0 + (abs(den)>=1e-12)*(0.5*(y1-y3)/den);
    f_subbin = f_roi(ksel) + max(min(delta, 1), -1) * df;
end
info = struct('prom_db', max_val, 'width_hz', 0, 'score', max_val);
end
function comps = mask_to_components_local(mask)
mask = logical(mask(:)); d = diff([false; mask; false]);
comps = [find(d == 1), find(d == -1) - 1];
end
function n = make_odd_local(n)
n = max(3, round(n)); if mod(n,2)==0, n = n + 1; end
end
function f_refined = refine_tone_freq_time_local(x, Fs, f_init_hz, lpf_bw_hz, iters, fir_order)
x = x(:); N = length(x); n = (0:N-1).'; f_refined = f_init_hz;
for it = 1:max(1, iters)
    z = x .* exp(-1j * 2*pi * f_refined * n / Fs);
    Wn = min(0.95, max(1e-5, (lpf_bw_hz / (2^(it-1))) / (Fs/2)));
    b = fir1(fir_order, Wn); zf = filtfilt(b, 1, z);
    trim = min(length(zf)-2, max(8, fir_order));
    z_use = zf; if trim*2 < length(zf)-1, z_use = zf(trim+1:end-trim); end
    if length(z_use) < 4, break; end
    r = conj(z_use(1:end-1)) .* z_use(2:end);
    df_hz = angle(sum(r)) * Fs / (2*pi);
    f_refined = f_refined + df_hz;
    if abs(df_hz) < 5, break; end
end
end

%% ===== 4QAM 盲解调核心函数 =====
function [x_corr, slope_best, phi_best, cost_best] = blind_phase_slope_comp_4qam(x_full, y_valid, k_valid, k_full)
y_valid = y_valid(:); k_valid = k_valid(:);
if numel(y_valid) < 8
    slope_best = 0; phi_best = 0; cost_best = inf; x_corr = x_full; return;
end
y_norm = normalize_to_4qam_power(y_valid);

slope_grid = linspace(-0.35, 0.35, 141); 
phi_grid   = linspace(-pi, pi, 181);
[slope_best, phi_best, cost_best] = search_cost_4qam(y_norm, k_valid, slope_grid, phi_grid);

slope_fine = linspace(slope_best - 0.02, slope_best + 0.02, 121);
phi_fine   = linspace(phi_best - pi/12, phi_best + pi/12, 121);
[slope_best, phi_best, cost_best] = search_cost_4qam(y_norm, k_valid, slope_fine, phi_fine);

phase_full = slope_best * k_full + phi_best;
x_corr = x_full .* exp(-1j * phase_full);
end

function [slope_best, phi_best, cost_best] = search_cost_4qam(y_norm, k_valid, slope_grid, phi_grid)
cost_best = inf; slope_best = 0; phi_best = 0;
for is = 1:length(slope_grid)
    y_s = y_norm .* exp(-1j * slope_grid(is) * k_valid);
    for ip = 1:length(phi_grid)
        p = phi_grid(ip);
        z = y_s * exp(-1j * p);
        z_hat = slicer_4qam(z);
        d = z - z_hat;
        e2 = abs(d).^2;
        thr = prctile(e2, 90);
        c = mean(e2(e2 <= thr));
        if c < cost_best
            cost_best = c; slope_best = slope_grid(is); phi_best = p;
        end
    end
end
phi_best = mod(phi_best + pi, 2*pi) - pi;
end

function z_hat = slicer_4qam(z)
I = sign(real(z)); I(I==0) = 1;
Q = sign(imag(z)); Q(Q==0) = 1;
z_hat = I + 1j*Q;
end

function y = normalize_to_4qam_power(x)
Es_target = 2;
g = sqrt(Es_target / (mean(abs(x).^2) + eps));
y = x * g;
end