% batch_16qam_vote_173.m
% 批量 173文件 SSS 16QAM 解调扫描与跨文件不变比特投票统计
% 1. -5 到 30 offset遍历，根据 16QAM Cost 函数取最佳 offset
% 2. 提取 4 种相位候选角度（0, 90, 180, 270）并映射为 Bit
% 3. 跨 173 个样本统计出高度一致的比特并输出位置

clear; clc; close all;

%% 1. 全局与批处理参数
num_files = 173;        % 循环读取到 sigtest173.iq
file_prefix = 'sigtest';
file_ext = '.iq';

offset_range = -5:30;   % Requirement 3: 从-5扫到30
read_length = 6992*3+1000;
sss_decode_start_idx = 1024+48+1024+48;

fs_source = 409.6e6;
fs_target = 60e6;
sro_ppm = 0;
N_fft = 1024;
demod_all_subcarriers = true; 

% ========================================================
% 建立存储变量
if demod_all_subcarriers
    payload_len_bits = 1024 * 4; % 16QAM = 4 bits per carrier
else
    % 如果有效载波不定长，可根据实际截取。这里以固定维数组容纳
    payload_len_bits = 800 * 4; 
end

% 维度: [文件数 x 4个候选相位 x 比特数]
all_cands_bits = NaN(num_files, 4, payload_len_bits);
valid_file_flags = false(num_files, 1);

%% 2. 批量解调与候选比特提取
for f_idx = 1:num_files
    inFile = sprintf('%s%d%s', file_prefix, f_idx, file_ext);
    if ~isfile(inFile)
        fprintf('File NOT found: %s, skipping...\n', inFile);
        continue;
    end
    
    fprintf('\n===== Processing File [%d/%d]: %s =====\n', f_idx, num_files, inFile);
    
    % --- 2.1 盲估中心频率 & 匹配滤波 ---
    sync_scan_start = 0;
    sync_scan_len   = 5e6;
    try
        [x_scan, meta_scan] = iq_read_int16_le(inFile, sync_scan_start, sync_scan_len);
    catch
        fprintf('Failed to read %s\n', inFile);
        continue;
    end
    if meta_scan.numSamplesRead <= 0, continue; end
    
    x_scan = double(x_scan(:));
    x_scan = x_scan - mean(x_scan);
    x_scan = x_scan / (mean(abs(x_scan)) + eps);

    [freq_shift_hz, ~] = estimate_center_pilot_blind_from_signal( ...
        x_scan, fs_source, 60e6, 0.35, 2.0e6, 1.2e6, 12.0, 20e6, ...
        6e6, 0.9e6, 4.5, 0.25e6, false, 0.6e6, 2, 256, 120e3);
    
    fc_sync = 63e6;
    fs_symbol = fs_target;
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
    
    corr_norm_pct = (abs(corr_out).^2 ./ (sum(abs(pss_local).^2) .* max(1e-10, filter(ones(1, length(pss_local)), 1, abs(x_scan).^2)))) * 100;
    [~, locs_all] = findpeaks(corr_norm_pct, 'MinPeakHeight', 50, 'MinPeakDistance', length(pss_local));
    if isempty(locs_all), [~, locs_all] = max(corr_norm_pct); end
    [~, max_pk_idx] = max(corr_norm_pct(locs_all));
    read_start_sample = locs_all(max_pk_idx) - length(pss_local);
    
    % --- 2.2 读取目标帧数据并预处理 ---
    [x_raw, ~] = iq_read_int16_le(inFile, read_start_sample, read_length);
    x_raw = double(x_raw);
    x_raw = (x_raw - mean(x_raw)) / mean(abs(x_raw));

    t_vec = (0:length(x_raw)-1) / fs_source;
    x_shifted = x_raw .* exp(-1j * 2 * pi * freq_shift_hz * t_vec).';
    x_cfo = x_shifted .* exp(-1j * 2 * pi * 0 * t_vec).'; % cfo_hz = 0 since estimated blindly
    
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
    idx_base = idx_base(valid_mask); mu = mu(valid_mask);
    
    h0 = -(mu - 1) .* (mu - 2) .* mu / 6;
    h1 =  (mu - 1) .* (mu + 1) .* (mu - 2) / 2;
    h2 = -(mu + 1) .* mu .* (mu - 2) / 2;
    h3 =  (mu + 1) .* (mu - 1) .* mu / 6;
    x_sro = h0 .* x_cfo_filtered(idx_base - 1) + h1 .* x_cfo_filtered(idx_base) + h2 .* x_cfo_filtered(idx_base + 1) + h3 .* x_cfo_filtered(idx_base + 2);
    x_sro = x_sro(:);

    % --- 2.3 扫偏 (Offset) 提取最优 Constellation ---
    best_cost = inf;
    best_syms_payload = [];
    
    for offset = offset_range
        idx_start = sss_decode_start_idx + offset;
        if idx_start < 1 || idx_start + N_fft > length(x_sro), continue; end
        
        x_sss_freq = fft(x_sro(idx_start : idx_start + N_fft - 1), N_fft) / sqrt(N_fft);
        
        pwr_bins = abs(x_sss_freq).^2;
        valid_mask_power = pwr_bins > (mean(pwr_bins) * 0.1);
        valid_mask_power(1:3) = false; valid_mask_power(N_fft-2:N_fft) = false;
        valid_idxs = find(valid_mask_power);
        if isempty(valid_idxs), continue; end

        rel_freq = zeros(N_fft, 1);
        rel_freq(1:N_fft/2) = (0:N_fft/2 - 1).'; rel_freq(N_fft/2+1:N_fft) = (-N_fft/2:-1).';
        freq_indices = rel_freq(valid_idxs);
        syms_valid = x_sss_freq(valid_idxs);
        [freq_indices, sort_idx] = sort(freq_indices);
        syms_valid = syms_valid(sort_idx);
        
        full_freq_indices = [(0:N_fft/2-1), (-N_fft/2:-1)].';
        [x_sss_freq_corr, ~, ~, cost_best] = blind_phase_slope_comp_16qam(x_sss_freq, syms_valid, freq_indices, full_freq_indices);
        
        if cost_best < best_cost
            best_cost = cost_best;
            if demod_all_subcarriers
                best_syms_payload = x_sss_freq_corr(:);
            else
                best_syms_payload = x_sss_freq_corr(valid_idxs(sort_idx));
            end
        end
    end
    
    if isempty(best_syms_payload)
        fprintf('Failed to find valid offset constellation.\n');
        continue; 
    end
    
    best_syms_payload = normalize_to_16qam_power(best_syms_payload);

    % --- 2.4 生成 4 个候选相位的比特流 (90度相位模糊) ---
    phase_shifts = [0, pi/2, pi, 3*pi/2];
    for cand_idx = 1:4
        rotated_syms = best_syms_payload * exp(1j * phase_shifts(cand_idx));
        bits = qam16_demap(rotated_syms);
        % 将比特流统一裁切填入矩阵中
        store_len = min(length(bits), payload_len_bits);
        all_cands_bits(f_idx, cand_idx, 1:store_len) = bits(1:store_len);
    end
    valid_file_flags(f_idx) = true;
    fprintf('File %d demodulated. Best cost: %.4f\n', f_idx, best_cost);
end

%% 3. 跨 173 样本的 “模糊解缠与少数一致比特抽取代”
fprintf('\n===== 开始统计跨文件的小部分相同比特位置 =====\n');
valid_files = find(valid_file_flags);
num_valid = length(valid_files);
if num_valid < 2
    error('有效解析文件数少于2，无法统计公用比特。');
end

% 由于“大部分是不一样的，小部分是一样的”，传统的直接求相似度可能被占多数的随机比特淹没。
% 为规避此问题，我们基于多数表决的方式强行挑出第一条样本为基准的相对最佳匹配：
usable_bits_len = size(all_cands_bits, 3);
mapped_cands_bits = NaN(num_valid, usable_bits_len);

% [算法] 2 pass 对齐与发现
% Pass 1: 强行拿基准进行全尺寸汉明距离比对，找到大致对其的序列。
ref_bits = squeeze(all_cands_bits(valid_files(1), 1, :)).';
mapped_cands_bits(1, :) = ref_bits;

for i = 2:num_valid
    f_idx = valid_files(i);
    best_dist = inf; best_cand = 1;
    for cand = 1:4
        cb = squeeze(all_cands_bits(f_idx, cand, :)).';
        dist = sum(cb ~= ref_bits, 'omitnan');
        if dist < best_dist
            best_dist = dist; best_cand = cand;
        end
    end
    mapped_cands_bits(i, :) = squeeze(all_cands_bits(f_idx, best_cand, :)).';
end

% Pass 1 的比特均值（出现1的概率）
bit_prob_pass1 = mean(mapped_cands_bits, 1, 'omitnan');
% 判断哪些bit位是高度“确信”的 (概率极高或极低)
stable_thresh = 0.90; % 90% 以上一致则认定为是不变部分
stable_mask = (bit_prob_pass1 >= stable_thresh) | (bit_prob_pass1 <= (1-stable_thresh));

% Pass 2: 如果由于大量的随机噪声导致 Pass1 锁错相位，我们仅利用提取出来的稳定小段进行二次重对齐
fprintf('Found %d potentially stable bits in Pass 1. Re-aligning candidates...\n', sum(stable_mask));
if sum(stable_mask) > 10
    ref_bits_stable = ref_bits(stable_mask);
    for i = 2:num_valid
        f_idx = valid_files(i);
        best_dist = inf; best_cand = 1;
        for cand = 1:4
            cb = squeeze(all_cands_bits(f_idx, cand, :)).';
            dist = sum(cb(stable_mask) ~= ref_bits_stable, 'omitnan');
            if dist < best_dist
                best_dist = dist; best_cand = cand;
            end
        end
        mapped_cands_bits(i, :) = squeeze(all_cands_bits(f_idx, best_cand, :)).';
    end
end

% --- 计算最终的一致性分布 ---
final_bit_prob = mean(mapped_cands_bits, 1, 'omitnan');

% 最终筛选出高度完全相同的比特
final_thresh = 0.95; % 严格门限，95%的文件在此位输出相同比特
identical_zeros = find(final_bit_prob <= (1 - final_thresh));
identical_ones  = find(final_bit_prob >= final_thresh);

fprintf('\n==== 统计结果汇总 (%d 个有效文件) ====\n', num_valid);
fprintf('发现高度一致为 "0" 的比特位置 共 %d 个\n', length(identical_zeros));
if ~isempty(identical_zeros)
    % 控制台只打印一部分防刷屏
    disp(identical_zeros(1:min(30, end)));
end

fprintf('\n发现高度一致为 "1" 的比特位置 共 %d 个\n', length(identical_ones));
if ~isempty(identical_ones)
    disp(identical_ones(1:min(30, end)));
end

% 把结果写出为文件供查阅
fid_out = fopen('common_bits_positions.txt', 'w');
fprintf(fid_out, '===== Common Bits across %d samples =====\n\n', num_valid);
fprintf(fid_out, 'Always 0 locations (Total %d):\n', length(identical_zeros));
fprintf(fid_out, '%d ', identical_zeros);
fprintf(fid_out, '\n\nAlways 1 locations (Total %d):\n', length(identical_ones));
fprintf(fid_out, '%d ', identical_ones);
fclose(fid_out);
fprintf('\n完整的一致性位置列表已写入: common_bits_positions.txt\n');


%% ================= 本地函数区 =================

function bits = qam16_demap(z)
    % 标准 16QAM 的硬判决到 Bit 流转换 (每个符号 4 个比特)
    % 假设 IQ 实部虚部分别按 [-3, -1, 1, 3] 分布
    I = real(z);
    Q = imag(z);
    
    b1 = I > 0;
    b2 = abs(I) < 2;
    b3 = Q > 0;
    b4 = abs(Q) < 2;
    
    N = length(z);
    bits = zeros(N*4, 1);
    bits(1:4:end) = b1;
    bits(2:4:end) = b2;
    bits(3:4:end) = b3;
    bits(4:4:end) = b4;
end

function [f_center_pilot_hz, info] = estimate_center_pilot_blind_from_signal( ...
    x, Fs, target_bw_hz, bw_tol_ratio, occ_bg_win_hz, occ_smooth_hz, occ_thresh_db, min_component_bw_hz, ...
    max_expected_cfo_hz, pilot_bg_win_hz, pilot_min_prom_db, pilot_width_ref_hz, ...
    refine_enable, refine_lpf_bw_hz, refine_iters, refine_fir_order, refine_max_delta_hz)
    
    % 省略了复杂的子模块，直接保留原始精简逻辑（具体可复用你前置代码）
    nfft = 65536; segLen = 32768; 
    N = length(x); if N < segLen * 2, segLen = min(2^floor(log2(floor(N/2))), N); nfft = max(nfft, 2^nextpow2(segLen)); end
    
    w = ones(N, 1); X = fftshift(fft(x .* w, nfft));
    P_out = abs(X); P_out = P_out / max(P_out); Pdb = 20*log10(P_out + eps);
    f = ((-nfft/2):(nfft/2-1)).' * (Fs / nfft);
    
    % 这里强制简化直接拿最高峰 (等同你之前的单音定位核心)
    [~, max_idx] = max(Pdb);
    f_center_pilot_hz = f(max_idx);
    info = struct();
end

function [x_corr, slope_best, phi_best, cost_best] = blind_phase_slope_comp_16qam(x_full, y_valid, k_valid, k_full)
    y_valid = y_valid(:); k_valid = k_valid(:);
    if numel(y_valid) < 16, slope_best = 0; phi_best = 0; cost_best = inf; x_corr = x_full; return; end
    y_norm = normalize_to_16qam_power(y_valid);
    
    % -------- 粗搜索 --------
    slope_grid = linspace(-0.35, 0.35, 141);     % rad/bin
    phi_grid   = linspace(-pi, pi, 181);         % rad
    [slope_best, phi_best, cost_best] = search_cost(y_norm, k_valid, slope_grid, phi_grid);
    
    % -------- 细搜索 --------
    slope_fine = linspace(slope_best - 0.02, slope_best + 0.02, 121);
    phi_fine   = linspace(phi_best - pi/12, phi_best + pi/12, 121);
    [slope_best, phi_best, cost_best] = search_cost(y_norm, k_valid, slope_fine, phi_fine);
    
    phase_full = slope_best * k_full + phi_best;
    x_corr = x_full .* exp(-1j * phase_full);
end

function [slope_best, phi_best, cost_best] = search_cost(y_norm, k_valid, slope_grid, phi_grid)
    cost_best = inf; slope_best = 0; phi_best = 0;
    for s = slope_grid
        y_s = y_norm .* exp(-1j * s * k_valid);
        for p = phi_grid
            z = y_s * exp(-1j * p);
            d = abs(z - slicer_16qam(z)).^2;
            c = mean(d(d <= prctile(d, 90)));
            if c < cost_best, cost_best = c; slope_best = s; phi_best = p; end
        end
    end
    phi_best = mod(phi_best + pi, 2*pi) - pi;
end

function z_hat = slicer_16qam(z)
    q = @(x) -3.*(x<-2) - 1.*(x>=-2&x<0) + 1.*(x>=0&x<2) + 3.*(x>=2);
    z_hat = q(real(z)) + 1j*q(imag(z));
end

function y = normalize_to_16qam_power(x)
    y = x * sqrt(10 / (mean(abs(x).^2) + eps));
end

function [x, meta] = iq_read_int16_le(filePath, startIdx, count)
    % 本地兼容旧读取函数
    f = fopen(filePath, 'r');
    if f < 0, meta.numSamplesRead = 0; x = []; return; end
    fseek(f, startIdx * 4, 'bof');
    raw = fread(f, count * 2, 'int16'); fclose(f);
    if isempty(raw), meta.numSamplesRead = 0; x = []; return; end
    x = raw(1:2:end) + 1j * raw(2:2:end);
    meta.numSamplesRead = length(x);
end
