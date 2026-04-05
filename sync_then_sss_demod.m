% sync_then_sss_demod.m
% 同步-解调一体脚本：
% 1) 先复刻 frame_sync_routine 的归一化匹配找峰，自动得到 burst 起始点
% 2) 将该起始点作为 SSS 解调的 read_start_sample
% 3) 执行 SSS 单点解调（支持可选 CFO 补偿）

clear; clc; close all;

%% ===================== 0. 全局参数 =====================
inFile = 'sigtest8.iq';

% --- 同步扫描参数（复刻 frame_sync_routine） ---
sync_scan_start = 0;      % 0-based
sync_scan_len   = 20e6;   % 最多读取多少点用于找峰
min_peak_height = 50;     % 峰值门限（%）
peak_select_mode = 'first'; % 'first' or 'strongest'

% --- SSS 解调参数 ---
read_length = 6992*3 + 1000;   % 从 read_start_sample 开始读取长度
sss_decode_start_idx = 1024 + 48;
target_offset = 0;

fs_source = 409.6e6;
fs_target = 60e6;
freq_shift_hz = 63.5e6;   % DDC

% CFO 补偿参数（可选）
enable_cfo_comp = true;   % true: 启用固定 CFO 补偿；false: 跳过
cfo_hz = -367188;               % 频偏(Hz)，可回填来自 gardner/其他估计

sro_ppm = 0;              % 若已知可回填

N_fft = 1024;
demod_all_subcarriers = true;
decision_rotate_deg = 45;

%% ===================== 1. 构造本地 PSS 模板 =====================
fprintf('Step 1: 生成本地 PSS 模板...\n');
fs_symbol = 60e6;
fc = 63e6;

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

[P, Q] = rat(fs_source / fs_symbol);
pss_up = resample(pss_base_with_cp, P, Q);
t_local = (0:length(pss_up)-1) / fs_source;
pss_local = pss_up .* exp(1j * 2 * pi * fc * t_local);
pss_local = pss_local / max(abs(pss_local));

%% ===================== 2. 全局同步找峰 =====================
fprintf('Step 2: 全局归一化匹配找峰...\n');

fInfo = dir(inFile);
if isempty(fInfo)
    error('找不到文件: %s', inFile);
end

totalSamples = floor(fInfo.bytes / 4);
if sync_scan_start < 0 || sync_scan_start >= totalSamples
    error('sync_scan_start 越界: %d / total=%d', sync_scan_start, totalSamples);
end

sync_scan_len = min(round(sync_scan_len), totalSamples - sync_scan_start);
[x_sync, meta_sync] = iq_read_int16_le(inFile, sync_scan_start, sync_scan_len);
if meta_sync.numSamplesRead <= 0
    error('同步扫描阶段未读取到有效数据。');
end

x_sync = double(x_sync(:));
x_sync = x_sync - mean(x_sync);
x_sync = x_sync / (mean(abs(x_sync)) + eps);

% 归一化匹配滤波
h_matched = fliplr(conj(pss_local));
corr_out = filter(h_matched, 1, x_sync);
corr_mag_sq = abs(corr_out).^2;

E_local = sum(abs(pss_local).^2);
win_ones = ones(1, length(pss_local));
E_rx_moving = filter(win_ones, 1, abs(x_sync).^2);
E_rx_moving(E_rx_moving < 1e-10) = 1e-10;

corr_norm_pct = (corr_mag_sq ./ (E_local .* E_rx_moving)) * 100;

min_peak_dist = length(pss_local);
[pks, locs] = findpeaks(corr_norm_pct, 'MinPeakHeight', min_peak_height, 'MinPeakDistance', min_peak_dist);

if isempty(locs)
    error('未检测到有效同步峰。请降低 min_peak_height 或检查输入信号。');
end

L_pss = length(pss_local);
end_pos_global = sync_scan_start + (locs - 1);
start_pos_global = end_pos_global - (L_pss - 1);

fprintf('  检测到 %d 个峰：\n', length(locs));
for i = 1:length(locs)
    fprintf('    Peak %02d: %.2f%%, 末端=%d, 起始=%d\n', i, pks(i), end_pos_global(i), start_pos_global(i));
end

switch lower(peak_select_mode)
    case 'first'
        pick_idx = 1;
    case 'strongest'
        [~, pick_idx] = max(pks);
    otherwise
        error('未知 peak_select_mode: %s', peak_select_mode);
end

read_start_sample = start_pos_global(pick_idx);
if read_start_sample < 0
    error('选中的起始点为负数: %d（请检查峰选择）', read_start_sample);
end

fprintf('\n>>> 选中峰 #%d, 匹配度 %.2f%%\n', pick_idx, pks(pick_idx));
fprintf('>>> 自动设置 read_start_sample = %d\n', read_start_sample);

% 同步结果图
figure('Position', [80, 80, 1400, 520], 'Name', 'Sync Peaks (Normalized Correlation)');
sample_axis_sync = sync_scan_start + (0:length(corr_norm_pct)-1);
plot(sample_axis_sync, corr_norm_pct, 'b', 'LineWidth', 1.1); hold on; grid on;
plot(sync_scan_start + (locs - 1), pks, 'r*', 'MarkerSize', 8, 'LineWidth', 1.2);
yline(min_peak_height, 'g--', 'LineWidth', 1.2);
xlabel('复采样点索引 (Complex Sample Index)');
ylabel('归一化相关系数 Match (%)');
title(sprintf('Sync Peaks (detected=%d, selected start=%d)', length(locs), read_start_sample));

%% ===================== 3. SSS 解调（使用自动起始点） =====================
fprintf('\nStep 3: 使用自动 read_start_sample 进行 SSS 解调...\n');
fprintf('Loading file: %s from %d, len %d...\n', inFile, read_start_sample, read_length);
[x_raw, ~] = iq_read_int16_le(inFile, read_start_sample, read_length);
x_raw = double(x_raw);

if isempty(x_raw)
    error('SSS 解调阶段读取为空。');
end

x_raw = x_raw - mean(x_raw);
x_raw = x_raw / (mean(abs(x_raw)) + eps);

%% 3.1 DDC + 可选 CFO 补偿
t_vec = (0:length(x_raw)-1) / fs_source;
fprintf('Shifting spectrum left by %.1f MHz at 409.6MHz...\n', freq_shift_hz/1e6);
x_shifted = x_raw .* exp(-1j * 2 * pi * freq_shift_hz * t_vec).';

if enable_cfo_comp
    fprintf('Applying CFO correction: %.2f Hz...\n', cfo_hz);
    x_base = x_shifted .* exp(-1j * 2 * pi * cfo_hz * t_vec).';
else
    fprintf('CFO correction disabled.\n');
    x_base = x_shifted;
end

%% 3.2 SRO 修正（Farrow）
fprintf('Applying SRO Correction (Farrow): %.2f ppm...\n', sro_ppm);
fs_eff = fs_source * (1 + sro_ppm/1e6);

Wn = 35e6 / (fs_source/2);
lpf_order = 30;
b_lpf = fir1(lpf_order, Wn);
x_filt = filtfilt(b_lpf, 1, x_base);

T_in = 1 / fs_eff;
T_out = 1 / fs_target;
t_out = 0 : T_out : (length(x_filt)-3)*T_in;

idx_frac = t_out / T_in + 1;
idx_base = floor(idx_frac);
mu = idx_frac - idx_base;

valid_mask = (idx_base >= 2) & (idx_base <= length(x_filt)-2);
idx_base = idx_base(valid_mask);
mu = mu(valid_mask);

h0 = -(mu - 1) .* (mu - 2) .* mu / 6;
h1 =  (mu - 1) .* (mu + 1) .* (mu - 2) / 2;
h2 = -(mu + 1) .* mu .* (mu - 2) / 2;
h3 =  (mu + 1) .* (mu - 1) .* mu / 6;

x_sro = h0 .* x_filt(idx_base - 1) + ...
        h1 .* x_filt(idx_base) + ...
        h2 .* x_filt(idx_base + 1) + ...
        h3 .* x_filt(idx_base + 2);
x_sro = x_sro(:);

%% 3.3 单点 SSS 频域解调
base_start_idx = sss_decode_start_idx;
offset = target_offset;
sss_start_idx_60 = base_start_idx + offset;

if sss_start_idx_60 < 1 || sss_start_idx_60 + N_fft > length(x_sro)
    error('SSS 起始点越界: idx=%d, len=%d', sss_start_idx_60, length(x_sro));
end

x_sss_time = x_sro(sss_start_idx_60 : sss_start_idx_60 + N_fft - 1);
x_sss_freq = fft(x_sss_time, N_fft) / sqrt(N_fft);

% 动态有效载波检测
pwr_bins = abs(x_sss_freq).^2;
mean_pwr = mean(pwr_bins);
threshold_pwr = mean_pwr * 0.1;
valid_mask_power = pwr_bins > threshold_pwr;

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

% 分段 unwrap + 分段 polyfit
syms_pow4 = syms_valid.^4;
df = diff(freq_indices);
seg_start = [1; find(df > 1) + 1];
seg_end = [find(df > 1); length(freq_indices)];

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
    slope_list(end+1,1) = p_seg(1); %#ok<SAGROW>
    w_list(end+1,1) = numel(idx_seg); %#ok<SAGROW>
end

if isempty(slope_list)
    error('有效载波不足，无法完成相位拟合。');
end

slope4 = sum(slope_list .* w_list) / sum(w_list);
pow4_detrended = syms_pow4 .* exp(-1j * slope4 * freq_indices);
phase0_4 = angle(mean(pow4_detrended));

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

%% 3.4 QPSK 4-way 候选
rot_angles_deg = decision_rotate_deg + [0, 90, 180, 270];
hex_candidates = cell(1,4);
for r = 1:4
    syms_payload_rot = syms_payload .* exp(1j * decision_rotate_deg * pi/180) .* exp(-1j * (r-1) * pi/2);
    bits_I = real(syms_payload_rot) < 0;
    bits_Q = imag(syms_payload_rot) < 0;

    demod_bits = zeros(length(syms_payload_rot)*2, 1);
    demod_bits(1:2:end) = bits_I;
    demod_bits(2:2:end) = bits_Q;

    full_hex_str = '';
    for i_hex = 1:4:length(demod_bits)
        if i_hex + 3 > length(demod_bits)
            chunk = demod_bits(i_hex:end);
            val = 0;
            for b = 1:length(chunk)
                val = val + chunk(b) * 2^(length(chunk)-b);
            end
            full_hex_str = [full_hex_str, dec2hex(val)]; %#ok<AGROW>
            break;
        end
        chunk = demod_bits(i_hex:i_hex+3);
        val = chunk(1)*8 + chunk(2)*4 + chunk(3)*2 + chunk(4);
        full_hex_str = [full_hex_str, dec2hex(val)]; %#ok<AGROW>
    end
    hex_candidates{r} = full_hex_str;
end

%% 3.5 可视化与输出
figure('Position', [300, 300, 640, 560], 'Name', 'SSS Demod (Sync-driven read_start_sample)');
subplot(1,2,1);
syms_plot = x_sss_freq_corr(:) .* exp(1j * decision_rotate_deg * pi/180);
plot(real(syms_plot), imag(syms_plot), 'b.', 'MarkerSize', 8); hold on; grid on; axis square;
th = 0:0.01:2*pi;
plot(cos(th), sin(th), 'k:', 'LineWidth', 0.5);
title(sprintf('Constellation (Rot %d°)', decision_rotate_deg));
xlabel('I'); ylabel('Q'); xlim([-2 2]); ylim([-2 2]);

subplot(1,2,2); axis off;
text(0.02, 0.82, sprintf('Auto read_start = %d\nOffset=%d\n%s', read_start_sample, offset, payload_label), 'FontSize', 11);
text(0.02, 0.54, sprintf('0°  : %s...', hex_candidates{1}(1:min(22,end))), 'FontSize', 10, 'Interpreter','none');
text(0.02, 0.42, sprintf('90° : %s...', hex_candidates{2}(1:min(22,end))), 'FontSize', 10, 'Interpreter','none');
text(0.02, 0.30, sprintf('180°: %s...', hex_candidates{3}(1:min(22,end))), 'FontSize', 10, 'Interpreter','none');
text(0.02, 0.18, sprintf('270°: %s...', hex_candidates{4}(1:min(22,end))), 'FontSize', 10, 'Interpreter','none');

fprintf('\n>>> [Auto Sync + SSS Demod] 完成 <<<\n');
fprintf('Auto read_start_sample = %d\n', read_start_sample);
fprintf('Extracted bits per candidate = %d (%s)\n', length(syms_payload)*2, payload_label);
fprintf('====== HEX CANDIDATES (QPSK 4-way rotation ambiguity) ======\n');
for r = 1:4
    fprintf('[%3d deg] %s\n', rot_angles_deg(r), hex_candidates{r});
end
fprintf('============================================================\n');
