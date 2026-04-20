% sss_demodulation_sweep_16qam.m
% SSS 偏移量遍历试验脚本（16QAM版本，仅用于 sweep/观察）
% 说明：参数保持与原脚本一致；不输出误码率/比特序列，仅观察星座图

clear; clc; close all;

%% 1. 参数设置
inFile = 'sigtest1.iq';

% ===== 手动填写的 3 个关键参数 =====
read_start_sample = 14472; % [参数1] 文件中开始读取的点（原始采样点）
read_length = 6992*3+1000;                  % [参数2] 读取长度（原始采样点个数）
sss_decode_start_idx = 1024+48+1024+48;           % [参数3] 从重采样后序列 x_sro 的第几个点开始作为 SSS 解调基准

fs_source = 409.6e6;
fs_target = 60e6;

sro_ppm = 0;
cfo_hz = 17556;

N_fft = 1024;
freq_shift_hz = 63e6;

offset_range = -5:5;   % 遍历范围
pause_each = false;     % 每个 offset 显示后是否暂停 (挑选模式设为 false 以便快速扫完)

% 解调模式：true=解调全部1024子载波；false=仅解调有效子载波
demod_all_subcarriers = true;

%% 1.5 盲估中心频率 & 匹配滤波寻找 PSS 同步点
fprintf('\n--- 自动盲估 & 同步阶段 ---\n');
sync_scan_start = 0;
sync_scan_len   = 5e6;
[x_scan, meta_scan] = iq_read_int16_le(inFile, sync_scan_start, sync_scan_len);
if meta_scan.numSamplesRead <= 0
    error('同步扫描数据读取失败');
end
x_scan = double(x_scan(:));
x_scan = x_scan - mean(x_scan);
x_scan = x_scan / (mean(abs(x_scan)) + eps);

fprintf('1. 正在通过大点数FFT盲估中心单音导频频率...\n');
[est_f_center, ~] = estimate_center_pilot_blind_from_signal( ...
    x_scan, fs_source, 60e6, 0.35, 2.0e6, 1.2e6, 12.0, 20e6, ...
    6e6, 0.9e6, 4.5, 0.25e6, false, 0.6e6, 2, 256, 120e3);
freq_shift_hz = est_f_center; % 覆盖顶部的手动输入
cfo_hz = 0;                   % 盲估一次性找准了，所以外加CFO设0
fprintf('-> 确认中心导频频率 = %.6f MHz\n', freq_shift_hz / 1e6);

fprintf('2. 正在生成本地 PSS 模板并执行归一化匹配滤波...\n');
fs_symbol = fs_target;
fc_sync = 63e6; % 保持和原脚本完全一致的 PSS 同步载频 63e6
hex_str = 'BD3CD0148871751F84CED8C1BE32AC96';
hex_chars = char(hex_str);
bits_dummy = zeros(1, 128);
for i = 1:length(hex_chars)
    bits_dummy((i-1)*4 + 1 : i*4) = dec2bin(hex2dec(hex_chars(i)), 4) - '0';
end
d_phi = (bits_dummy == 1) * (-pi/2) + (bits_dummy == 0) * (pi/2);
m_syms = exp(1j * cumsum(d_phi));
pss_base = [-m_syms, m_syms, m_syms, m_syms, m_syms, m_syms, m_syms, m_syms];
cp_syms = -pss_base(end - 47 : end);
pss_base_with_cp = [cp_syms, pss_base];

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
    fprintf('警告：未找到有效 PSS 同步峰，回退到全局最高点。\n');
    [~, max_idx] = max(corr_norm_pct);
    locs_all = max_idx;
end
[~, max_pk_idx] = max(corr_norm_pct(locs_all));
best_peak_loc = locs_all(max_pk_idx);

% 校正绝对采样点坐标（0-based）
read_start_sample = best_peak_loc - length(pss_local);
fprintf('-> 同步成功！锁定信号起播点 (read_start_sample) = %d\n', read_start_sample);
fprintf('------------------------\n\n');

%% 2. 读取原始数据
fprintf('Loading file: %s from %d, len %d...\n', inFile, read_start_sample, read_length);
[x_raw, ~] = iq_read_int16_le(inFile, read_start_sample, read_length);
x_raw = double(x_raw);

x_raw = x_raw - mean(x_raw);
x_raw = x_raw / mean(abs(x_raw));

%% 3. DDC + CFO
t_vec = (0:length(x_raw)-1) / fs_source;

fprintf('Shifting spectrum left by %.1f MHz at 409.6MHz...\n', freq_shift_hz/1e6);
x_shifted = x_raw .* exp(-1j * 2 * pi * freq_shift_hz * t_vec).';

fprintf('Applying separately CFO Correction: %.2f Hz...\n', cfo_hz);
x_cfo = x_shifted .* exp(-1j * 2 * pi * cfo_hz * t_vec).';

%% 4. SRO/Farrow
fprintf('Applying SRO Correction (Farrow Interpolation): %.2f ppm...\n', sro_ppm);
fs_eff = fs_source * (1 + sro_ppm/1e6);

fprintf('Applying Zero-Phase Anti-aliasing LPF (35MHz)...\n');
Wn = 35e6 / (fs_source / 2);
lpf_order = 30;
b_lpf = fir1(lpf_order, Wn);
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

%% 5. 遍历 offset（16QAM 盲相位斜率补偿 + 星座观察）
fprintf('Analyzing Time Domain Signal for 16QAM...\n');
base_start_idx = sss_decode_start_idx;

fig_const = figure('Position', [320, 220, 760, 660], 'Name', 'SSS Demodulation Sweep (16QAM)');

% 记录最佳解结果
best_offset = NaN;
best_cost = inf;
best_syms_plot = [];
best_slope_hat = 0;
best_phi_hat = 0;

for offset = offset_range
    sss_start_idx_60 = base_start_idx + offset;

    if sss_start_idx_60 < 1 || sss_start_idx_60 + N_fft > length(x_sro)
        warning('Offset %d 导致下标越界，跳过。', offset);
        continue;
    end

    x_sss_time = x_sro(sss_start_idx_60 : sss_start_idx_60 + N_fft - 1);
    x_sss_freq = fft(x_sss_time, N_fft) / sqrt(N_fft);

    % 有效载波筛选（沿用原门限逻辑）
    pwr_bins = abs(x_sss_freq).^2;
    mean_pwr = mean(pwr_bins);
    threshold_pwr = mean_pwr * 0.1;
    valid_mask_power = pwr_bins > threshold_pwr;

    dc_guard = 3;
    valid_mask_power(1:dc_guard) = false;
    valid_mask_power(N_fft-dc_guard+1:N_fft) = false;

    % ===== [新增] 将载波分离为 4QAM(导频) 和 16QAM(数据) =====
    % 假设 488~495 和 528~535 是 MATLAB 的 1-based 索引编号
    % 如果你的编号是 0-based，可以相应改成 489:496 和 529:536
    pilot_indices = [488:495, 528:535]; 
    is_pilot = false(N_fft, 1);
    is_pilot(pilot_indices) = true;

    valid_data_mask = valid_mask_power & ~is_pilot;
    valid_pilot_mask = valid_mask_power & is_pilot;

    data_idxs = find(valid_data_mask);
    pilot_idxs = find(valid_pilot_mask);

    if isempty(data_idxs)
        warning('Offset %d 无有效的数据子载波，跳过。', offset);
        continue;
    end

    rel_freq = zeros(N_fft, 1);
    rel_freq(1:N_fft/2) = (0:N_fft/2 - 1).';
    rel_freq(N_fft/2+1:N_fft) = (-N_fft/2:-1).';

    % 提取 16QAM 数据载波
    freq_indices_data = rel_freq(data_idxs);
    syms_valid_data = x_sss_freq(data_idxs);

    [freq_indices_data, sort_idx] = sort(freq_indices_data);
    syms_valid_data = syms_valid_data(sort_idx);
    data_idxs = data_idxs(sort_idx);

    % 提取 4QAM 导频载波
    freq_indices_pilot = rel_freq(pilot_idxs);
    syms_valid_pilot = x_sss_freq(pilot_idxs);

    [freq_indices_pilot, sort_idx_p] = sort(freq_indices_pilot);
    syms_valid_pilot = syms_valid_pilot(sort_idx_p);

    full_freq_indices = [(0:N_fft/2-1), (-N_fft/2:-1)].';

    % ===== 4QAM 盲补偿：网格搜索“线性相位斜率 + 常量相位” =====
    % 注意：此时转为仅将 16个 "4QAM导频子载波" 送入针对4QAM的盲搜评估函数
    [x_sss_freq_corr, slope_hat, phi_hat, cost_best] = blind_phase_slope_comp_4qam( ...
        x_sss_freq, syms_valid_pilot, freq_indices_pilot, full_freq_indices);

    if demod_all_subcarriers
        syms_payload_data = x_sss_freq_corr(:);
        payload_label = sprintf('All carriers: %d', N_fft);
    else
        syms_payload_data = x_sss_freq_corr(data_idxs);
        payload_label = sprintf('Data carriers: %d', length(data_idxs));
    end

    syms_payload_pilot = x_sss_freq_corr(pilot_idxs);

    % 仅为了可视化把幅度规范到标准16QAM平均功率（Es=10）
    % 对 4QAM 也可以用相同的增益去乘，方便在同一个图里观测
    Es_target = 10;
    g = sqrt(Es_target / (mean(abs(syms_payload_data).^2) + eps));
    
    syms_plot_data = syms_payload_data * g;
    syms_plot_pilot = syms_payload_pilot * g;

    figure(fig_const); clf;
    % 绘制 16QAM 数据点 (蓝色大圆点)
    plot(real(syms_plot_data), imag(syms_plot_data), 'b.', 'MarkerSize', 9, 'DisplayName', '16QAM Data'); hold on; grid on; axis square;
    
    % 绘制 4QAM 导频点 (绿色星号)
    if ~isempty(syms_plot_pilot)
        plot(real(syms_plot_pilot), imag(syms_plot_pilot), 'g*', 'MarkerSize', 6, 'LineWidth', 1.0, 'DisplayName', '4QAM Pilot');
    end

    % 叠加理想16QAM参考点
    lv = [-3 -1 1 3];
    [Iref, Qref] = meshgrid(lv, lv);
    plot(Iref(:), Qref(:), 'ro', 'MarkerSize', 5, 'LineWidth', 1.0);

    xlabel('I'); ylabel('Q');
    title(sprintf(['16QAM Constellation (Offset=%d)\\n', ...
        'phase slope=%.5f rad/bin, phi=%.3f rad, cost=%.4f | %s'], ...
        offset, slope_hat, phi_hat, cost_best, payload_label));
    xlim([-4.2 4.2]); ylim([-4.2 4.2]);

    fprintf('Offset=%d | slope=%.6f rad/bin | phi=%.4f rad | cost=%.5f | %s\n', ...
        offset, slope_hat, phi_hat, cost_best, payload_label);

    % 记录最佳（成本最小）的 offset
    if cost_best < best_cost
        best_cost = cost_best;
        best_offset = offset;
        best_slope_hat = slope_hat;
        best_phi_hat = phi_hat;
        best_syms_plot_data = syms_plot_data;
        best_syms_plot_pilot = syms_plot_pilot;
    end

    if pause_each
        fprintf('按 Enter 继续下一个 offset...\n');
        pause;
    else
        drawnow;
    end
end

fprintf('\nSweep 完成。最佳 Offset = %d (此时 Cost 最小: %.5f)\n', best_offset, best_cost);

% 最终再弹出一个高亮显示的单独图像展示最佳星座图
if ~isnan(best_offset)
    figure('Position', [400, 300, 700, 600], 'Name', '16QAM 最佳 Offset 星座图');
    plot(real(best_syms_plot_data), imag(best_syms_plot_data), 'b.', 'MarkerSize', 9, 'DisplayName', '16QAM Data'); hold on; grid on; axis square;
    if ~isempty(best_syms_plot_pilot)
        plot(real(best_syms_plot_pilot), imag(best_syms_plot_pilot), 'g*', 'MarkerSize', 6, 'LineWidth', 1.0, 'DisplayName', '4QAM Pilot');
    end
    [Iref, Qref] = meshgrid([-3 -1 1 3], [-3 -1 1 3]);
    plot(Iref(:), Qref(:), 'ro', 'MarkerSize', 5, 'LineWidth', 1.0);
    xlabel('I'); ylabel('Q');
    title(sprintf('最佳 16QAM 星座图 (Offset=%d)\nphase slope=%.5f, phi=%.3f, cost=%.4f', ...
        best_offset, best_slope_hat, best_phi_hat, best_cost));
    xlim([-4.2 4.2]); ylim([-4.2 4.2]);
end

%% ================= 本地函数区 =================

%% ===== 自动中心频率与PSS同步相关助手函数 =====
function [f_center_pilot_hz, info] = estimate_center_pilot_blind_from_signal( ...
    x, Fs, target_bw_hz, bw_tol_ratio, occ_bg_win_hz, occ_smooth_hz, occ_thresh_db, min_component_bw_hz, ...
    max_expected_cfo_hz, pilot_bg_win_hz, pilot_min_prom_db, pilot_width_ref_hz, ...
    refine_enable, refine_lpf_bw_hz, refine_iters, refine_fir_order, refine_max_delta_hz)

segLen = 32768;
overlapRatio = 0.75;
nfft = 65536;

N = length(x);
if N < 2048
    error('样本过短。');
end
if N < segLen * 2
    segLen = 2^floor(log2(max(1024, floor(N/2))));
    segLen = min(segLen, N);
    nfft = max(nfft, 2^nextpow2(segLen));
end

[P_med, ~, f] = robust_welch_psd_local(x, Fs, segLen, overlapRatio, nfft);
Pdb = 20*log10(P_med + eps);

bin_hz = Fs / length(f);
occ_bg_bins = make_odd_local(max(9, round(occ_bg_win_hz / bin_hz)));
occ_sm_bins = make_odd_local(max(5, round(occ_smooth_hz / bin_hz)));

P_occ_base = movmedian(Pdb, occ_bg_bins);
P_occ_env = movmean(P_occ_base, occ_sm_bins);
noise_floor_db = prctile(P_occ_env, 20);
occ_thr_db = noise_floor_db + occ_thresh_db;
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
    oL = max(1, iL-guard):max(1, iL-1);
    oR = min(length(P_occ_env), iR+1):min(length(P_occ_env), iR+guard);
    if isempty(oL), outL = noise_floor_db; else, outL = mean(P_occ_env(oL)); end
    if isempty(oR), outR = noise_floor_db; else, outR = mean(P_occ_env(oR)); end
    edge_drop = in_mean - 0.5*(outL + outR);

    score = 3.0*(1 - min(bw_err,1)) + 0.12*(in_mean - noise_floor_db) + 0.18*max(edge_drop,0);
    if bw_err <= bw_tol_ratio
        cand = [cand; iL, iR, score]; %#ok<AGROW>
    end
end

if isempty(cand)
    Nf = length(f);
    bw_bins = max(9, round(target_bw_hz / bin_hz));
    half_bw = floor(bw_bins/2);
    out_half = min(floor(1.5*bw_bins), floor((Nf-1)/2));
    score_sw = -inf(Nf,1);
    for k = (out_half+2):(Nf-out_half-1)
        i1 = k-half_bw; i2 = k+half_bw;
        if i1 < 1 || i2 > Nf, continue; end
        in_band = P_occ_env(i1:i2);
        l1 = max(1, k-out_half); l2 = max(1, i1-1);
        r1 = min(Nf, i2+1);      r2 = min(Nf, k+out_half);
        if l2 < l1 || r2 < r1, continue; end
        out_vals = [P_occ_env(l1:l2); P_occ_env(r1:r2)];
        contrast = mean(in_band) - mean(out_vals);
        score_sw(k) = contrast - 0.15*std(in_band);
    end
    [~, k_best] = max(score_sw);
    f_band_center = f(k_best);
else
    [~, id] = max(cand(:,3));
    iL = cand(id,1); iR = cand(id,2);
    f_band_center = (f(iL) + f(iR)) / 2;
end

search_half_hz = max(1.5e6, max_expected_cfo_hz);
search_half_hz = min(search_half_hz, target_bw_hz/2 - 0.5e6);

[~, f_subbin, pinfo] = detect_center_pilot_once_local( ...
    Pdb, f, f_band_center, search_half_hz, pilot_bg_win_hz, pilot_min_prom_db, pilot_width_ref_hz);

f_center_pilot_hz = f_subbin;

if refine_enable
    f_refined = refine_tone_freq_time_local(x, Fs, f_subbin, refine_lpf_bw_hz, refine_iters, refine_fir_order);
else
    f_refined = f_subbin;
end

refine_delta = f_refined - f_subbin;
if abs(refine_delta) > refine_max_delta_hz
    % 防止时域精修在非理想单音条件下“跑飞”，沿用原脚本这块的保护逻辑
    f_center_pilot_hz = f_subbin;
    refine_used = 0;
else
    f_center_pilot_hz = f_refined;
    refine_used = 1;
end
info = pinfo;
info.f_band_center = f_band_center;
info.occ_threshold_db = occ_thr_db;
info.f_fft_subbin_hz = f_subbin;
info.f_refined_hz = f_refined;
info.refine_delta_hz = refine_delta;
info.refine_used = refine_used;
info.f_axis = f;
info.Pdb = Pdb;
info.P_occ_env = P_occ_env;
info.occ_thr_db = occ_thr_db;
end

function [P_med, P_mean, f] = robust_welch_psd_local(x, Fs, segLen, overlapRatio, nfft)
x = x(:);
N = length(x);
nfftEff = 2^nextpow2(N);
w = ones(N, 1);
X = fftshift(fft(x .* w, nfftEff));
P_out = abs(X); 
P_out = (P_out / max(P_out));
P_med = P_out;
P_mean = P_out;
f = ((-nfftEff/2):(nfftEff/2-1)).' * (Fs / nfftEff);
end

function [f_coarse, f_subbin, info] = detect_center_pilot_once_local(Pdb, f, f_center, search_half_hz, bg_win_hz, prom_min_db, width_ref_hz)
df = abs(f(2)-f(1));
idx = abs(f - f_center) <= search_half_hz;
f_roi = f(idx); 
y = Pdb(idx);
if isempty(y), error('搜索区间为空。'); end
[max_val, ksel] = max(y);
f_coarse = f_roi(ksel);
if ksel <= 1 || ksel >= length(y)
    f_subbin = f_coarse;
else
    y1 = y(ksel-1); y2 = y(ksel); y3 = y(ksel+1);
    den = (y1 - 2*y2 + y3);
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

%% ===== 16QAM解调相关本地函数 =====
function [x_corr, slope_best, phi_best, cost_best] = blind_phase_slope_comp_16qam(x_full, y_valid, k_valid, k_full)
% 在有效子载波上，盲搜索线性相位模型：phase(k)=slope*k+phi
% 代价函数：旋转后到16QAM硬判决点的均方距离

y_valid = y_valid(:);
k_valid = k_valid(:);

if numel(y_valid) < 16
    slope_best = 0;
    phi_best = 0;
    cost_best = inf;
    x_corr = x_full;
    return;
end

% 规范到标准16QAM平均功率，避免幅度尺度对代价函数造成偏置
y_norm = normalize_to_16qam_power(y_valid);

% -------- 粗搜索 --------
slope_grid = linspace(-0.35, 0.35, 141);     % rad/bin
phi_grid   = linspace(-pi, pi, 181);          % rad

[slope_best, phi_best, cost_best] = search_cost(y_norm, k_valid, slope_grid, phi_grid);

% -------- 细搜索 --------
slope_fine = linspace(slope_best - 0.02, slope_best + 0.02, 121);
phi_fine   = linspace(phi_best - pi/12, phi_best + pi/12, 121);

[slope_best, phi_best, cost_best] = search_cost(y_norm, k_valid, slope_fine, phi_fine);

% 应用到全频点
phase_full = slope_best * k_full + phi_best;
x_corr = x_full .* exp(-1j * phase_full);
end

function [slope_best, phi_best, cost_best] = search_cost(y_norm, k_valid, slope_grid, phi_grid)
cost_best = inf;
slope_best = 0;
phi_best = 0;

for is = 1:length(slope_grid)
    s = slope_grid(is);
    y_s = y_norm .* exp(-1j * s * k_valid);

    for ip = 1:length(phi_grid)
        p = phi_grid(ip);
        z = y_s * exp(-1j * p);

        z_hat = slicer_16qam(z);
        d = z - z_hat;

        % 鲁棒代价：保留90%较小误差，降低离群点影响
        e2 = abs(d).^2;
        thr = prctile(e2, 90);
        e2 = e2(e2 <= thr);
        c = mean(e2);

        if c < cost_best
            cost_best = c;
            slope_best = s;
            phi_best = p;
        end
    end
end

% 相位归一化到 [-pi, pi]
phi_best = mod(phi_best + pi, 2*pi) - pi;
end

function z_hat = slicer_16qam(z)
% 硬判决到标准16QAM点集（I/Q 电平：±1, ±3）
I = real(z);
Q = imag(z);

Ih = quant_level_4pam(I);
Qh = quant_level_4pam(Q);

z_hat = Ih + 1j*Qh;
end

function q = quant_level_4pam(x)
% 4-PAM 硬判决门限：[-inf,-2),[-2,0),[0,2),[2,inf)
q = zeros(size(x));
q(x < -2) = -3;
q(x >= -2 & x < 0) = -1;
q(x >= 0 & x < 2) = 1;
q(x >= 2) = 3;
end

function y = normalize_to_16qam_power(x)
% 标准16QAM平均符号功率 Es = 10（未归一化点集±1/±3）
Es_target = 10;
g = sqrt(Es_target / (mean(abs(x).^2) + eps));
y = x * g;
end

%% ===== 4QAM 盲解调核心函数 =====
function [x_corr, slope_best, phi_best, cost_best] = blind_phase_slope_comp_4qam(x_full, y_valid, k_valid, k_full)
% 仅在 4QAM 导频子载波上运行二维盲搜，门限基于 ±1
y_valid = y_valid(:);
k_valid = k_valid(:);

if numel(y_valid) < 8
    warning('4QAM 导频数量不足，跳过相位补偿。');
    slope_best = 0;
    phi_best = 0;
    cost_best = inf;
    x_corr = x_full;
    return;
end

% 规范化到 4QAM (Es=2)
y_norm = normalize_to_4qam_power(y_valid);

% -------- 粗搜索 --------
slope_grid = linspace(-0.35, 0.35, 141);     % rad/bin
phi_grid   = linspace(-pi, pi, 181);          % rad

[slope_best, phi_best, cost_best] = search_cost_4qam(y_norm, k_valid, slope_grid, phi_grid);

% -------- 细搜索 --------
slope_fine = linspace(slope_best - 0.02, slope_best + 0.02, 121);
phi_fine   = linspace(phi_best - pi/12, phi_best + pi/12, 121);

[slope_best, phi_best, cost_best] = search_cost_4qam(y_norm, k_valid, slope_fine, phi_fine);

% 应用基于 4QAM 估算出的最优斜率和常量相位到全频点
phase_full = slope_best * k_full + phi_best;
x_corr = x_full .* exp(-1j * phase_full);
end

function [slope_best, phi_best, cost_best] = search_cost_4qam(y_norm, k_valid, slope_grid, phi_grid)
cost_best = inf;
slope_best = 0;
phi_best = 0;

for is = 1:length(slope_grid)
    s = slope_grid(is);
    y_s = y_norm .* exp(-1j * s * k_valid);

    for ip = 1:length(phi_grid)
        p = phi_grid(ip);
        z = y_s * exp(-1j * p);

        z_hat = slicer_4qam(z);
        d = z - z_hat;

        % 因为点数很少(16个)，我们依然去掉最差的那个以抗突发脉冲干扰
        e2 = abs(d).^2;
        thr = prctile(e2, 90);
        e2 = e2(e2 <= thr);
        c = mean(e2);

        if c < cost_best
            cost_best = c;
            slope_best = s;
            phi_best = p;
        end
    end
end

phi_best = mod(phi_best + pi, 2*pi) - pi;
end

function z_hat = slicer_4qam(z)
% 硬判决到标准 4QAM 点集（I/Q 电平仅有：±1）
I = sign(real(z)); I(I==0) = 1;
Q = sign(imag(z)); Q(Q==0) = 1;
z_hat = I + 1j*Q;
end

function y = normalize_to_4qam_power(x)
% 标准 4QAM 平均符号功率 Es = 2（未归一化点集±1）
Es_target = 2;
g = sqrt(Es_target / (mean(abs(x).^2) + eps));
y = x * g;
end
