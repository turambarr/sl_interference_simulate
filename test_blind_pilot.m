% test_blind_pilot.m
% 独立测试自动盲估中心单载波频率的脚本

clear; clc; close all;

%% 1. 参数设置
inFile = 'sigtest1.iq'; % 或 '8e6.iq'
input_format = 'auto';
dat_header_bytes = 0; 
sync_scan_start = 0;      
sync_scan_len   = 5e6;    % 盲估取 5M 个点一般足够了
fs_source = 409.6e6;
center_nominal_hz = 63.5e6;

% 盲估模式参数
blind_target_bw_hz = 60e6;
blind_bw_tol_ratio = 0.35;
blind_occ_bg_win_hz = 2.0e6;
blind_occ_smooth_hz = 1.2e6;
blind_occ_thresh_db = 12.0; % [修改] 从5.5提高到12，确保能切出信号的高原台阶，避开整个底噪带
blind_min_component_bw_hz = 20e6;
blind_max_expected_cfo_hz = 6e6;
blind_pilot_bg_win_hz = 0.9e6;
blind_pilot_min_prom_db = 4.5;
blind_pilot_width_ref_hz = 0.25e6;

% 盲估后时域精修参数
blind_refine_enable = false;
blind_refine_lpf_bw_hz = 0.6e6;   
blind_refine_iters = 2;           
blind_refine_fir_order = 256;     
blind_refine_max_delta_hz = 120e3; 

%% 2. 加载数据
fprintf('读取数据: %s...\n', inFile);
[x_sync, meta_sync] = read_iq_auto_local(inFile, sync_scan_start, sync_scan_len, input_format, dat_header_bytes);
if meta_sync.numSamplesRead <= 0
    error('未读取到有效数据。');
end
x_sync = double(x_sync(:));
x_sync = x_sync - mean(x_sync);
x_sync = x_sync / (mean(abs(x_sync)) + eps);

%% 3. 执行盲估
fprintf('开始盲估中心频率...\n');
[f_center_pilot_hz, blind_info] = estimate_center_pilot_blind_from_signal( ...
    x_sync, fs_source, blind_target_bw_hz, blind_bw_tol_ratio, ...
    blind_occ_bg_win_hz, blind_occ_smooth_hz, blind_occ_thresh_db, blind_min_component_bw_hz, ...
    blind_max_expected_cfo_hz, blind_pilot_bg_win_hz, blind_pilot_min_prom_db, blind_pilot_width_ref_hz, ...
    blind_refine_enable, blind_refine_lpf_bw_hz, blind_refine_iters, blind_refine_fir_order, ...
    blind_refine_max_delta_hz);

freq_shift_hz_used = f_center_pilot_hz;

fprintf('\n>>> 盲估结果 <<<\n');
fprintf('Freq mode: BLIND_PILOT, center=%.6f MHz (equiv CFO vs nominal=%.3f kHz)\n', ...
    freq_shift_hz_used/1e6, (freq_shift_hz_used-center_nominal_hz)/1e3);

%% 4. 画图查看盲估推算过程
fprintf('绘制推算频谱图...\n');
figure('Name', 'Blind Pilot Center Estimation', 'Position', [100, 100, 900, 500]);
f_mhz = blind_info.f_axis / 1e6;

plot(f_mhz, blind_info.Pdb, 'Color', [0.7 0.7 0.7], 'LineWidth', 1, 'DisplayName', '原始 PSD (Welch)'); hold on; grid on;
plot(f_mhz, blind_info.P_occ_env, 'b', 'LineWidth', 1.5, 'DisplayName', '平滑包络 (P\_occ\_env)');
yline(blind_info.occ_threshold_db, 'g--', 'LineWidth', 1.5, 'DisplayName', '占用判定阈值 (occ\_thr\_db)');

% 标注宽带粗估中心
xline(blind_info.f_band_center / 1e6, 'c-', 'LineWidth', 1.5, 'DisplayName', '粗估带宽中心');

% 标注最终定位的单载波峰 (Pilot Center)
plot(f_center_pilot_hz / 1e6, blind_info.P_occ_env(find(blind_info.f_axis >= f_center_pilot_hz, 1)), ...
    'r*', 'MarkerSize', 12, 'LineWidth', 1.5, 'DisplayName', sprintf('精确定位中心 (%.3f MHz)', f_center_pilot_hz/1e6));

xlabel('频率 (MHz)');
ylabel('功率谱密度 PSD (dB)');
title(sprintf('盲估“中心单载波频率”过程全景图\n(估计中心 = %.6f MHz, 相对偏移 = %.3f kHz)', ...
    f_center_pilot_hz/1e6, (f_center_pilot_hz - center_nominal_hz)/1e3));
legend('Location', 'northeast');
xlim([min(f_mhz) max(f_mhz)]);
ylim([min(blind_info.Pdb)-5, max(blind_info.Pdb)+10]);

%% ================= 本地函数区 =================

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
% 注意：此时 P_med 已经是对齐 plot_psd.m 的单量纲振幅比值了
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
    if bw < min_component_bw_hz
        continue;
    end

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
    % 滑窗兜底
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
% 将核心观测数据打包，供外层画图分析
info.f_axis = f;
info.Pdb = Pdb;
info.P_occ_env = P_occ_env;
info.occ_thr_db = occ_thr_db;
end

function [P_med, P_mean, f] = robust_welch_psd_local(x, Fs, segLen, overlapRatio, nfft)
% 完全抛弃 Welch 分段 PSD，直接1:1复刻 plot_psd.m 暴力大点数FFT，赚取单音极大的相干增益
x = x(:);
N = length(x);

% 直接取 2^18 或全长（看齐 plot_psd.m 的 nfftEff = 2^nextpow2(N)）
nfftEff = 2^nextpow2(N);
w = ones(N, 1);

X = fftshift(fft(x .* w, nfftEff));
% 直接用振幅幅度计算dB（不除以 Fs，刻意让单音靠着 N 的成倍暴涨拔高），与 plot_psd.m 完全一致
P_out = abs(X); 
P_out = (P_out / max(P_out)); % 归一化，使得那根最强的针顶格

% 强制把这全尺寸频谱输送出去
P_med = P_out;
P_mean = P_out;
f = ((-nfftEff/2):(nfftEff/2-1)).' * (Fs / nfftEff);
end

function [f_coarse, f_subbin, info] = detect_center_pilot_once_local(Pdb, f, f_center, search_half_hz, bg_win_hz, prom_min_db, width_ref_hz)
df = abs(f(2)-f(1));
idx = abs(f - f_center) <= search_half_hz;
f_roi = f(idx); 
y = Pdb(idx);

if isempty(y)
    error('搜索区间为空。');
end

% 【核心修改】：去掉了所有鲁棒性的多维度记分（抛弃对比背景底噪和峰宽），
% 直接暴力寻找指定带宽窗口内，PSD 幅值最高的那一根单音主峰！
[max_val, ksel] = max(y);
f_coarse = f_roi(ksel);

% 保留三点抛物线插值，用于在离散频率格点之间再提取一步小数精确度
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

function [x, meta] = read_iq_auto_local(filename, startSample, numSamples, input_format, dat_header_bytes)
if nargin < 5 || isempty(dat_header_bytes), dat_header_bytes = 0; end
if nargin < 4 || isempty(input_format), input_format = 'auto'; end
[~, ~, ext] = fileparts(filename); ext = lower(ext); fmt = lower(input_format);
if strcmp(fmt, 'auto')
    if strcmp(ext, '.dat'), fmt = 'dat'; else, fmt = 'iq'; end
end
if strcmp(fmt, 'iq'), [x, meta] = iq_read_int16_le(filename, startSample, numSamples);
elseif strcmp(fmt, 'dat'), [x, meta] = iq_read_int16_le(filename, startSample, numSamples, dat_header_bytes);
else, error('不支持的 input_format'); end
end