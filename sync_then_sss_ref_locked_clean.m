% sync_then_sss_ref_locked_clean.m
% 完整精简版流程：
%   1) 生成本地 PSS 模板，在 409.6 MHz 原始流中做归一化匹配同步；
%   2) 根据 PSS 起点独立读取 SSS 解调流；
%   3) 在 409.6 MHz 域构造 [SSS; SSS; SSS] 干净上下文；
%   4) 用 SSS 局部窗口做鲁棒归一化，并限幅非保护区余量；
%   5) 盲估中心单音，DDC 到基带，32 MHz 低通，Farrow 重采样到 60 MHz；
%   6) 调用 sss_ref_locked_demod_core 做参考辅助补偿、offset 搜索、carrierShift 锁定；
%   7) 输出 HEX 候选与完整 BER 数字。

clear; clc; close all;

%% ===================== 参数区 =====================
inFile = "sigtest1.iq";
input_format = 'auto';       % 'auto' | 'iq' | 'dat'
dat_header_bytes = 0;

fs_source = 409.6e6;
fs_symbol = 60e6;
fs_target = 60e6;
center_nominal_hz = 63.5e6;

sync_scan_start = 0;         % 0-based
sync_scan_len = 20e6;
min_peak_height = 50;        % PSS normalized correlation percent
peak_select_mode = 'first';  % 'first' | 'strongest'

N_fft = 1024;
sss_decode_start_idx = 1024 + 48; % PSS 起点后的 SSS 标称位置，60 MHz 样点
target_offset = 0;

isolated_sss_pre_guard_60 = 0;
isolated_sss_post_guard_60 = 2048;

raw_norm_guard_60 = 256;
raw_margin_clip_factor = 1.0;
raw_margin_protect_guard_60 = 128;

blind_est_pre_start_samples = round(5000 * 409.6 / 600);
blind_est_len = round(40000 * 409.6 / 600);
blind_search_half_hz = 6e6;

lpf_bw_hz = 32e6;
lpf_order = 30;
sro_ppm = 0;

decision_rotate_deg = 45;
fixed_valid_rel_freq_indices = [-512:-10, 7:511];
offset_search = -8:8;
carrier_shift_search = -32:32;
carrier_shift_min_gain = 0.02;
accept_ber_thresh = 0.03;

ber_ref_hex_list = {
    '80C46B5267759FD6505644FC1F0052124876B17F9FA9EA13D3C2A8E274700BE9F511FBCDF1B237B1E6E46FB1A80A762952FCA4FCDB45612E5C6641BC200043442D23B538172067919D89A3D640557B004A5507460C3DC60DEDE45BB18D5FC9D5E944683265F591D6585671A960FF43002E13B478182063925276497C15006A10178064119F89AE831FAABF442C33B9F812205B9272769C29EAFC2E57F1E270347E7388504C853A85EF97AA0E79A3FDD5CC122832DE4EA2E7191A815164CC8CE96EEF4BC01B139D2E5D75B280250380B83D208AC7BD8ACD3A263D2B583218CAA1EB37998E6EF65A2982FFCE5624BB8FCD6F5F4483E3EC5384ED75493BE9CD690B', ...
    '159DC2F4CEEF3ABCF5FCDDA97A55F474D1EC27EA3A038076B6940184EDE55283AF77A29BA7246E278C8DCA270150EC43F4A90DA9B2DFC748F9CCD7294555D6DD4B462F617E45CE373B1306BCD5FFE255D0FF5EDC596B9C5B8B8DF2271BFA93BF83DDC164CFAF37BCF1FCE703C5AAD65548762DE17145C634F4ECD3E97F55C0757E15CD773A1308167A002ADD496623A17445F234E4EC394380A948FEA784E56DE8E611F5D91F601F8A3E0058E306ABBF99744164B8D8048E737017F7CD991983C88AD29572763B48FBEF24154F5615216B45109E2B109B604C6B42F164719007826E3318C8ACF04314AA98FC4D221A9BCAFADD168689F61D8BEFD362839AC352', ...
    '7F3B94AD988A6029AFA9BB03E0FFADEDB7894E80605615EC2C3D571D8B8FF4160AEE04320E4DC84E191B904E57F589D6AD035B0324BA9ED1A399BE43DFFFBCBBD2DC4AC7E8DF986E62765C29BFAA84FFB5AAF8B9F3C239F2121BA44E72A0362A16BB97CD9A0A6E29A7A98E569F00BCFFD1EC4B87E7DF9C6DAD89B683EAFF95EFE87F9BEE6076517CE05540BBD3CC4607EDDFA46D8D8963D61503D1A80E1D8FCB818C77AFB37AC57A106855F1865C022A33EDD7CD21B15D18E6E57EAE9B3373169110B43FE4EC62D1A28A4D7FDAFC7F47C2DF7538427532C5D9C2D4A7CDE7355E14C866719109A5D67D0031A9DB44703290A0BB7C1C13AC7B128AB6C4163096F4', ...
    'EA623D0B3110C5430A03225685AA0B8B2E13D815C5FC7F89496BFE7B121AAD7C50885D6458DB91D8737235D8FEAF13BC0B56F2564D2038B7063328D6BAAA2922B4B9D09E81BA31C8C4ECF9432A001DAA2F00A123A69463A474720DD8E4056C407C223E9B3050C8430E0318FC3A5529AAB789D21E8EBA39CB0B132C1680AA3F8A81EA3288C5ECF7E985FFD522B699DC5E8BBA0DCB1B13C6BC7F56B701587B1A921719EE0A26E09FE075C1FFA71CF95440668BBE9B4727FB718C8FE8083266E67C37752D6A8D89C4B70410DBEAB0A9EADE94BAEF61D4EF649FB394BD0E9B8E6FF87D91CCE737530FBCEB556703B2DDE564350522E9797609E274102C9D7C643CAD'
};
unstable_bit_positions = [];
unstable_nibble_positions = [1, 2, 3, 4, 257, 510, 511, 512];

%% ===================== Step 1: PSS 同步 =====================
fprintf('Step 1: 生成本地 PSS 模板...\n');
pss_local = make_pss_template_local(fs_source, fs_symbol);

fprintf('Step 2: 全局归一化匹配找峰...\n');
file_info = dir(inFile);
if isempty(file_info)
    error('找不到文件: %s', inFile);
end

data_bytes = file_info.bytes;
if strcmpi(input_format, 'dat') || (strcmpi(input_format, 'auto') && endsWith(lower(inFile), '.dat'))
    data_bytes = max(0, data_bytes - dat_header_bytes);
end
total_samples = floor(data_bytes / 4);
sync_scan_len = min(round(sync_scan_len), total_samples - sync_scan_start);

[x_sync, meta_sync] = read_iq_auto_local(inFile, sync_scan_start, sync_scan_len, input_format, dat_header_bytes);
if meta_sync.numSamplesRead <= 0
    error('同步扫描阶段未读取到有效数据。');
end

x_sync = double(x_sync(:));
x_sync = x_sync - mean(x_sync);
x_sync = x_sync / (mean(abs(x_sync)) + eps);

[read_start_sample, sync_info] = find_pss_start_local( ...
    x_sync, pss_local, sync_scan_start, min_peak_height, peak_select_mode);

fprintf('  检测到 %d 个峰。\n', sync_info.num_peaks);
for k = 1:sync_info.num_peaks
    fprintf('    Peak %02d: %.2f%%, 末端=%d, 起始=%d\n', ...
        k, sync_info.pks(k), sync_info.end_pos(k), sync_info.start_pos(k));
end
fprintf('>>> 选中同步起点 read_start_sample = %d\n', read_start_sample);

%% ===================== Step 2: 独立读取 SSS 解调流 =====================
demod_stream_start_60 = max(1, sss_decode_start_idx - isolated_sss_pre_guard_60);
demod_stream_end_60 = sss_decode_start_idx + N_fft - 1 + isolated_sss_post_guard_60;
demod_stream_start_offset_raw = floor((demod_stream_start_60 - 1) * fs_source / fs_target);
demod_stream_len_raw = ceil((demod_stream_end_60 - demod_stream_start_60 + 1) * fs_source / fs_target) + 8;
demod_read_start_sample = read_start_sample + demod_stream_start_offset_raw;
demod_read_length = min(round(demod_stream_len_raw), total_samples - demod_read_start_sample);
sss_decode_start_idx_eff = max(1, sss_decode_start_idx - demod_stream_start_60 + 1);

fprintf('Loading SSS demod stream: start=%d, len=%d, sss_idx_eff=%d\n', ...
    demod_read_start_sample, demod_read_length, sss_decode_start_idx_eff);
[x_raw, ~] = read_iq_auto_local(inFile, demod_read_start_sample, demod_read_length, input_format, dat_header_bytes);
x_raw = double(x_raw(:));
if isempty(x_raw)
    error('SSS 解调阶段读取为空。');
end

%% ===================== Step 3: 409.6 MHz 域构造干净 SSS 上下文 =====================
raw_center_start_60 = sss_decode_start_idx_eff + target_offset;
raw_center_end_60 = raw_center_start_60 + N_fft - 1;
raw_center_start = max(1, floor((raw_center_start_60 - 1) * fs_source / fs_target) + 1);
raw_center_end = min(length(x_raw), ceil((raw_center_end_60 - 1) * fs_source / fs_target) + 1);
if raw_center_end <= raw_center_start
    error('raw SSS 中心窗口越界: raw=%d:%d, len=%d', raw_center_start, raw_center_end, length(x_raw));
end

x_raw_center = x_raw(raw_center_start:raw_center_end);
raw_center_len = length(x_raw_center);
x_raw = [x_raw_center; x_raw_center; x_raw_center];
sss_decode_start_idx_eff = floor(raw_center_len * fs_target / fs_source) + 1;
fprintf('Raw SSS repeat context: center_raw=%d:%d (%d), new_len=%d, sss_idx_eff=%d\n', ...
    raw_center_start, raw_center_end, raw_center_len, length(x_raw), sss_decode_start_idx_eff);

%% ===================== Step 4: 局部鲁棒归一化与余量限幅 =====================
sss_start_idx_60_for_norm = sss_decode_start_idx_eff + target_offset;
norm_start_60 = max(1, sss_start_idx_60_for_norm - raw_norm_guard_60);
norm_end_60 = sss_start_idx_60_for_norm + N_fft - 1 + raw_norm_guard_60;
norm_start_raw = max(1, floor((norm_start_60 - 1) * fs_source / fs_target) + 1);
norm_end_raw = min(length(x_raw), ceil((norm_end_60 - 1) * fs_source / fs_target) + 1);
[x_raw, norm_info] = normalize_iq_local(x_raw, norm_start_raw, norm_end_raw);
fprintf('Raw normalization: win=%d:%d, dc=%.4g%+.4gj, scale=%.4g\n', ...
    norm_info.start_idx, norm_info.end_idx, real(norm_info.dc), imag(norm_info.dc), norm_info.scale);

clip_protect_start_60 = max(1, sss_start_idx_60_for_norm - raw_margin_protect_guard_60);
clip_protect_end_60 = sss_start_idx_60_for_norm + N_fft - 1 + raw_margin_protect_guard_60;
clip_protect_start_raw = max(1, floor((clip_protect_start_60 - 1) * fs_source / fs_target) + 1);
clip_protect_end_raw = min(length(x_raw), ceil((clip_protect_end_60 - 1) * fs_source / fs_target) + 1);
[x_raw, clip_info] = clip_margin_local(x_raw, raw_margin_clip_factor, clip_protect_start_raw, clip_protect_end_raw);
fprintf('Raw margin clip: clip=%.3g, protect=%d:%d, margin=%d, clipped=%d\n', ...
    raw_margin_clip_factor, clip_info.protect_start, clip_info.protect_end, clip_info.n_margin, clip_info.n_clipped);

%% ===================== Step 5: 盲估中心单音 + DDC + LPF + Farrow =====================
blind_read_start = max(0, read_start_sample - blind_est_pre_start_samples);
blind_read_len = min(round(blind_est_len), total_samples - blind_read_start);
[x_blind, meta_blind] = read_iq_auto_local(inFile, blind_read_start, blind_read_len, input_format, dat_header_bytes);
if meta_blind.numSamplesRead < 2048
    error('盲估中心频率读取样本过短。');
end
x_blind = double(x_blind(:));
x_blind = x_blind - mean(x_blind);
x_blind = x_blind / (mean(abs(x_blind)) + eps);
freq_shift_hz = estimate_center_pilot_simple_local(x_blind, fs_source, center_nominal_hz, blind_search_half_hz);
fprintf('Freq estimate: center=%.6f MHz, equiv CFO vs nominal=%.3f kHz\n', ...
    freq_shift_hz/1e6, (freq_shift_hz-center_nominal_hz)/1e3);

t = (0:length(x_raw)-1).' / fs_source;
x_base = x_raw .* exp(-1j * 2*pi * freq_shift_hz * t);

Wn = lpf_bw_hz / (fs_source/2);
b_lpf = fir1(lpf_order, Wn);
x_filt = filtfilt(b_lpf, 1, x_base);
fprintf('Main LPF: filtfilt FIR order=%d, cutoff=%.1f MHz\n', lpf_order, lpf_bw_hz/1e6);

x_sro = farrow_resample_local(x_filt, fs_source * (1 + sro_ppm/1e6), fs_target);
fprintf('Farrow resample: %.1f MHz -> %.1f MHz, len=%d\n', fs_source/1e6, fs_target/1e6, length(x_sro));

%% ===================== Step 6: 参考锁定 SSS 解调 =====================
opts = struct();
opts.N_fft = N_fft;
opts.target_offset = target_offset;
opts.offset_search = offset_search;
opts.decision_rotate_deg = decision_rotate_deg;
opts.fixed_valid_rel_freq_indices = fixed_valid_rel_freq_indices;
opts.ref_hex_list = ber_ref_hex_list;
opts.unstable_bit_positions = unstable_bit_positions;
opts.unstable_nibble_positions = unstable_nibble_positions;
opts.enable_ref_carrier_shift_lock = true;
opts.ref_carrier_shift_search = carrier_shift_search;
opts.ref_carrier_shift_min_gain = carrier_shift_min_gain;
opts.accept_ber_thresh = accept_ber_thresh;
opts.verbose = true;

result = sss_ref_locked_demod_core(x_sro, sss_decode_start_idx_eff, opts);

fprintf('\n>>> [Clean Sync + Ref-Locked SSS Demod] 完成 <<<\n');
fprintf('read_start_sample=%d, SSS offset=%+d, accepted=%d\n', ...
    read_start_sample, result.offset, result.accepted);
fprintf('Compensation=%s\n', result.comp_mode);
fprintf('Carrier shift: selected=%+d, rawBER=%.6g, shiftedBER=%.6g\n', ...
    result.carrier_shift, result.carrier_shift_raw_ber, result.carrier_shift_ber);
fprintf('Best BER: %.6g (err=%d/%d), Ref=%d, Rot=%d deg\n', ...
    result.best_ber, result.best_err, result.best_n, result.best_ref_idx, result.best_rot_deg);

rot_angles_deg = decision_rotate_deg + [0, 90, 180, 270];
fprintf('====== HEX CANDIDATES (QPSK 4-way rotation ambiguity) ======\n');
for r = 1:4
    fprintf('[%3d deg] %s\n', rot_angles_deg(r), result.hex_candidates{r});
end
fprintf('============================================================\n');

print_ber_matrix_local(result.demod_bits_candidates, ber_ref_hex_list, ...
    unstable_bit_positions, unstable_nibble_positions, rot_angles_deg);

figure('Name', 'Clean Ref-Locked SSS Constellation', 'Position', [300, 300, 560, 520]);
syms_plot = result.x_freq_corr(:) .* exp(1j * decision_rotate_deg * pi/180);
plot(real(syms_plot), imag(syms_plot), 'b.', 'MarkerSize', 8); hold on; grid on; axis square;
th = linspace(0, 2*pi, 360);
plot(cos(th), sin(th), 'k:', 'LineWidth', 0.5);
title(sprintf('SSS Constellation | offset=%+d, carrierShift=%+d', result.offset, result.carrier_shift));
xlabel('I'); ylabel('Q'); xlim([-2 2]); ylim([-2 2]);

%% ===================== 本脚本必要辅助函数 =====================
function pss_local = make_pss_template_local(fs_source, fs_symbol)
fc = 63e6;
hex_str = 'BD3CD0148871751F84CED8C1BE32AC96';
bits = zeros(1, 128);
for i = 1:length(hex_str)
    bin_str = dec2bin(hex2dec(hex_str(i)), 4);
    bits((i-1)*4 + 1 : i*4) = bin_str - '0';
end

d_phi = (bits == 1) * (-pi/2) + (bits == 0) * (pi/2);
phase_accum = cumsum(d_phi);
m_syms = exp(1j * phase_accum);
pss_base = [-m_syms, m_syms, m_syms, m_syms, m_syms, m_syms, m_syms, m_syms];
cp_syms = -pss_base(end-48+1:end);
pss_base_with_cp = [cp_syms, pss_base];

[P, Q] = rat(fs_source / fs_symbol);
pss_up = resample(pss_base_with_cp, P, Q);
t_local = (0:length(pss_up)-1) / fs_source;
pss_local = pss_up .* exp(1j * 2*pi * fc * t_local);
pss_local = pss_local / max(abs(pss_local));
end

function [read_start_sample, info] = find_pss_start_local(x_sync, pss_local, sync_scan_start, min_peak_height, peak_select_mode)
h = fliplr(conj(pss_local));
corr_out = fftfilt(h(:), x_sync(:));
corr_mag_sq = abs(corr_out).^2;
E_local = sum(abs(pss_local).^2);
E_rx = fftfilt(ones(length(pss_local),1), abs(x_sync(:)).^2);
E_rx(E_rx < 1e-10) = 1e-10;
corr_norm_pct = (corr_mag_sq ./ (E_local .* E_rx)) * 100;

[pks, locs] = findpeaks(corr_norm_pct, ...
    'MinPeakHeight', min_peak_height, 'MinPeakDistance', length(pss_local));
if isempty(locs)
    error('未检测到 PSS 同步峰。');
end

end_pos = sync_scan_start + (locs - 1);
start_pos = end_pos - (length(pss_local) - 1);
switch lower(peak_select_mode)
    case 'first'
        pick = 1;
    case 'strongest'
        [~, pick] = max(pks);
    otherwise
        error('未知 peak_select_mode: %s', peak_select_mode);
end

read_start_sample = start_pos(pick);
if read_start_sample < 0
    error('选中的 PSS 起始点为负: %d', read_start_sample);
end

info = struct('num_peaks', numel(locs), 'pks', pks, ...
    'end_pos', end_pos, 'start_pos', start_pos, 'picked', pick);
end

function f_center = estimate_center_pilot_simple_local(x, Fs, f_nominal, search_half_hz)
x = x(:);
nfft = 2^nextpow2(length(x));
X = fftshift(fft(x, nfft));
f = ((-nfft/2):(nfft/2-1)).' * (Fs / nfft);
roi = abs(f - f_nominal) <= search_half_hz;
if ~any(roi)
    error('中心单音搜索区间为空。');
end
y = abs(X(roi));
f_roi = f(roi);
[~, k] = max(y);
if k <= 1 || k >= length(y)
    f_center = f_roi(k);
else
    y1 = 20*log10(y(k-1) + eps);
    y2 = 20*log10(y(k) + eps);
    y3 = 20*log10(y(k+1) + eps);
    den = y1 - 2*y2 + y3;
    if abs(den) < 1e-12
        delta = 0;
    else
        delta = max(min(0.5*(y1-y3)/den, 1), -1);
    end
    f_center = f_roi(k) + delta * (Fs / nfft);
end
end

function x_sro = farrow_resample_local(x, fs_in, fs_out)
x = x(:);
T_in = 1 / fs_in;
T_out = 1 / fs_out;
t_out = 0 : T_out : (length(x)-3)*T_in;

idx_frac = t_out / T_in + 1;
idx_base = floor(idx_frac);
mu = idx_frac - idx_base;
valid = (idx_base >= 2) & (idx_base <= length(x)-2);
idx_base = idx_base(valid);
mu = mu(valid);

h0 = -(mu - 1) .* (mu - 2) .* mu / 6;
h1 =  (mu - 1) .* (mu + 1) .* (mu - 2) / 2;
h2 = -(mu + 1) .* mu .* (mu - 2) / 2;
h3 =  (mu + 1) .* (mu - 1) .* mu / 6;

x_sro = h0(:) .* x(idx_base(:)-1) + ...
        h1(:) .* x(idx_base(:)) + ...
        h2(:) .* x(idx_base(:)+1) + ...
        h3(:) .* x(idx_base(:)+2);
x_sro = x_sro(:);
end

function [x_norm, info] = normalize_iq_local(x, win_start, win_end)
x = x(:);
n = length(x);
win_start = max(1, min(n, round(win_start)));
win_end = max(win_start, min(n, round(win_end)));
x_win = x(win_start:win_end);

dc = median(real(x_win)) + 1j * median(imag(x_win));
scale = median(abs(x_win - dc));
if scale < eps
    scale = mean(abs(x_win - dc));
end
if scale < eps
    scale = 1;
end

x_norm = (x - dc) / scale;
info = struct('start_idx', win_start, 'end_idx', win_end, 'dc', dc, 'scale', scale);
end

function [x_out, info] = clip_margin_local(x, clip_level, protect_start, protect_end)
x_out = x(:);
n = length(x_out);
protect_start = max(1, min(n, round(protect_start)));
protect_end = max(protect_start, min(n, round(protect_end)));
margin = true(n, 1);
margin(protect_start:protect_end) = false;

mag = abs(x_out);
idx_margin = find(margin);
clip_mask = mag(idx_margin) > clip_level;
idx_clip = idx_margin(clip_mask);
x_out(idx_clip) = x_out(idx_clip) .* (clip_level ./ (mag(idx_clip) + eps));

info = struct('protect_start', protect_start, 'protect_end', protect_end, ...
    'n_margin', nnz(margin), 'n_clipped', nnz(clip_mask));
end

function [x, meta] = read_iq_auto_local(filename, start_sample, num_samples, input_format, dat_header_bytes)
if nargin < 5
    dat_header_bytes = 0;
end
[~, ~, ext] = fileparts(filename);
fmt = lower(input_format);
if strcmp(fmt, 'auto')
    if strcmpi(ext, '.dat')
        fmt = 'dat';
    else
        fmt = 'iq';
    end
end

switch fmt
    case 'iq'
        [x, meta] = iq_read_int16_le(filename, start_sample, num_samples);
    case 'dat'
        [x, meta] = iq_read_int16_le(filename, start_sample, num_samples, dat_header_bytes);
    otherwise
        error('不支持的 input_format: %s', input_format);
end
end

function print_ber_matrix_local(demod_bits_candidates, ref_hex_list, unstable_bit_positions, unstable_nibble_positions, rot_angles_deg)
fprintf('\n================ BER EVALUATION ================\n');
fprintf('Unstable drop: bit_pos=%d, nibble_pos=%d\n', ...
    numel(unstable_bit_positions), numel(unstable_nibble_positions));

for i = 1:numel(demod_bits_candidates)
    best_ber = nan;
    best_ref = 0;
    bits_est = demod_bits_candidates{i};
    for j = 1:numel(ref_hex_list)
        ref_hex = strtrim(ref_hex_list{j});
        if isempty(ref_hex)
            continue;
        end
        bits_ref = hex_to_bits_local(ref_hex);
        n_use = min(length(bits_est), length(bits_ref));
        keep = make_ber_keep_mask_local(n_use, unstable_bit_positions, unstable_nibble_positions);
        err = nnz((bits_est(1:n_use) ~= bits_ref(1:n_use)) & keep);
        n_keep = nnz(keep);
        ber = err / n_keep;
        fprintf('Cand[%3d deg] vs Ref[%d]: BER=%.6g (err=%d/%d)\n', ...
            rot_angles_deg(i), j, ber, err, n_keep);
        if isnan(best_ber) || ber < best_ber
            best_ber = ber;
            best_ref = j;
        end
    end
    fprintf('Best for Cand[%3d deg] -> Ref[%d], BER=%.6g\n', ...
        rot_angles_deg(i), best_ref, best_ber);
end
fprintf('================================================\n');
end

function keep = make_ber_keep_mask_local(n_use, unstable_bit_positions, unstable_nibble_positions)
keep = true(n_use, 1);
ub = unstable_bit_positions(:);
ub = ub(ub >= 1 & ub <= n_use);
keep(ub) = false;

un = unstable_nibble_positions(:);
un = un(un >= 1);
for k = 1:numel(un)
    b1 = 4*(un(k)-1) + 1;
    b2 = min(4*un(k), n_use);
    if b1 <= n_use
        keep(b1:b2) = false;
    end
end
end

function bits = hex_to_bits_local(hex_str)
hex_str = upper(regexprep(hex_str, '\s+', ''));
hex_str = regexprep(hex_str, '[^0-9A-F]', '');
bits = zeros(4*length(hex_str), 1);
for i = 1:length(hex_str)
    v = hex2dec(hex_str(i));
    bits(4*(i-1)+1) = bitget(v, 4);
    bits(4*(i-1)+2) = bitget(v, 3);
    bits(4*(i-1)+3) = bitget(v, 2);
    bits(4*(i-1)+4) = bitget(v, 1);
end
end
