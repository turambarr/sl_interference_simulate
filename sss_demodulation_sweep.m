% sss_demodulation_sweep.m
% SSS 偏移量遍历试验脚本（仅用于 sweep/观察，不用于最终单点导出）

clear; clc; close all;

%% 1. 参数设置
inFile = 'sigtest1.iq';

% ===== 手动填写的 3 个关键参数 =====
read_start_sample = 16386-874*2; % [参数1] 文件中开始读取的点（原始采样点）
read_length = 6992*3+1000;                  % [参数2] 读取长度（原始采样点个数）
sss_decode_start_idx = 1048;           % [参数3] 从重采样后序列 x_sro 的第几个点开始作为 SSS 解调基准

fs_source = 409.6e6;
fs_target = 60e6;

sro_ppm = 0;
cfo_hz = 17556;

N_fft = 1024;
freq_shift_hz = 63e6;

offset_range = 0:5;   % 遍历范围
pause_each = true;      % 每个 offset 显示后是否暂停

% 解调模式：true=解调全部1024子载波；false=仅解调有效子载波
demod_all_subcarriers = true;

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

%% 5. 遍历 offset
fprintf('Analyzing Time Domain Signal...\n');
base_start_idx = sss_decode_start_idx;

fig_const = figure('Position', [300, 300, 700, 600], 'Name', 'SSS Demodulation Sweep');

for offset = offset_range
    sss_start_idx_60 = base_start_idx + offset;

    if sss_start_idx_60 < 1 || sss_start_idx_60 + N_fft > length(x_sro)
        warning('Offset %d 导致下标越界，跳过。', offset);
        continue;
    end

    x_sss_time = x_sro(sss_start_idx_60 : sss_start_idx_60 + N_fft - 1);
    x_sss_freq = fft(x_sss_time, N_fft) / sqrt(N_fft);

    pwr_bins = abs(x_sss_freq).^2;
    mean_pwr = mean(pwr_bins);
    threshold_pwr = mean_pwr * 0.1;
    valid_mask_power = pwr_bins > threshold_pwr;

    dc_guard = 3;
    valid_mask_power(1:dc_guard) = false;
    valid_mask_power(N_fft-dc_guard+1:N_fft) = false;

    valid_idxs = find(valid_mask_power);

    rel_freq = zeros(N_fft, 1);
    rel_freq(1:N_fft/2) = (0:N_fft/2 - 1).';
    rel_freq(N_fft/2+1:N_fft) = (-N_fft/2:-1).';

    freq_indices = rel_freq(valid_idxs);
    syms_valid = x_sss_freq(valid_idxs);

    [freq_indices, sort_idx] = sort(freq_indices);
    syms_valid = syms_valid(sort_idx);
    valid_idxs = valid_idxs(sort_idx);

    % ===== 抗 DC 缺口相位估计：分段 unwrap + 分段 polyfit（更稳） =====
    syms_pow4 = syms_valid.^4;
    df = diff(freq_indices);

    % 以 df>1 的位置作为频率断点（例如跨越 DC 缺口）
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

        slope_list(end+1,1) = p_seg(1); %#ok<SAGROW>
        w_list(end+1,1) = numel(idx_seg); %#ok<SAGROW>
    end

    if isempty(slope_list)
        warning('Offset %d 连续有效载波不足，跳过。', offset);
        continue;
    end

    % 用段长加权平均估计 4 次方域斜率
    slope4 = sum(slope_list .* w_list) / sum(w_list);

    % 已知 slope4 后，用复均值估计 4 次方域常量相位项（对噪声更稳）
    pow4_detrended = syms_pow4 .* exp(-1j * slope4 * freq_indices);
    phase0_4 = angle(mean(pow4_detrended));

    full_freq_indices = [(0:N_fft/2-1), (-N_fft/2:-1)].';
    phase_correction = (slope4 * full_freq_indices + phase0_4) / 4;
    x_sss_freq_corr = x_sss_freq .* exp(-1j * phase_correction);
    syms_all = x_sss_freq_corr(:);

    subplot(1, 2, 1);
    cla;
    plot(real(syms_all), imag(syms_all), 'b.', 'MarkerSize', 8);
    hold on; grid on; axis square;
    th = 0:0.01:2*pi;
    plot(cos(th), sin(th), 'k:', 'LineWidth', 0.5);
    title(sprintf('SSS Constellation\nOffset=%d, STO slope=%.4f', offset, slope4/4));
    xlabel('I'); ylabel('Q');
    xlim([-2 2]); ylim([-2 2]);

    if demod_all_subcarriers
        syms_payload = x_sss_freq_corr(:);  % 全部 1024 子载波都参与硬判决
        payload_label = sprintf('All carriers: %d', N_fft);
    else
        syms_payload = x_sss_freq_corr(valid_idxs);
        payload_label = sprintf('Valid carriers: %d', length(valid_idxs));
    end
    bits_I = real(syms_payload) < 0;
    bits_Q = imag(syms_payload) < 0;

    demod_bits = zeros(length(syms_payload)*2, 1);
    demod_bits(1:2:end) = bits_I;
    demod_bits(2:2:end) = bits_Q;

    full_hex_str = '';
    for i_hex = 1:4:length(demod_bits)
        if i_hex+3 > length(demod_bits)
            chunk = demod_bits(i_hex:end);
            val = 0;
            for b = 1:length(chunk)
                val = val + chunk(b) * 2^(length(chunk)-b);
            end
            full_hex_str = [full_hex_str, dec2hex(val)];
            break;
        end
        chunk = demod_bits(i_hex:i_hex+3);
        val = chunk(1)*8 + chunk(2)*4 + chunk(3)*2 + chunk(4)*1;
        full_hex_str = [full_hex_str, dec2hex(val)];
    end

    hex_snippet = full_hex_str(1:min(32, length(full_hex_str)));
    subplot(1, 2, 2);
    cla; axis off;
    title('Demodulated Hard Bits (Snippet)');
    text(0.08, 0.70, sprintf('Offset = %d\nFirst 32 hex chars:\n%s...', offset, hex_snippet), 'FontSize', 12, 'Interpreter', 'none');
    text(0.08, 0.42, sprintf('Total bits: %d\n%s * 2 bit/sym', length(demod_bits), payload_label), 'FontSize', 12);

    fprintf('Offset=%d | %s | first 32 hex=%s\n', offset, payload_label, hex_snippet);

    if pause_each
        fprintf('按 Enter 继续下一个 offset...\n');
        pause;
    else
        drawnow;
    end
end

fprintf('\nSweep 完成。\n');


