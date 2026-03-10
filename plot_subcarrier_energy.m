% plot_subcarrier_energy.m
% 读出划定范围的信号，进行降采样(Farrow)并对齐提取后，
% 绘制 1024 个 OFDM 子载波的能量分布，找出有效子载波位置。

clear; clc; close all;

%% 1. 参数设置
inFile = 'sigtest1.iq';

% 你请求的：设定提取处理的 起点 和 终点 (或长度)
% 这里的数值单位是 原始采样率 409.6MHz 的采样点，从 0 开始计数
extract_start_sample = 21876+6992;              % 处理的起点标号 (比如改成 16386)
extract_end_sample   = 28789+6992;         % 处理的终点标号 (设为 Inf 可以读到文件尾)

fs_source = 409.6e6;
fs_target = 60e6;

sro_ppm  = 0;       
cfo_hz   = 3079075.69;    

N_fft = 1024;
target_offset = -4; 
freq_shift_hz = 63e6; % 频谱向左搬移 63MHz 以归基带

%% 2. 读取原始数据并绘制提取范围的时域图
fprintf('Loading data for subcarrier check...\n');
d = dir(inFile);
full_len = floor(d.bytes / 4);

% 安全限制终点
if isinf(extract_end_sample) || extract_end_sample >= full_len
    extract_end_sample = full_len - 1;
end
read_len = extract_end_sample - extract_start_sample + 1;

% 正式读取设定的目标片段用于分析
[x_raw, ~] = iq_read_int16_le(inFile, extract_start_sample, read_len);
x_raw = double(x_raw);

figure('Name', 'Extracted Range Time Domain', 'Position', [100, 600, 1000, 300]);
plot(extract_start_sample : extract_end_sample, abs(x_raw));
title(sprintf('提取范围波形: %d ~ %d', extract_start_sample, extract_end_sample));
xlabel('Sample Index (409.6MHz)'); ylabel('Amplitude');

x_raw = x_raw - mean(x_raw);
x_raw = x_raw / mean(abs(x_raw));

%% 3. 频谱下变频与 CFO 修正 (分步独立处理)
t_vec = (0:length(x_raw)-1) / fs_source;

% 第一步：在 409.6MHz 原始采样率下，向左搬移 63MHz (DDC 归基带)
fprintf('Shifting spectrum left by %.1f MHz at 409.6MHz...\n', freq_shift_hz/1e6);
x_shifted = x_raw .* exp(-1j * 2 * pi * freq_shift_hz * t_vec).'; 

% 第二步：单独进行 CFO 频偏补偿
fprintf('Applying separately CFO Correction: %.2f Hz...\n', cfo_hz);
x_cfo = x_shifted .* exp(-1j * 2 * pi * cfo_hz * t_vec).'; 

%% 4. SRO 修正与下采样 (Farrow)
fs_eff = fs_source * (1 + sro_ppm/1e6);

% 低通滤波抗混叠
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

%% 5. 对分块提取 1024 点 FFT 并求平均能量
% 将降采样后的整段 60MHz 信号切成 1024 段的 Block
M = floor(length(x_sro) / N_fft);
x_sro_trunc = x_sro(1 : M * N_fft);
x_sro_matrix = reshape(x_sro_trunc, N_fft, M);

% 对所有符号进行并行的 1024 点 FFT
x_sss_freq_all = fft(x_sro_matrix, N_fft, 1) / sqrt(N_fft);  

% 在时间刻度上求各频点能量的平均值
mean_power = mean(abs(x_sss_freq_all).^2, 2);
subcarrier_power_dB = 10 * log10(mean_power + 1e-12);

%% 6. 绘制子载波能量分布
figure('Name', 'Downsampled Subcarrier Energy', 'Position', [150, 150, 900, 450]);
plot(1:N_fft, subcarrier_power_dB, 'b-', 'LineWidth', 1.5, 'DisplayName', 'FFT Bins Power');
hold on; grid on;

max_power = max(subcarrier_power_dB);
yline(max_power - 15, 'r--', 'Threshold (Signal vs Noise)', 'DisplayName', '-15dB Threshold');

% 用色块标出预期代码中写死的正负有效载波区间 (正半轴 2:413, 负半轴 613:1024)
x_fill1 = [2 413 413 2];
x_fill2 = [613 1024 1024 613];
y_fill = [max_power-60 max_power-60 max_power+10 max_power+10];

patch(x_fill1, y_fill, 'm', 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'DisplayName', 'Expected Pos Active (2:413)');
patch(x_fill2, y_fill, 'g', 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'DisplayName', 'Expected Neg Active (613:1024)');

title(sprintf('指定范围(%d~%d)降采样后 %d 点 OFDM 平均子载波能量分布', extract_start_sample, extract_end_sample, N_fft));
xlabel('MATLAB FFT Index (1 to 1024)');
ylabel('Subcarrier Power (dB)');
xlim([1 1024]);
ylim([max_power - 50, max_power + 5]);
legend('Location', 'south');

fprintf('分析完成，请查看图像。通过观察能量凸起区域，可以校验代码中 412 的载波个数是否准确匹配。\n');