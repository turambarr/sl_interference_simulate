% sss_demodulation.m
% SSS 符号解调脚本 (OFDM + 4QAM)
% 基于 GARDNER 模块估计出的 SRO 和 CFO 进行预修正

clear; clc; close all;

%% 1. 参数设置
inFile = 'sigtest8.iq';
startSample = 23003; 
readLen = 10000;     % 稍微多读一点，防止重采样后点数不足 1024

% 原始采样率
fs_source = 409.6e6;
% 目标符号/采样率 (OFDM基带)
fs_target = 60e6;

% --- 之前估计的同步参数 ---
% 请根据 gardner_farrow_timing_recovery.m 的输出回填这里
% 示例值：
sro_ppm  = 0;       % 采样率偏差 (ppm)
cfo_hz   = 3079075.69;    % 载波频偏 (Hz)

% SSS 长度: 1个 OFDM 符号 = 1024 点 (假设无 CP 或与 pss_test.m 一致)
N_fft = 1024;

%% 2. 读取原始数据
fprintf('Loading file: %s from %d, len %d...\n', inFile, startSample, readLen);
[x_raw, ~] = iq_read_int16_le(inFile, startSample, readLen);
x_raw = double(x_raw);

% 归一化
x_raw = x_raw - mean(x_raw);
x_raw = x_raw / mean(abs(x_raw));

% --- 调试：全文件概览 (Visual Check) ---
fprintf('Visualizing read window on whole file context...\n');
try
    d = dir(inFile);
    full_len = d.bytes / 4; % int16 complex
    [x_full_vis, ~] = iq_read_int16_le(inFile, 0, full_len); % 0-based start
    x_full_vis = double(x_full_vis);
    
    figure('Position', [100, 50, 1200, 400], 'Name', 'Whole File Context');
    plot(abs(x_full_vis)); hold on;
    title(['Whole File Amplitude (Start=' num2str(startSample) ', Len=' num2str(readLen) ')']);
    xlabel('Sample Index (1-based in Matlab, 0-based in file+1)'); 
    ylabel('Amplitude');
    
    % 绘制读取区域 (Matlab Index = File Index + 1)
    % Start: startSample + 1
    % End:   startSample + readLen
    high_x = [startSample+1, startSample+readLen, startSample+readLen, startSample+1];
    y_lim = ylim;
    high_y = [y_lim(1), y_lim(1), y_lim(2), y_lim(2)];
    patch(high_x, high_y, 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'r', 'LineWidth', 2);
    
    text(startSample+1, y_lim(2)*0.8, 'Current Read (SSS Analysis)', 'Color', 'r', 'FontWeight', 'bold');
    
    legend('Signal Amplitude', 'Current Read Window');
catch ME
    warning('Visualization failed: %s', ME.message);
end
% -------------------------------------------

%% 3. CFO 修正 (Carrier Frequency Offset Correction)
fprintf('Applying CFO Correction: %.2f Hz...\n', cfo_hz);
t_vec = (0:length(x_raw)-1) / fs_source;
% 补偿项: exp(-j * 2pi * f * t)
x_cfo = x_raw .* exp(-1j * 2 * pi * cfo_hz * t_vec).'; 

%% 4. SRO 修正 (采用 Farrow 分数阶延迟插值消除边界效应)
fprintf('Applying SRO Correction (Farrow Interpolation): %.2f ppm...\n', sro_ppm);
fs_eff = fs_source * (1 + sro_ppm/1e6);

% 使用三次 Lagrange (Cubic Lagrange) Farrow 结构进行重采样
T_in = 1 / fs_eff;
T_out = 1 / fs_target;
% 避免两端由于缺少足量插值基准点而越界，舍弃最后微小的尾部
t_out = 0 : T_out : (length(x_cfo)-3)*T_in; 

% 映射到输入序列的虚拟索引 (1-based)
idx_frac = t_out / T_in + 1; 
idx_base = floor(idx_frac);
mu = idx_frac - idx_base;

% 确保 base 索引不会越界 (Farrow Cubic 需要 base-1 到 base+2 共4个点)
valid_mask = (idx_base >= 2) & (idx_base <= length(x_cfo)-2);
idx_base = idx_base(valid_mask);
mu = mu(valid_mask);

% Farrow 立方插值滤波器系数 (Cubic Lagrange)
h0 = -(mu - 1) .* (mu - 2) .* mu / 6;
h1 =  (mu - 1) .* (mu + 1) .* (mu - 2) / 2;
h2 = -(mu + 1) .* mu .* (mu - 2) / 2;
h3 =  (mu + 1) .* (mu - 1) .* mu / 6;

% 卷积求和（无边界失真效应）
x_sro = h0 .* x_cfo(idx_base - 1) + ...
        h1 .* x_cfo(idx_base) + ...
        h2 .* x_cfo(idx_base + 1) + ...
        h3 .* x_cfo(idx_base + 2);
x_sro = x_sro(:); % 转为列向量

%% 5. 定位 SSS 符号与时域分析
fprintf('Analyzing Time Domain Signal...\n');

% 绘制时域幅度，帮助确认信号位置
figure('Position', [100, 100, 1000, 400], 'Name', 'Time Domain Amplitude');
plot(abs(x_sro)); grid on;
title('Resampled Signal Amplitude (Time Domain)');
xlabel('Sample Index (60MHz)'); ylabel('Amplitude');
xlim([1 length(x_sro)]);

% --- 滑动窗口搜索最佳 FFT 位置 (用户要求暂时去掉时偏补偿/自动搜索) ---
% search_range = 1 : min(3000, length(x_sro) - N_fft);
% ... (Code commented out) ...
% sss_start_idx_60 = best_start_sample;

fprintf('Time Offset Compensation (Auto-Search) DISABLED.\n');
sss_start_idx_60 = 1; % 直接使用起始位置

% 截取 SSS
extract_len = N_fft;
if sss_start_idx_60 + extract_len > length(x_sro)
    error('SSS extraction out of bounds.');
end

x_sss_time = x_sro(sss_start_idx_60 : sss_start_idx_60 + extract_len - 1);

%% 6. OFDM 解调 (FFT)
fprintf('Demodulating SSS (FFT size %d)...\n', N_fft);
x_sss_freq = fft(x_sss_time, N_fft) / sqrt(N_fft);  % 能量归一化 or not

% 提取有效子载波
% 协议明确: 仅 k=0, 1, N-2, N-1 这4个子载波为0，其余均为非零数据
% 对应 MATLAB 里的 1-based 索引为 1, 2, 1023, 1024 为 0
% 因此有效非零子载波的索引就是 3 到 1022 (共1020个点)
idx_nonzero = 3 : (N_fft - 2);
syms_all = x_sss_freq(idx_nonzero); 
syms_all = syms_all(:); % 列向量

%% 7. 绘图与分析
figure('Position', [300, 300, 600, 600], 'Name', 'SSS Demodulation Result');
plot(real(syms_all), imag(syms_all), 'b.', 'MarkerSize', 8);
grid on; axis square;
title(['SSS Constellation (OFDM 4QAM)\nCFO=' num2str(cfo_hz) 'Hz, SRO=' num2str(sro_ppm) 'ppm']);
xlabel('I'); ylabel('Q');
xlim([-2 2]); ylim([-2 2]);

% 绘制参考圆/点
hold on;
th = 0:0.01:2*pi;
plot(cos(th), sin(th), 'k:', 'LineWidth', 0.5); % 单位圆
% 4QAM 理想点 (1+j, -1+j...) / sqrt(2)?
% pss_test 里用的 exp(1j*(...)) 也就是模为1
ref_pts = [1+1j, -1+1j, -1-1j, 1-1j] / sqrt(2); 
% 或者是 exp(1j*pi/4)... 模是1
% gen_ofdm_symbol: exp(1j * (bits * pi/2 + pi/4)) -> Modulus 1
ref_pts_rot = exp(1j * (0:3)*pi/2 + 1j*pi/4);
plot(real(ref_pts_rot), imag(ref_pts_rot), 'rx', 'MarkerSize', 12, 'LineWidth', 2);
legend('Rx Symbols', 'Ideal Reference');

% 自动相位修正 (去除残留相位旋转)
% 就算 CFO 去除了，相位初相可能还是乱的 (Global Phase Rotation)
% 我们可以尝试去旋转
% 找一个使得点最聚拢的角度? 或者简单地看均值(如果是4QAM, E[x^4] = -1)
% % --- 暂时注释掉盲相位补偿的四次方模糊度 --- 
% mean_pow4 = mean(syms_all.^4); 
% est_phase_bias = angle(mean_pow4) / 4; 
% % 4次方会把 pi/4 映射到 pi, 所以...
% % QPSK 4次方都是 -1 (exp(j*pi)=-1)
% % 比如 exp(j*pi/4)^4 = exp(j*pi) = -1
% % 如果有一个偏差 phi, 则 (exp(j(pi/4+phi)))^4 = -1 * exp(j*4phi)
% % angle(mean) = angle(-1 * exp(j4phi)) = pi + 4phi
% % 4phi = angle - pi
% % phi = (angle - pi) / 4
% est_phase_bias = (angle(mean_pow4) - pi) / 4;
% 
% fprintf('Estimated Global Phase Offset: %.2f degrees\n', est_phase_bias * 180/pi);
% 
% % 补偿相位
% syms_all_rot = syms_all * exp(-1j * est_phase_bias);
% 
% figure('Position', [950, 300, 600, 600], 'Name', 'SSS Constellation (Phase Corrected)');
% plot(real(syms_all_rot), imag(syms_all_rot), 'b.', 'MarkerSize', 8);
% grid on; axis square;
% title('SSS Constellation (Global Phase Corrected)');
% xlabel('I'); ylabel('Q');
% xlim([-2 2]); ylim([-2 2]);
% hold on;
% plot(real(ref_pts_rot), imag(ref_pts_rot), 'rx', 'MarkerSize', 12, 'LineWidth', 2);

