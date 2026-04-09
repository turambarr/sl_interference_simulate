% sss_demodulation.m
% SSS 符号解调脚本 (OFDM + 4QAM)
% 基于 GARDNER 模块估计出的 SRO 和 CFO 进行预修正

clear; clc; close all;

%% 1. 参数设置
inFile = '20250912222305_part1.iq';
startSample = 307700 - 874 * 8;%16386+874*6+218; 
readLen = 7000 * 8;     % 稍微多读一点，防止重采样后点数不足 1024

% 原始采样率
fs_source = 409.6e6;
% 目标符号/采样率 (OFDM基带)
fs_target = 60e6;

% --- 之前估计的同步参数 ---
% 请根据 gardner_farrow_timing_recovery.m 的输出回填这里
% 示例值：
sro_ppm  = 0;       % 采样率偏差 (ppm)
cfo_hz   = 17556;    % 载波频偏 (Hz)

% SSS 长度: 1个 OFDM 符号 = 1024 点 (假设无 CP 或与 pss_test.m 一致)
N_fft = 1024;
target_offset = -3; % 指定你想要测试的偏移量
freq_shift_hz = 63e6; % 频谱归基带向左搬移 63MHz

%% 2. 读取原始数据
fprintf('Loading file: %s from %d, len %d...\n', inFile, startSample, readLen);
[x_raw, ~] = iq_read_int16_le(inFile, startSample, readLen);
x_raw = double(x_raw);

% 归一化
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

%% 4. SRO 修正 (采用 Farrow 分数阶延迟插值消除边界效应)
fprintf('Applying SRO Correction (Farrow Interpolation): %.2f ppm...\n', sro_ppm);
fs_eff = fs_source * (1 + sro_ppm/1e6);

% --- 新增：使用零相移(Zero-Phase)低通滤波器抗混叠 ---
% 目标下降到 60MHz 采样率，放宽截止频率为 35MHz
fprintf('Applying Zero-Phase Anti-aliasing LPF (35MHz)...\n');
Wn = 35e6 / (fs_source / 2);
lpf_order = 30;
b_lpf = fir1(lpf_order, Wn);
% 消除群时延(相位偏移)，过滤 CFO 修正后的信号
x_cfo_filtered = filtfilt(b_lpf, 1, x_cfo);

% 使用三次 Lagrange (Cubic Lagrange) Farrow 结构进行重采样
T_in = 1 / fs_eff;
T_out = 1 / fs_target;
% 避免两端由于缺少足量插值基准点而越界，舍弃最后微小的尾部
t_out = 0 : T_out : (length(x_cfo_filtered)-3)*T_in; 

% 映射到输入序列的虚拟索引 (1-based)
idx_frac = t_out / T_in + 1; 
idx_base = floor(idx_frac);
mu = idx_frac - idx_base;

% 确保 base 索引不会越界 (Farrow Cubic 需要 base-1 到 base+2 共4个点)
valid_mask = (idx_base >= 2) & (idx_base <= length(x_cfo_filtered)-2);
idx_base = idx_base(valid_mask);
mu = mu(valid_mask);

% Farrow 立方插值滤波器系数 (Cubic Lagrange)
h0 = -(mu - 1) .* (mu - 2) .* mu / 6;
h1 =  (mu - 1) .* (mu + 1) .* (mu - 2) / 2;
h2 = -(mu + 1) .* mu .* (mu - 2) / 2;
h3 =  (mu + 1) .* (mu - 1) .* mu / 6;

% 卷积求和（无边界失真效应）
x_sro = h0 .* x_cfo_filtered(idx_base - 1) + ...
        h1 .* x_cfo_filtered(idx_base) + ...
        h2 .* x_cfo_filtered(idx_base + 1) + ...
        h3 .* x_cfo_filtered(idx_base + 2);
x_sro = x_sro(:); % 转为列向量

%% 绘制循环前缀同步相关峰
% 求取接收信号各样点处的定时度量结果
N = 1024;
L = 32;

%构造干扰信号
%干扰循环前缀同步算法
disturb_list = [];
W = 4;  %窗函数带宽，可以调整干扰能量分布形状
delay_L = 0;
disturb_L = L + delay_L * 2;
for i = 1 : 2 : disturb_L
    disturb_list = cat(2, disturb_list, exp(1j * 2 * pi * i / disturb_L) * sinc((i - disturb_L / 2) / W));
    disturb_list = cat(2, disturb_list, exp(1j * 2 * pi * (i + 1) / disturb_L) * sinc((i - disturb_L / 2) / W));
end
disturb_power = 0;       %干扰信号平均功率
for i = 1 : disturb_L
    disturb_power = disturb_power + abs(disturb_list(i)) ^ 2; 
end
disturb_power = 2 * disturb_power / (disturb_L + N);
% 统计目标符号平均功率
signal_power = 0;
for i = N : (N + L) * 3
    signal_power = signal_power + abs(x_sro(i)) ^ 2;
end
signal_power = signal_power / (N + L) / 3;

%施加干扰信号
SJR = 0;
x = -8 : 2 : 28;    %横坐标，定时偏移量
maxPeakList = zeros(1, length(x));     %存储干扰后OFDM符号附近的最高相关峰峰值
for l = 1 : 2
    if l == 2
        SJR = 7;
    end
    currentPeakList = [];
    for offset = x
        disturbed_list = x_sro(1 : (10 * (N + L)));
        % offset = 16;      %干扰信号定时偏移
        %施加循环前缀压制性干扰
        for i = (N - delay_L + offset) : (N + L) : (N + L) * 4
            for j = 0: disturb_L - 1
                if mod(j, 2) == 0
                    disturbed_list(i + j) = x_sro(i + j) + disturb_list(j + 1) / sqrt(disturb_power) * sqrt(signal_power / (10 ^ (SJR / 10)));
                else
                    disturbed_list(i + j) = x_sro(i + j) - disturb_list(j + 1) / sqrt(disturb_power) * sqrt(signal_power / (10 ^ (SJR / 10)));
                end
                disturbed_list(i + N + j) = x_sro(i + N + j) + disturb_list(j + 1) / sqrt(disturb_power) * sqrt(signal_power / (10 ^ (SJR / 10)));
            end
        end

        % %施加加性高斯白噪声干扰
        % x_sro = awgn(x_sro(1:((N + L) * 8)), SJR, 'measured');

        result = [];
        mysum_list = [];
        power_list = [];
        for i = 1 : (N + L) * 7
            mysum = 0;
            power = 0;
            for j = 0: L - 1
                mysum = mysum + disturbed_list(i + j) * conj(disturbed_list(i + j + N));
                power = power + (abs(disturbed_list(i + j)) ^ 2) + (abs(disturbed_list(i + j + N)) ^ 2);
            end
            mysum_list = cat(2, mysum_list, abs(mysum));
            power_list = cat(2, power_list, power);
        %             result = cat(2, result, abs(mysum));
            result = cat(2, result, 2 * abs(mysum) / power);
        end
    %     figure(1);
    %     plot(result);
    %     title("添加信干比为0dB、定时滞后16个样点的循环前缀干扰时的循环前缀相关峰");
    %     figure(2);
    %     plot(mysum_list);
    %     title("添加信干比为0dB的循环前缀干扰时PSS符号的循环前缀自相关值");
    %     figure(3);
    %     plot(power_list);
    %     title("添加信干比为0dB的循环前缀干扰时PSS符号的相关窗口功率");

        % 记录当前干扰后的最高相关峰
        currentPeakList = cat(2, currentPeakList, max(result(1500: 4000)));
    end
    maxPeakList = cat(1, maxPeakList, currentPeakList);
end

figure(4);
hold on;
plot(x, maxPeakList(2, :), 'r-');
plot(x, maxPeakList(3, :), 'g-');
plot(x, zeros(1, length(x)) + 0.5, 'b--');
legend("0dB CP压制", "7dB CP压制", 'Location', 'southeast','Fontsize', 6);
xlim([-8, 28]);
xticks(x);
xticklabels(x);
ylim([0, 1]);
yticks(0: 0.1: 1);
yticklabels(0: 0.1: 1);
xlabel("循环前缀压制性干扰定时偏移");
ylabel("残存归一化最高峰值（%）");
title("循环前缀压制性干扰定时敏感度对比（0dB vs 7dB）");
hold off;


% %% 5. 定位 SSS 符号与时域分析
% fprintf('Analyzing Time Domain Signal...\n');
% 
% % 为了能往前看，我们把基础点往后挪20，也就是从重采样后数据的第21个点开始作为我们的 "0 延时基准"
% % 这是因为 Farrow 插值等操作如果不舍弃前面的数据会导致索引报错，我们强制给出一个安全 Buffer
% base_start_idx = 21; 
% 
% fig_const = figure('Position', [300, 300, 600, 600], 'Name', 'SSS Demodulation Result');
% 
% 
% % 遍历偏移量从 -20 到 +20
% for offset = -5 : 0
%     sss_start_idx_60 = base_start_idx + offset;
%         
%     % 防止越界
%     if sss_start_idx_60 < 1 || sss_start_idx_60 + N_fft > length(x_sro)
%         warning('Offset %d 导致下标越界，跳过。', offset);
%         continue;
%     end
%     
%     x_sss_time = x_sro(sss_start_idx_60 : sss_start_idx_60 + N_fft - 1);
%     
%     %% 6. OFDM 解调 (FFT) 与 小数定时偏差/相位纠正
%     x_sss_freq = fft(x_sss_time, N_fft) / sqrt(N_fft);  
%     
%     % 在这里我们先提取出有效子载波用来做相位拟合 (去除空载波对拟合的影响)
%     half_sc = 412;
%     idx_pos = 2 : (half_sc + 1);
%     idx_neg = (N_fft - half_sc + 1) : N_fft;
%     
%     % 将正、负频率轴拼接并映射出一条统一的 -half_sc 到 +half_sc 逻辑坐标轴
%     syms_valid = [x_sss_freq(idx_neg); x_sss_freq(idx_pos)]; 
%     freq_indices = [-half_sc:-1, 1:half_sc].';
% 
%     % 1. 消除调制相位的四次方操作
%     syms_pow4 = syms_valid.^4;
%     
%     % 2. 提取包裹后的相位，并展开 (unwrap)，因为相位随频率可能是线性变化且跨越平面的
%     phase_pow4 = unwrap(angle(syms_pow4));
%     
%     % 3. Polyfit 线性拟合: Phase = P(1)*f + P(2)
%     % 其中 P(1) 代表由于小数定时偏差引起的一阶斜率，P(2) 代表系统初相
%     P = polyfit(freq_indices, phase_pow4, 1);
%     
%     % 4. 计算补偿用相位：除以 4 恢复出原始偏移
%     % 这里我们构建全体1024点的补偿向量以作用于所有的子载波
%     full_freq_indices = [(0:N_fft/2-1), (-N_fft/2:-1)].';
%     phase_correction = (P(1) * full_freq_indices + P(2)) / 4;
%     
%     % 5. 应用补偿
%     x_sss_freq_corr = x_sss_freq .* exp(-1j * phase_correction);
% 
%     % 提取所有子载波用于最终显示
%     syms_all = x_sss_freq_corr(:);
%     
%     %% 7. 绘图
%     clf(fig_const); % 清除上一帧的图
%     plot(real(syms_all), imag(syms_all), 'b.', 'MarkerSize', 8);
%     hold on; grid on; axis square;
%     
%     th = 0:0.01:2*pi;
%     plot(cos(th), sin(th), 'k:', 'LineWidth', 0.5); % 单位圆
%     
%     title(sprintf('SSS Constellation\nStart Offset: %d 点\nResidual STO slope: %.4f\n(Press Enter for next)', offset, P(1)/4));
%     xlabel('I'); ylabel('Q');
%     xlim([-2 2]); ylim([-2 2]);
%     
%     fprintf('\n显示 Offset = %d 的结果。按 Enter 继续...\n', offset);
%     pause;
% end
% 
% fprintf('\n所有偏移点遍历结束。\n');

% 下面是原本的一些被注释掉的相位修正逻辑 (跳过执行以防报错)
return;

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

