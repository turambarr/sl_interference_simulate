% estimate_sss_cp_60mhz.m
% 先降采样到60MHz，再在60MHz基带寻找精准 CP 长度与起点
clear; clc; close all;

inFile = 'sigtest1.iq';
% 稍微往前多读一点，防止重采样边缘效应吃掉开头的 CP
startSample = 16386+874*6-200; 
readLen = 10000; 

fprintf('Loading and Resampling data...\n');
[x_raw, ~] = iq_read_int16_le(inFile, startSample, readLen);
x_raw = double(x_raw);
x_raw = x_raw - mean(x_raw);

% --- 新增：频谱下变频 (DDC 归基带) ---
freq_shift_hz = 63e6; % 频谱归基带向左搬移 63MHz
fprintf('Shifting spectrum left by %.1f MHz at 409.6MHz...\n', freq_shift_hz/1e6);
t_vec = (0:length(x_raw)-1) / 409.6e6;
x_raw = x_raw .* exp(-1j * 2 * pi * freq_shift_hz * t_vec).';

% --- 新增：使用零相移(Zero-Phase)低通滤波器抗混叠 ---
% 目标下降到 60MHz 采样率，奈奎斯特频率为 30MHz。
% 原始采样率为 409.6MHz，归一化截止频率 Wn = 30 / (409.6 / 2)
fprintf('Applying Zero-Phase Anti-aliasing LPF...\n');
Wn = 30e6 / (409.6e6 / 2);
% 设计一个适中阶数 (例如 30 阶) 的 FIR 滤波器
lpf_order = 30;
b_lpf = fir1(lpf_order, Wn);
% filtfilt 执行正反双向滤波：彻底消除群时延(相位偏移)，完美保护 CP 在时域上的绝对位置！
x_raw_filtered = filtfilt(b_lpf, 1, x_raw);

% 1. 重采样到 60MHz (采用 Farrow 分数阶延迟插值消除 resample 带来的边界失真)
fs_source = 409.6e6;
fs_target = 60e6;

T_in = 1 / fs_source;
T_out = 1 / fs_target;
% 避免两端由于缺少足量插值基准点而越界，舍弃最后微小的尾部
t_out = 0 : T_out : (length(x_raw_filtered)-3)*T_in; 

% 映射到输入序列的虚拟索引 (1-based)
idx_frac = t_out / T_in + 1; 
idx_base = floor(idx_frac);
mu = idx_frac - idx_base;

% 确保 base 索引不会越界 (Farrow Cubic 需要 base-1 到 base+2 共4个点)
valid_mask = (idx_base >= 2) & (idx_base <= length(x_raw_filtered)-2);
idx_base = idx_base(valid_mask);
mu = mu(valid_mask);

% Farrow 立方插值滤波器系数 (Cubic Lagrange)
h0 = -(mu - 1) .* (mu - 2) .* mu / 6;
h1 =  (mu - 1) .* (mu + 1) .* (mu - 2) / 2;
h2 = -(mu + 1) .* mu .* (mu - 2) / 2;
h3 =  (mu + 1) .* (mu - 1) .* mu / 6;

% 卷积求和（无边界失真效应），此处的源数据已替换为抗混叠滤波后的结果
x_60 = h0 .* x_raw_filtered(idx_base - 1) + ...
       h1 .* x_raw_filtered(idx_base) + ...
       h2 .* x_raw_filtered(idx_base + 1) + ...
       h3 .* x_raw_filtered(idx_base + 2);
x_60 = x_60(:); % 转为列向量

N = length(x_60);
% 60MHz下的有效数据体长度是确切的 1024 点
delay_D = 1024; 

% 2. 逐点共轭乘积
corr_prod = x_60(1 : N - delay_D) .* conj(x_60(delay_D + 1 : end));

% 3. 遍历确切的基带 CP 长度 (通常在 10~150 之间，例如常见的 72, 144 等)
W_test = 10 : 1 : 150;
peak_vals = zeros(size(W_test));
for i = 1:length(W_test)
    W = W_test(i);
    y = filter(ones(W, 1), 1, corr_prod);
    peak_vals(i) = max(abs(y));
end

% 4. Kneedle 算法找数学拐点
x = W_test(:); y = peak_vals(:);
xn = (x - min(x)) / (max(x) - min(x));
yn = (y - min(y)) / (max(y) - min(y));
dist = yn - xn;
[~, knee_idx] = max(dist);
exact_cp_60 = x(knee_idx);

% 5. 利用上述讨论的二次平滑机制寻找极其稳定的起点
% 整段窗求和
y_optimal = filter(ones(exact_cp_60, 1), 1, corr_prod);
% 局部平滑，去除山顶毛刺
y_smooth = smoothdata(abs(y_optimal), 'movmean', 5); 
[~, max_idx] = max(y_smooth);
cp_start_offset_60 = max_idx - exact_cp_60 + 1;

% 6. 核心：计算可以直接填进 sss_demodulation.m 的 startSample
% SSS 有效数据的起点(60MHz下) = CP的起点 + CP的长度
sss_data_start_60 = cp_start_offset_60 + exact_cp_60;
% 将偏移量按比例换算回 409.6MHz 的原生采样点位置
offset_4096 = round((sss_data_start_60 - 1) * (409.6 / 60));
global_startSample = startSample + offset_4096;

fprintf('\n=== 60MHz 基带 CP 分析完成 ===\n');
fprintf('60MHz 下精确 CP 长度 : %d 个点\n', exact_cp_60);
fprintf('CP 相对截取区(60MHz)起点偏移 : %d 个点\n', cp_start_offset_60);
fprintf('\n>>> 请将以下值填入 sss_demodulation.m 的 startSample:\n');
fprintf('startSample = %d;\n', global_startSample);
fprintf('==============================\n');

figure;
plot(x, y, 'b-', 'LineWidth', 2); hold on; grid on;
plot([x(1), x(end)], [y(1), y(end)], 'k--');
plot(exact_cp_60, y(knee_idx), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
xline(exact_cp_60, 'r--'); yline(y(knee_idx), 'r--');
title('60MHz基带下相关峰值上升曲线及拐点');