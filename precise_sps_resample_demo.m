% precise_sps_resample_demo.m
% 目的：
% 1. 利用 PSS 重复结构，通过线性回归极其精确地推算 SPS (Samples Per Symbol)。
% 2. 使用精确的 SPS 进行分数倍插值，消除采样率误差(SRO)导致的长期漂移。
% 3. 绘制星座图验证。

clear; clc; close all;

%% 1. 参数设置
filename = 'sigtest8.iq'; % 请确保此文件在您的路径中
start_offset = 19934 - 874 * 5;   % 用户指定的PSS起始位置
read_len = 874 * 12;      % 读取长度 (略多于8个块以确保足够)

% 粗略参数 (409.6MHz / (60MHz/128) )
% 理论值: 409.6 / 0.46875 = 873.8133...
block_len_approx = 874;   
num_blocks_to_track = 8;  % 多少个重复块用来估算斜率

% Constellation Interp / Costas 参数
freq_offset_init = 1000; % 如果已知大概频偏可填，否则 Costas 会自适应

%% 2. 读取数据
if ~isfile(filename)
    error('文件 %s 不存在，请检查路径。', filename);
end
fprintf('正在读取文件: %s ...\n', filename);
[x_raw, ~] = iq_read_int16_le(filename, start_offset, read_len);
x_raw = double(x_raw);
x_raw = x_raw - mean(x_raw);          % 去直流
x_raw = x_raw / mean(abs(x_raw));     % 归一化

%% 3. 第一步：确定突发起始点 (burst_start)
% 用户已指定读取位置为 PSS 起始点，因此 buffer 中的索引为 1
burst_start = 1;

fprintf('使用指定起始索引 (相对于读取缓冲): %d\n', burst_start);

%% 4. 第二步：精密 SPS 估算 (Linear Regression)
% 1. 提取第一个块作为模板
% 确保不越界
if burst_start + block_len_approx > length(x_raw)
    error('数据长度不足以提取模板');
end
template = x_raw(burst_start : burst_start + block_len_approx - 1);

% 2. 追踪后续块的位置
search_range = 20; % 搜索范围 (样本)
peak_indices_abs = zeros(1, num_blocks_to_track);
peak_indices_abs(1) = burst_start; 

fprintf('正在执行精密同步 (追踪 %d 个块)...\n', num_blocks_to_track);

for k = 2:num_blocks_to_track
    % 预测大概位置
    est_pos = burst_start + (k-1) * block_len_approx;
    
    % 截取搜索窗口
    s_start = max(1, est_pos - search_range);
    s_end   = min(length(x_raw), est_pos + block_len_approx + search_range * 2); % 稍微多读点
    
    if s_end - s_start < length(template)
        break; 
    end
    
    segment = x_raw(s_start : s_end);
    
    % 互相关找精确偏移
    % 方式: conv(segment, flip(conj(template)))
    cc = abs(conv(segment, flipud(conj(template))));
    [~, local_peak] = max(cc);
    
    % conv 的峰值位置换算：
    % length(segment) + length(template) - 1
    % 峰值索引对应 segment 中与 template 对齐的结束点
    % 换算回 segment 的起始对齐点：
    match_start_in_segment = local_peak - length(template) + 1;
    
    % 绝对位置
    abs_pos = s_start + match_start_in_segment - 1;
    peak_indices_abs(k) = abs_pos;
end

% 3. 线性回归算出斜率
valid_mask = peak_indices_abs > 0;
block_indices = (find(valid_mask) - 1)'; % 0, 1, 2...
detected_positions = peak_indices_abs(valid_mask)';

if length(block_indices) < 3
    error('无法追踪到足够的块来进行回归分析。请检查信号质量或 block_len_approx 参数。');
end

% 拟合: Pos = slope * block_index + intercept
p = polyfit(block_indices, detected_positions, 1);
precise_block_len_samples = p(1);

% 4. 计算 SPS
symbols_per_block = 128; % 已知 PSS 结构参数
SPS_precise = precise_block_len_samples / symbols_per_block;

fprintf('------------------------------------------------\n');
fprintf('参考 Block 长度: %d\n', block_len_approx);
fprintf('实测 Block 长度: %.6f (Linear Fit Slope)\n', precise_block_len_samples);
fprintf('粗略 SPS: %.6f\n', block_len_approx / symbols_per_block);
fprintf('精密 SPS: %.6f\n', SPS_precise);
fprintf('采样率偏差 (估算): %.2f ppm\n', (precise_block_len_samples - 873.81333)/873.81333 * 1e6);
fprintf('------------------------------------------------\n');

%% 5. 第三步：分数倍插值重采样
% 从 burst_start 开始重采样后续数据
data_to_resample = x_raw(burst_start : end);

% 生成非整数的采样时刻
num_symbols_available = floor(length(data_to_resample) / SPS_precise);
% 时刻点: 0.5*SPS, 1.5*SPS, 2.5*SPS ... (采样在符号中心)
sample_instants = ((0 : num_symbols_available-1) + 0.5) * SPS_precise + 1;

% 确保不越界
valid_idx = sample_instants <= length(data_to_resample);
sample_instants = sample_instants(valid_idx);

% 插值 (Spline 精度较高)
fprintf('正在利用精密 SPS 进行插值重采样 (%d symbols)...\n', length(sample_instants));
x_syms = interp1(1:length(data_to_resample), data_to_resample, sample_instants, 'spline').';

%% 6. 第四步：星座图验证 (带 Costas 环去频偏)
% 这里的代码直接复用之前验证过的 Costas 逻辑

% Costas 环参数
BL_T = 0.05; 
zeta = 0.707;
w_n = BL_T / (zeta + 1/(4*zeta));
alpha = 2 * zeta * w_n;
beta  = (w_n)^2;

phase_est = -pi/4;  % 初始相位，针对对角线锁定优化
freq_est  = 0;

N_sym = length(x_syms);
x_final = zeros(size(x_syms));
freq_log = zeros(1, N_sym);

for n = 1:N_sym
    val = x_syms(n);
    z = val * exp(-1j * phase_est);
    x_final(n) = z;
    
    % PED
    I = real(z); Q = imag(z);
    err = Q * sign(I) - I * sign(Q);
    
    % Loop Filter
    freq_est = freq_est + beta * err;
    phase_est = phase_est + freq_est + alpha * err;
    
    freq_log(n) = freq_est;
end

% 旋转回坐标轴
x_final = x_final * exp(-1j * pi/4);

%% 7. 绘图
figure('Position', [100, 100, 1000, 500], 'Name', 'Precise SPS Resampling Result');

subplot(1, 2, 1);
plot(block_indices, detected_positions, 'o', 'LineWidth', 1.5); hold on;
plot(block_indices, polyval(p, block_indices), 'r-', 'LineWidth', 1);
title('线性回归测量 SPS');
xlabel('Block Index'); ylabel('Sample Index');
legend('Measured Peaks', 'Linear Fit');
grid on;

subplot(1, 2, 2);
plot(real(x_final), imag(x_final), 'k.', 'MarkerSize', 8);
axis square; grid on; xline(0); yline(0);
title({sprintf('重采样后星座图 (SPS=%.5f)', SPS_precise); 'After Costas Correction'});
xlim([-2 2]); ylim([-2 2]);
xlabel('I'); ylabel('Q');
