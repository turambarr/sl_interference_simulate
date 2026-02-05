% constellation_interp.m
% 使用插值法 (interp1) 进行非整数倍降采样并绘制星座图
% 目的：替代 resample 函数，更精确地控制采样时刻（相位），解决采样点对不准的问题

clear; clc; close all;

%% 1. 参数设置
inFile = 'sigtest8.iq'; 
startSample = 19934;   % 读取起始点
readLen     = 874;           % 读取长度 (适当读长一点，保证足够插值)

D = 6.3975;                 % 降采样倍率 (每个符号占用的原始样点数)
interpolation_method = 'spline'; % 插值方法: 'linear', 'spline', 'pchip'

%% 2. 读取原始数据
fprintf('读取文件: %s (Start=%d, Len=%d)\n', inFile, startSample, readLen);
[x_raw, ~] = iq_read_int16_le(inFile, startSample, readLen);
x_raw = double(x_raw);
x_raw = x_raw - mean(x_raw); % 去直流

% 归一化幅度
x_raw = x_raw / mean(abs(x_raw));

%% 3. 插值核心逻辑
% 原始时间轴: 0, 1, 2, ..., N-1
t_raw = 0:(length(x_raw)-1);

% 我们不再搜索最佳的采样相位，而是使用固定的 Offset
% 如果需要微调，请修改这里的 offset_ratio (0 ~ 1 之间)
offset_ratio = 0; 
current_offset = offset_ratio * D;

fprintf('使用固定 Offset: %.4f (Ratio=%.2f)\n', current_offset, offset_ratio);

%% 4. 执行单次插值与绘图
figure('Position', [100, 100, 800, 800], 'Name', 'Constellation (Fixed Offset)');

% --- 构造新的采样时间点 ---
% 从 current_offset 开始，每次步进 D，直到不超过原始数据长度
t_new = current_offset : D : (length(x_raw)-1);

% --- 执行插值 (interp1) ---
x_resampled = interp1(t_raw, x_raw, t_new, interpolation_method);

% --- 绘图 ---
plot(x_resampled, '.', 'Color', 'b', 'MarkerSize', 6);
axis square; grid on;
title(sprintf('插值法重采样星座图\nD = %.4f, Offset = %.2f * D', D, offset_ratio));
xlabel('In-Phase (I)'); ylabel('Quadrature (Q)');
xlim([-2 2]); ylim([-2 2]);

fprintf('\n=== 结果说明 ===\n');
fprintf('点数: %d\n', length(x_resampled));
