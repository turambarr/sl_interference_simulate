% check_modulation_type.m
% 功能：基于插值后的符号，分析差分域特性 (Y_k * Y_{k-1}^*)
% 目的：区分 QPSK (4点) 和 对称 DPSK (差分域2点: +/- j)

clear; clc; close all;

%% === 1. 参数设置 (保持与 constellation_interp.m 一致) ===
inFile = 'sigtest8.iq'; 
startSample = 19934-5*874;   % 最佳信号起始点
readLen     = 874*8;     % 读取长度

D = 6.3975;            % 最佳降采样倍率
interpolation_method = 'spline';

% --- 频率校正参数 ---
fs = 409.6e6;          % 采样率
f_offset = 1000;       % 假设频偏 1kHz

%% === 2. 读取与预处理 ===
fprintf('读取文件: %s (Start=%d, Len=%d)\n', inFile, startSample, readLen);
[x_raw, ~] = iq_read_int16_le(inFile, startSample, readLen);
x_raw = double(x_raw);
x_raw = x_raw - mean(x_raw); 

% 归一化幅度
x_raw = x_raw / mean(abs(x_raw));

% --- 去除频偏 ---
if f_offset ~= 0
    fprintf('执行去频偏 (Offset=%.1f Hz)...\n', f_offset);
    if isrow(x_raw), x_raw = x_raw(:); end
    t_vec = (0 : length(x_raw)-1)' / fs;
    x_raw = x_raw .* exp(-1i * 2 * pi * f_offset * t_vec);
    x_raw = x_raw.'; 
end

%% === 3. 插值降采样提取符号 ===
t_raw = 0:(length(x_raw)-1);
offset_ratio = 0; 
current_offset = offset_ratio * D;

% 构造采样时间点
t_new = current_offset : D : (length(x_raw)-1);

% 执行插值
x_syms = interp1(t_raw, x_raw, t_new, interpolation_method);

% 再次幅度归一化 (针对符号层面)
x_syms = x_syms / mean(abs(x_syms));

%% === 4. 计算差分域符号 ===
% Z[k] = Y[k] * conj(Y[k-1])
y_curr = x_syms(2:end);
y_prev = x_syms(1:end-1);

z_diff = y_curr .* conj(y_prev);

% 对差分符号进行归一化以便观察
z_diff = z_diff / mean(abs(z_diff));

%% === 5. 绘图对比 ===
figure('Position', [100, 100, 1000, 500], 'Name', 'Modulation Analysis');

% 子图1: 原始接收星座图 (Static)
subplot(1, 2, 1);
plot(x_syms, '.', 'Color', 'b', 'MarkerSize', 8);
hold on;
xline(0, 'k:'); yline(0, 'k:');
title(sprintf('原始符号 Y_k (Static)\n(呈现 4 个簇说明是 QPSK 或 对称DPSK)'));
xlabel('I'); ylabel('Q');
axis square; grid on;
xlim([-2 2]); ylim([-2 2]);

% 子图2: 差分域星座图 (Dynamic)
subplot(1, 2, 2);
plot(z_diff, '.', 'Color', 'r', 'MarkerSize', 8);
hold on;
xline(0, 'k:'); yline(0, 'k:');
title(sprintf('差分符号 Z_k = Y_k \\cdot Y_{k-1}^*\n(Dynamic / Differential Domain)'));
xlabel('Diff I'); ylabel('Diff Q');
axis square; grid on;
xlim([-2 2]); ylim([-2 2]);

% --- 辅助判断说明 ---
annotation('textbox', [0.15, 0.02, 0.7, 0.1], 'String', ...
    {'判决辅助:', ...
     '1. 如果右图呈现 4 个点 ( 上下左右 ) -> 可能是 QPSK (相邻符号相位差可以是 0, 90, 180, 270)', ...
     '2. 如果右图呈现 2 个点 ( 仅在虚轴 +/- j ) -> 极可能是 对称DPSK (pi/2-BPSK, 相位跳变仅 +/- 90度)', ...
     '3. 如果右图呈现 2 个点 ( 仅在实轴 +/- 1 ) -> 普通 BPSK/DPSK'}, ...
    'FitBoxToText','on', 'EdgeColor','none', 'BackgroundColor','w');

fprintf('\n=== 分析完成 ===\n');
fprintf('符号总数: %d\n', length(x_syms));
fprintf('请观察右侧差分星座图的分布情况。\n');
