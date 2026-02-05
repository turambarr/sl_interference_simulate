% cross_corr_interp.m
% 功能：先对两段信号进行基于插值的降采样 (interp1)，然后计算互相关
% 目的：分析在符号速率层面的相关性，去除过采样带来的冗余计算，以及验证重采样后的波形匹配度

clear; clc; close all;

%% === 1. 参数设置 ===

% --- 信号 1 (模板) ---
file1 = 'sigtest1.iq';
start1 = 15530-874+10; % 参考 constellation_interp 中的位置
% 注意：为了保证降采样后有足够的点，这里读取原始样点数要足够多
% 如果原来的模板长度是 874, 降采样 D=6.4, 则大约剩下 136 个符号
len1_raw   = 874;     

% --- 信号 2 (搜索区域) ---
file2 = 'sigtest26.iq';
start2 = 0;
len2_raw   = 30000;    % 搜索范围

% --- 降采样参数 ---
D = 6.398;            % 降采样倍率 (请使用 optimize_D_kmeans 算出的最佳值)
interp_method = 'spline';

%% === 2. 读取原始数据 ===
fprintf('正在读取信号 1: %s (Start=%d, Len=%d)...\n', file1, start1, len1_raw);
[s1_raw, ~] = iq_read_int16_le(file1, start1, len1_raw);
s1_raw = double(s1_raw);
s1_raw = s1_raw - mean(s1_raw);

fprintf('正在读取信号 2: %s (Start=%d, Len=%d)...\n', file2, start2, len2_raw);
[s2_raw, ~] = iq_read_int16_le(file2, start2, len2_raw);
s2_raw = double(s2_raw);
s2_raw = s2_raw - mean(s2_raw);

%% === 3. 插值降采样 (Interp1) ===
fprintf('执行插值降采样 (D=%.4f)...\n', D);

% 对信号 1 降采样
t_raw1 = 0 : (length(s1_raw)-1);
t_new1 = 0 : D : (length(s1_raw)-1);
s1_down = interp1(t_raw1, s1_raw, t_new1, interp_method);

% 对信号 2 降采样
t_raw2 = 0 : (length(s2_raw)-1);
t_new2 = 0 : D : (length(s2_raw)-1);
s2_down = interp1(t_raw2, s2_raw, t_new2, interp_method);

% 处理可能的 NaN (如果 D 步进导致最后一个点越界)
s1_down(isnan(s1_down)) = [];
s2_down(isnan(s2_down)) = [];

fprintf('  信号 1: 原始 %d -> 降采样后 %d\n', length(s1_raw), length(s1_down));
fprintf('  信号 2: 原始 %d -> 降采样后 %d\n', length(s2_raw), length(s2_down));

% 归一化幅度 (可选，为了相关结果更有物理意义)
% 这里不对波形做功率归一化，只在互相关时做能量归一化
% 但为了数值稳定，可以先简单缩放一下
scale_factor = 1 / mean(abs(s1_down));
s1_down = s1_down * scale_factor;
s2_down = s2_down * scale_factor;

%% === 4. 计算互相关 ===
% 注意：这里的 Lag 单位变成了 "符号" (Symbols) 而不是原始采样点
[xc, lags] = xcorr(s2_down, s1_down, 'none');

% 能量归一化 (除以模板 s1_down 的能量)
template_energy = sum(abs(s1_down).^2);
if template_energy > 0
    xc = xc / template_energy;
end

xc_abs = abs(xc);
[max_val, max_idx] = max(xc_abs);
best_lag_symbol = lags(max_idx);

%% === 5. 绘图展示 ===
figure('Position', [100, 100, 1000, 600], 'Name', 'Cross-Correlation (Downsampled)');

% 子图1: 互相关模值
subplot(2, 1, 1);
plot(lags, xc_abs, 'b', 'LineWidth', 1);
title(sprintf('降采样互相关模值 (D=%.4f)\n最大匹配度: %.4f @ Lag=%d (Symbols)', D, max_val, best_lag_symbol));
ylabel('归一化模值');
xlabel('Lag (符号数)');
grid on; axis tight;
xline(best_lag_symbol, 'r--');

% 标记最大值
hold on;
plot(best_lag_symbol, max_val, 'ro', 'MarkerSize', 6, 'LineWidth', 2);

% 子图2: 转换回原始样点域的估算
subplot(2, 1, 2);
% 将 Lag(符号) 换算回 原始样点 Lag
lags_raw_approx = lags * D;
plot(lags_raw_approx, xc_abs, 'k', 'LineWidth', 1);
title(sprintf('映射回原始样点域 (Lag * D)\n估计原始偏移约: %.1f 样点', best_lag_symbol * D));
ylabel('模值');
xlabel('Lag (近似原始样点数)');
grid on; axis tight;
xline(best_lag_symbol * D, 'r--');

fprintf('\n=== 分析结果 ===\n');
fprintf('最大互相关值: %.4f\n', max_val);
fprintf('最佳对齐位置 (降采样域): Lag = %d (符号)\n', best_lag_symbol);
fprintf('推算原始对齐位置: Lag * D ≈ %.2f (原始样点)\n', best_lag_symbol * D);
