% estimate_sss_cp.m
% 精确寻找 CP 长度 - 步长为1的高精度遍历 + Kneedle
clear; clc; close all;
inFile = "sigtest8.iq";
startSample = 15564 + 874*8-200;
readLen = 20000;
delay_D = 6992;
fprintf("Loading data...\\n");
[x_raw, ~] = iq_read_int16_le(inFile, startSample, readLen);
x_raw = double(x_raw);
x_raw = x_raw - mean(x_raw);
N = length(x_raw);
corr_prod = x_raw(1 : N - delay_D) .* conj(x_raw(delay_D + 1 : end));
W_test = 200 : 1 : 800;
peak_vals = zeros(size(W_test));
for i = 1:length(W_test), W = W_test(i); y = filter(ones(W, 1), 1, corr_prod); peak_vals(i) = max(abs(y)); end
x = W_test(:); y = peak_vals(:);
xn = (x - min(x)) / (max(x) - min(x));
yn = (y - min(y)) / (max(y) - min(y));
dist = yn - xn;
[~, knee_idx] = max(dist);
exact_cp_4096 = x(knee_idx);
exact_cp_60 = exact_cp_4096 * (60 / 409.6);

% --- 新增：计算精确的 CP 起始位置偏移量 ---
% 使用找到的高精度最优 CP 长度，计算一次该长度下的滑动相关序列
y_optimal = filter(ones(exact_cp_4096, 1), 1, corr_prod);
[~, max_idx] = max(abs(y_optimal));
% filter 给出的是滑动窗结束时的值，因此 CP 在当前读取的 x_raw 数据块中的起始索引是:
cp_start_offset = max_idx - exact_cp_4096 + 1;
% 换算为全局绝对采样点位置 (Global Sample Index)
global_start_idx = startSample + cp_start_offset - 1;

fprintf("\n=== 高精度 CP 长度分析完成 ===\n");
fprintf("409.6MHz 原生采样率 CP 长度 : %d 个点\n", exact_cp_4096);
fprintf("60MHz 基带等效 CP 长度      : %.2f 个点\n", exact_cp_60);
fprintf("\n=== 精确 CP 起点分析结果 ===\n");
fprintf("CP 相对当前截取区起点偏移量 : %d 个采样点\n", cp_start_offset);
fprintf("CP 在文件中的绝对全局起止点 : 第 %d 到第 %d 采样点\n", global_start_idx, global_start_idx + exact_cp_4096 - 1);
fprintf("==============================\\n");
figure;
plot(x, y, "b-", "LineWidth", 2); hold on; grid on;
plot([x(1), x(end)], [y(1), y(end)], "k--");
plot(exact_cp_4096, y(knee_idx), "ro", "MarkerSize", 10, "MarkerFaceColor", "r");
xline(exact_cp_4096, "r--"); yline(y(knee_idx), "r--");
title("高精度相关峰值上升曲线及精确数学拐点 (Kneedle)");
