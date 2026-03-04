% estimate_sss_cp_60mhz.m
% 先降采样到60MHz，再在60MHz基带寻找精准 CP 长度与起点
clear; clc; close all;

inFile = 'sigtest8.iq';
% 稍微往前多读一点，防止重采样边缘效应吃掉开头的 CP
startSample = 15564 + 874*8 - 1000; 
readLen = 10000; 

fprintf('Loading and Resampling data...\n');
[x_raw, ~] = iq_read_int16_le(inFile, startSample, readLen);
x_raw = double(x_raw);
x_raw = x_raw - mean(x_raw);

% 1. 重采样到 60MHz
fs_source = 409.6e6;
fs_target = 60e6;
[P, Q] = rat(fs_target / fs_source);
x_60 = resample(x_raw, P, Q);

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