% cross_corr_sig1_sig2.m
% 选取 signal1 的一段作为参考（Template）
% 选取 signal2 的一段作为搜索区域（Target）
% 计算滑动的互相关系数，寻找相似结构

clear; clc;

%% === 参数设置区域 ===

% --- 信号 1 (参考信号/模板) ---
file1 = 'sigtest1.iq';     % 文件名
start1 = 15514;                % 起始点（0-based）
len1   = 874;             % 长度 (例如选一个特定Burst或PSS长度)

% --- 信号 2 (待检测信号/长序列) ---
file2 = 'sigtest144.iq';     % 文件名
start2 = 0;                % 起始点（0-based）
len2   = 30000;             % 长度

% --- 分析参数 ---
use_normalized_xcorr = true; % 是否使用 'coeff' 归一化 (需要信号处理工具箱)
                             % 如果设为 true, 输出范围理论在 [-1, 1]
                             % 如果设为 false, 输出为原始能量积

%% === 读取数据 ===
fprintf('正在读取信号 1: %s (Start=%d, Len=%d)...\n', file1, start1, len1);
[s1, ~] = iq_read_int16_le(file1, start1, len1);
s1 = double(s1);

fprintf('正在读取信号 2: %s (Start=%d, Len=%d)...\n', file2, start2, len2);
[s2, ~] = iq_read_int16_le(file2, start2, len2);
s2 = double(s2);

% 简单预处理 (去直流)
s1 = s1 - mean(s1);
s2 = s2 - mean(s2);

%% === 计算互相关 (xcorr) ===
% calc: R(lag) = sum(s2(n) * conj(s1(n-lag)))
% 也就是 s1 在 s2 上滑动
fprintf('计算互相关...\n');

if use_normalized_xcorr
    % xcorr 的 'coeff' 选项要求两信号等长，且它是全局归一化（除以 sqrt(E1*E2)），会导致长信号匹配时峰值很小。
    % 这里我们改为：除以“模板信号的能量”。
    % 物理意义：如果 s2 中某一段完美匹配 s1（且幅度一致），则峰值为 1.0。
    [xc, lags] = xcorr(s2, s1, 'none'); 
    
    template_energy = sum(abs(s1).^2);
    if template_energy > 0
        xc = xc / template_energy;
    end
    
    ylabel_str = '归一化互相关 (相对于模板能量)';
    ylim_range = [-1.2, 1.2]; % 稍微留点余量
else
    [xc, lags] = xcorr(s2, s1);          % 原始互相关
    ylabel_str = '原始相关幅值';
    ylim_range = 'auto';
end

% 分离实部和幅度
% xc_real = real(xc); %不再需要实部
xc_abs  = abs(xc);

% 找出最大峰值即最佳对齐Lag
[max_val, max_idx] = max(xc_abs);
best_lag = lags(max_idx);

%%  绘图展示 
figure('Position', [100, 100, 1000, 500], 'Name', 'Cross-Correlation Analysis');

% 互相关模值 (Magnitude) - 代表“相似度强度”
plot(lags, xc_abs, 'b', 'LineWidth', 1); 
title(sprintf('互相关模值 |R| (Magnitude) - 仅表示匹配强度，忽略相位\n最大模值: %.4f @ Lag=%d', max_val, best_lag));
ylabel('模值 (Abs)');
xlabel('Lag (延迟采样数)');
grid on; axis tight;
xline(best_lag, 'r--'); % 标记最佳位置

% 标记最大值点
hold on;
plot(best_lag, max_val, 'ro', 'MarkerSize', 8, 'LineWidth', 2);
text(best_lag, max_val, sprintf(' Lag=%d\n Val=%.2f', best_lag, max_val), ...
    'VerticalAlignment', 'bottom');

fprintf('\n=== 分析结果 ===\n');
fprintf('最大相关模值: %.4f\n', max_val);
fprintf('对应 Lag: %d\n', best_lag);
fprintf('这意味着: 信号1 向右延迟 %d 个样点后与 信号2 最匹配。\n', best_lag);



