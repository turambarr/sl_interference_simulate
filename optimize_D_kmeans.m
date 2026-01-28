% optimize_D_kmeans.m
% 自动搜索最佳降采样倍率 D
% 算法：黄金分割搜索 (Golden Section Search) - 类似于二分法，用于查找单峰函数的极小值
% 评价指标：K-Means(k=4) 聚类后的簇内平均方差
% 目的：找到让星座图最聚拢的 D 值

clear; clc; close all;

%% 1. 基础参数
inFile = 'sigtest1.iq'; 
startSample = 15530-874; 
readLen     = 874 * 8; % 读取足够长的数据以获得稳定的统计结果

base_D = 6.398; % 当前 D 值 (3.199 * 2)
range_width = 0.01;
bounds = [base_D - range_width, base_D + range_width];

max_iter = 1000;
tol = 1e-5; % 停止阈值

%% 2. 读取数据
[x_raw, ~] = iq_read_int16_le(inFile, startSample, readLen);
x_raw = double(x_raw);
x_raw = x_raw - mean(x_raw); 
x_raw = x_raw / mean(abs(x_raw)); % 归一化
t_raw = 0:(length(x_raw)-1);

%% 3. 黄金分割搜索初始化
phi = (sqrt(5) - 1) / 2; % 黄金比例 0.618...
a = bounds(1);
b = bounds(2);
c = b - phi * (b - a);
d = a + phi * (b - a);

fprintf('开始搜索最佳 D 值...\n');
fprintf('初始范围: [%.4f, %.4f]\n', a, b);
fprintf('最大轮数: %d\n', max_iter);

% 计算初始点的评估值
fc = evaluate_variance(c, t_raw, x_raw);
fd = evaluate_variance(d, t_raw, x_raw);

history_D = zeros(max_iter, 1);
history_Var = zeros(max_iter, 1);

%% 4. 搜索循环
for k = 1:max_iter
    if abs(c - d) < tol
        break;
    end
    
    if fc < fd
        % 最小值在 [a, d] 之间
        b = d;
        d = c;
        fd = fc;
        c = b - phi * (b - a);
        fc = evaluate_variance(c, t_raw, x_raw);
        best_val = fc;
        best_D = c;
    else
        % 最小值在 [c, b] 之间
        a = c;
        c = d;
        fc = fd;
        d = a + phi * (b - a);
        fd = evaluate_variance(d, t_raw, x_raw);
        best_val = fd;
        best_D = d;
    end
    
    history_D(k) = best_D;
    history_Var(k) = best_val;
    
    fprintf('Iter %d: Range=[%.4f, %.4f], Best D=%.4f, Min Var=%.4f\n', k, a, b, best_D, best_val);
end

fprintf('\n搜索完成!\n');
fprintf('最佳 D 值: %.6f\n', best_D);
fprintf('最小方差: %.6f\n', best_val);

%% 5. 结果验证绘图
figure('Position', [100, 100, 1000, 500], 'Name', 'Optimization Result');

% 绘制收敛曲线
subplot(1, 2, 1);
valid_idx = find(history_D > 0);
yyaxis left;
plot(valid_idx, history_D(valid_idx), 'b.-'); ylabel('D Value');
yyaxis right;
plot(valid_idx, history_Var(valid_idx), 'r.-'); ylabel('Variance');
xlabel('Iteration'); title('Optimization Convergence');
grid on;

% 绘制最佳 D 下的星座图
subplot(1, 2, 2);
[~, best_x_resampled, best_offset] = evaluate_variance(best_D, t_raw, x_raw);
plot(best_x_resampled, '.', 'MarkerSize', 4);
axis square; grid on;
title(sprintf('Best D = %.6f (Offset=%.2f)', best_D, best_offset));
xlabel('I'); ylabel('Q');
xlim([-2 2]); ylim([-2 2]);


%% === 辅助函数：插值 + K-Means 评估 ===
function [min_var, best_samples, best_off] = evaluate_variance(D_val, t_raw, x_raw)
    % 评估给定 D_val 下的星座图质量
    % 遍历固定相位偏移 [0, 0.25, 0.5, 0.75] * D_val
    % 选取其中最小的方差作为该 D_val 的评分
    
    interpolation_method = 'spline';
    
    offsets_ratios = [0, 0.25, 0.5, 0.75];
    min_var = inf;
    best_samples = [];
    best_off = 0;
    
    for r = offsets_ratios
        offset = r * D_val;
        t_new = offset : D_val : (length(x_raw)-1);
        
        if length(t_new) < 50
             this_var = 1e9;
             x_resampled = [];
        else
            x_resampled = interp1(t_raw, x_raw, t_new, interpolation_method);
            
            % 将复数转为 [Real, Imag] 矩阵供 kmeans 使用
            X_cluster = [real(x_resampled(:)), imag(x_resampled(:))];
            
            try
                % k=4聚类
                [~, ~, sumd] = kmeans(X_cluster, 4, 'MaxIter', 10, 'Replicates', 1, 'Display', 'off');
                % 平均方差
                this_var = sum(sumd) / length(x_resampled);
            catch
                this_var = 1e9; 
            end
        end
        
        if this_var < min_var
            min_var = this_var;
            best_samples = x_resampled;
            best_off = offset;
        end
    end
end
