% find_best_match_across_files.m
% 遍历 sigtest2.iq 到 sigtest144.iq
% 与 sigtest1.iq 的指定模板段进行互相关
% 找出信噪比（互相关峰值）最高的那个文件和位置

clear; clc;

%% === 1. 模板设置 ===
template_file = 'sigtest1.iq';
temp_start = 15514;
temp_len   = 874;

fprintf('读取模板: %s (Start=%d, Len=%d)\n', template_file, temp_start, temp_len);
[s1, ~] = iq_read_int16_le(template_file, temp_start, temp_len);
s1 = double(s1);
s1 = s1 - mean(s1);
template_energy = sum(abs(s1).^2);

if template_energy == 0
    error('模板能量为0，请检查读取范围');
end

%% === 2. 遍历搜索 ===
file_indices = 2:144; % 目标文件范围
results = struct('file', {}, 'max_val', {}, 'lag', {});
cnt = 0;

fprintf('开始扫描 %d 个文件...\n', length(file_indices));

for k = file_indices
    target_file = sprintf('sigtest%d.iq', k);
    
    if ~isfile(target_file)
        fprintf('  [Skipped] %s 不存在\n', target_file);
        continue;
    end
    
    % 读取目标文件 (读取整个文件，设一个足够大的值，比如 50000)
    % 也可以先dir看大小，这里直接读50000
    [s2, meta] = iq_read_int16_le(target_file, 0, 100000); 
    if meta.numSamplesRead < temp_len
        continue; 
    end
    
    s2 = double(s2);
    s2 = s2 - mean(s2);
    
    % 计算互相关
    % 使用 'none' 模式然后手动归一化
    [xc, lags] = xcorr(s2, s1, 'none');
    xc_max_val = max(abs(xc)) / template_energy;
    
    cnt = cnt + 1;
    results(cnt).file = target_file;
    results(cnt).max_val = xc_max_val;
    
    % 找到最大值对应的 lag
    [~, max_idx] = max(abs(xc));
    results(cnt).lag = lags(max_idx);
    
    % 简单的进度条
    if mod(cnt, 20) == 0
        fprintf('  已按照处理 %d 个文件...\n', cnt);
    end
end

%% === 3. 结果汇总 ===
% table 方便排序
T = struct2table(results);
T_sorted = sortrows(T, 'max_val', 'descend');

fprintf('\n=== 扫描完成 ===\n');
fprintf('Top 5 匹配度最高的文件:\n');
disp(T_sorted(1:min(5, height(T_sorted)), :));

best_file = T_sorted.file{1};
best_val  = T_sorted.max_val(1);
best_lag  = T_sorted.lag(1);

fprintf('\n【最佳结果】\n文件: %s\n峰值: %.4f\n位置: Lag=%d\n', best_file, best_val, best_lag);

%% === 4. 画出最佳那个的图 ===
fprintf('\n正在绘制最佳结果...\n');
figure('Position', [100, 100, 1000, 600], 'Name', 'Best Match Analysis');

% 重读最佳文件
[s2_best, ~] = iq_read_int16_le(best_file, 0, 100000);
s2_best = double(s2_best) - mean(double(s2_best));

[xc_best, lags_best] = xcorr(s2_best, s1, 'none');
xc_abs_best = abs(xc_best) / template_energy;

plot(lags_best, xc_abs_best, 'b');
grid on; axis tight;
title(sprintf('最佳匹配: %s (Max Val = %.4f)', best_file, best_val));
xlabel('Lag'); ylabel('归一化互相关');
hold on;
plot(best_lag, best_val, 'ro', 'MarkerSize', 8, 'LineWidth', 2);
text(best_lag, best_val, sprintf(' File=%s\n Val=%.2f', best_file, best_val), ...
    'VerticalAlignment', 'bottom');
