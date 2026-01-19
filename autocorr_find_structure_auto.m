% 独立脚本：自动寻找满足特定位置和幅值要求的 W, D
% 目标：在样点索引 36300 之后，寻找自相关值 > 0.8 或 < -0.8 的结构
% 遍历范围：W = 64 ~ 1024 (步长1)，默认假设 D = W

clear; clc; close all;

%% 1. 参数设置
inFile = 'test2.iq';

startSample0 = 0;   
% 为了确保覆盖 36300 之后足够长的区域，建议读取长度适中
% 假设感兴趣的事件发生在 36300 之后不久
numSamples = 100000; 

normalizeToUnit = true; 
removeMean = true;      

scanRange = 64:2000; % 遍历 W 的范围
targetIndexThreshold = 30000; % 索引门限
corrThreshold = 0.8; % 相关值门限 (绝对值)

%% 2. 读取数据
[x, meta] = iq_read_int16_le(inFile, startSample0, local_num_samples(inFile, numSamples, startSample0));
L = meta.numSamplesRead;
if L == 0, error('未读取到数据。'); end
x = double(x);
if normalizeToUnit, x = x / 32768; end
if removeMean, x = x - mean(x); end

t_axis = (0:length(x)-1) + startSample0;

fprintf('数据已加载 (L=%d)。\n正在扫描 W=%d~%d, 寻找 Index > %d 且 |Corr| > %.1f 的结构...\n', ...
    length(x), min(scanRange), max(scanRange), targetIndexThreshold, corrThreshold);

%% 3. 遍历扫描
foundCandidates = []; % 存储结果: [W, D, MaxVal, MaxIdx, MinVal, MinIdx]

for W = scanRange
    D = W; % 默认 Schmidl & Cox 配置
    
    if length(x) <= D + W
        continue;
    end
    
    % --- 核心计算 (复用之前的高效逻辑) ---
    rx_delayed = x(1+D:end);
    rx_base    = x(1:end-D);
    
    conj_prod = conj(rx_base) .* rx_delayed;
    
    % 滑动求和
    P_metric = filter(ones(1, W), 1, conj_prod);
    
    % 能量归一化
    rx_power = abs(rx_base).^2;
    R_energy = filter(ones(1, W), 1, rx_power);
    
    % 计算实部度量
    M_complex = P_metric ./ (R_energy + 1e-10);
    M_real = real(M_complex);
    
    % 补齐长度以便与 t_axis 对齐
    M_full = [M_real; zeros(D, 1)];
    
    % --- 筛选条件 ---
    % 1. 找到所有 > 36300 的索引
    valid_mask = t_axis > targetIndexThreshold;
    
    if ~any(valid_mask)
        continue;
    end
    
    M_ROI = M_full(valid_mask);
    t_ROI = t_axis(valid_mask);
    
    % 2. 检查阈值
    [maxVal, maxLoc] = max(M_ROI);
    [minVal, minLoc] = min(M_ROI);
    
    isHit = false;
    if maxVal > corrThreshold
        candidate = struct('W', W, 'D', D, 'Type', 'Positive', 'Val', maxVal, 'Idx', t_ROI(maxLoc));
        foundCandidates = [foundCandidates; candidate];
        isHit = true;
    end
    
    if minVal < -corrThreshold
        candidate = struct('W', W, 'D', D, 'Type', 'Negative', 'Val', minVal, 'Idx', t_ROI(minLoc));
        foundCandidates = [foundCandidates; candidate];
        isHit = true;
    end
    
    if isHit
        % 可以在这里打印，或者存下来最后打印
        % fprintf('Found: W=%d, +Val=%.2f (@%d), -Val=%.2f (@%d)\n', W, maxVal, t_ROI(maxLoc), minVal, t_ROI(minLoc));
    end
end

%% 4. 结果汇总与展示
if isempty(foundCandidates)
    fprintf('未找到满足条件的 W。\n');
else
    fprintf('\n====== 搜索结果 (Top 20 by Amplitude) ======\n');
    % 按绝对值幅值排序
    [~, sortIdx] = sort(abs([foundCandidates.Val]), 'descend');
    sortedCandidates = foundCandidates(sortIdx);
    
    % 只显示前 20 个，避免刷屏
    dispCount = min(length(sortedCandidates), 20);
    fprintf('%-6s %-6s %-10s %-8s %-10s\n', 'W', 'D', 'Type', 'Value', 'Index');
    for k = 1:dispCount
        c = sortedCandidates(k);
        fprintf('%-6d %-6d %-10s %-8.4f %-10d\n', c.W, c.D, c.Type, c.Val, c.Idx);
    end
    
    % 绘制最佳结果
    best = sortedCandidates(1);
    W = best.W; D = best.D;
    
    fprintf('\n绘图展示最佳结果: W=%d, D=%d\n', W, D);
    
    % 重新计算最佳 W 的曲线用于绘图
    rx_delayed = x(1+D:end);
    rx_base    = x(1:end-D);
    P_metric = filter(ones(1, W), 1, conj(rx_base) .* rx_delayed);
    R_energy = filter(ones(1, W), 1, abs(rx_base).^2);
    M_best = real(P_metric ./ (R_energy + 1e-10));
    M_best = [M_best; zeros(D, 1)];
    
    figure('Name', 'Auto Search Result', 'Position', [100, 100, 1000, 600]);
    plot(t_axis, M_best, 'b'); hold on;
    xline(targetIndexThreshold, 'g--', 'Search Start (36300)');
    yline(corrThreshold, 'r--', 'Pos Thresh');
    yline(-corrThreshold, 'r--', 'Neg Thresh');
    
    % 标记找到的峰值
    plot(best.Idx, best.Val, 'ro', 'MarkerSize', 10, 'LineWidth', 2);
    text(best.Idx, best.Val, sprintf('  Val=%.3f\n  @%d', best.Val, best.Idx), 'Color', 'r');
    
    title(sprintf('最佳匹配结果: W=%d, D=%d, Index=%d', W, D, best.Idx));
    xlabel('Sample Index'); ylabel('Correlation (Real)');
    grid on;
    xlim([max(0, best.Idx - 2000), min(length(x)+startSample0, best.Idx + 2000)]);
end

%% ===== local helpers =====
function n = local_num_samples(inFile, numSamples, startSample0)
if ~isempty(numSamples)
    n = numSamples;
    return;
end
bytesPerComplexSample = 4;
info = dir(inFile);
if isempty(info)
    error('找不到文件: %s', inFile);
end
Ntotal = floor(info.bytes / bytesPerComplexSample);
if startSample0 >= Ntotal
    n = 0;
else
    n = Ntotal - startSample0;
end
end
