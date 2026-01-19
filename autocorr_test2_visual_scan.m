% 独立脚本：Autocorr Visual Scanner
% 目标：完全仿照 autocorr_test2_slide128 的处理方式（W=D, 计算 M_real），
% 遍历 W=64:1024。每 10 个 W 画一屏（每屏 10 张子图），方便肉眼筛选。

clear; clc; close all;

%% 1. 基础配置
inFile = 'test2.iq';

startSample0 = 0;   
numSamples = 200000; % 为了画图速度，建议只取前 20w 点即可看清头部结构
                     
normalizeToUnit = true; 
removeMean = true;      

% 扫描范围
scanRange = 64:1200; % 64 到 1024，步长 1

% 分页显示设置
plotsPerPage = 10;   % 每张 Figure 显示多少个 W 的结果

%% 2. 读取数据 (只读一次)
% 放在循环外，避免重复 I/O
[x, meta] = iq_read_int16_le(inFile, startSample0, local_num_samples(inFile, numSamples, startSample0));
L = meta.numSamplesRead;
if L == 0, error('未读取到数据。'); end
x = double(x);
if normalizeToUnit, x = x / 32768; end
if removeMean, x = x - mean(x); end

t_axis = (0:length(x)-1) + startSample0;

fprintf('数据已加载 (L=%d)。准备开始可视化扫描...\n', length(x));
fprintf('按 Ctrl+C 可随时终止。\n');

%% 3. 遍历与绘图
numTotal = length(scanRange);
numPages = ceil(numTotal / plotsPerPage);

for p = 1:numPages
    % 当前页要显示的 W 索引范围
    idxStart = (p-1) * plotsPerPage + 1;
    idxEnd   = min(p * plotsPerPage, numTotal);
    
    currentBatchW = scanRange(idxStart:idxEnd);
    
    % 创建新 Figure (全屏/大窗口)
    figName = sprintf('Scan Batch %d/%d (W: %d ~ %d)', p, numPages, currentBatchW(1), currentBatchW(end));
    f = figure('Name', figName, 'NumberTitle', 'off', 'Position', [50, 50, 1200, 900]);
    
    % 循环画子图
    for k = 1:length(currentBatchW)
        W = currentBatchW(k);
        D = W; % 假定 Delay = Window
        
        % --- 核心算法 (同 slide128) ---
        if length(x) <= D
            continue; 
        end

        rx_delayed = x(1+D:end);
        rx_base    = x(1:end-D);
        
        conj_prod = conj(rx_base) .* rx_delayed;
        
        % 滑动求和
        P_metric = filter(ones(1, W), 1, conj_prod);
        
        % 能量
        rx_power = abs(rx_base).^2;
        R_energy = filter(ones(1, W), 1, rx_power);
        
        % 计算实部度量
        M_complex = P_metric ./ (R_energy + 1e-10);
        M_real = real(M_complex);
        
        % 补齐长度
        M_plot = [M_real; zeros(D, 1)];
        
        % --- 绘图 ---
        % 计算子图位置
        subplot(plotsPerPage/2, 2, k); % 比如 5行2列 = 10张图
        
        plot(t_axis, M_plot, 'b', 'LineWidth', 1);
        hold on;
        yline(0, 'k-');
        yline(-0.5, 'r:', 'Possible Hit');
        
        ylim([-1.1, 1.1]);
        xlim([t_axis(1), t_axis(min(end, 20000))]); % 默认只显示前 20000 点，看头部即可
        title(sprintf('W = D = %d', W));
        grid on;
    end
    
    drawnow;
    
    % 交互暂停：询问是否继续
    fprintf('当前显示第 %d 页 (W=%d ~ %d)。\n', p, currentBatchW(1), currentBatchW(end));
    userIn = input('输入 Enter 继续下一页，输入 q 退出: ', 's');
    if strcmpi(userIn, 'q')
        fprintf('用户终止扫描。\n');
        break;
    end
    
    % 关闭旧图以防内存溢出（或者您可以选择保留，把 close f 注释掉）
    close(f); 
end

fprintf('扫描结束。\n');

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
