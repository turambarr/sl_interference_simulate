% 独立脚本：扫描降采样倍率 D 并批量绘制星座图
% 用途：在 D 周围精细扫描，寻找最佳星座图
% 交互：每页显示 6 张，回车翻页，q 退出

clear; clc; close all;

%% 1. 参数设置
inFile = 'sigtest1.iq'; 


% 采样范围设置 (尽量与原先保持一致或根据需要调整)
regionStart = 14475 + 181; 
L_read = 874;      % 读取长度
% regionEnd = regionStart + L_read; 

% D 扫描设置
D_center = 3.199;
% 扫描半径：比如 +/- 0.005，覆盖 0.01 的范围
% 100 张图，步长 = 0.01 / 100 = 0.0001
D_half_range = 0.005; 
num_steps = 100;

D_values = linspace(D_center - D_half_range, D_center + D_half_range, num_steps);

normalizeToUnit = true; 
removeMean = true;      

%% 2. 读取数据
[x, meta] = local_iq_read(inFile, regionStart, local_num_samples_limit(inFile, L_read, regionStart));
L = meta.numSamplesRead;
if L == 0
    error('未读取到数据，请检查起始点是否越界');
end

x = double(x);
if normalizeToUnit
    x = x / 32768;
end
if removeMean
    x = x - mean(x);
end

fprintf('已读取数据: %s, 区间[%d, %d], 长度 %d\n', inFile, regionStart, regionStart+L-1, L);
fprintf('准备扫描 D: %.6f ~ %.6f (中心=%.4f, 共 %d 点, 步长≈%.6f)\n', ...
    min(D_values), max(D_values), D_center, num_steps, D_values(2)-D_values(1));
fprintf('提示：弹出的图形窗口中将每页显示 6 张图。\n');
fprintf('需要先点击一下 Matlab 命令窗口使其获得焦点，然后按 Enter 翻页，输入 q 退出。\n');
fprintf('按任意键开始...\n');
pause;

%% 3. 循环扫描绘图
hFig = figure('Name', 'Constellation Scan', 'NumberTitle', 'off', 'Color', 'w');
set(hFig, 'Position', [100, 100, 1200, 800]); % 调整窗口大小

plots_per_page = 6;
num_pages = ceil(num_steps / plots_per_page);

for idx = 1:num_steps
    D = D_values(idx);
    
    % 执行重采样
    % D 意味着 InputRate / OutputRate = D
    % resample(x, P, Q) -> OutputRate = InputRate * P/Q
    % 所以 P/Q Should be approx 1/D
    [P, Q] = rat(1 / D, 1e-7); 
    
    % resample 可能会引入滤波器暂态，这里数据短，可能会有影响，但只能如此
    x_resampled = resample(x, P, Q);
    
    % 子图位置计算
    page_idx = mod(idx-1, plots_per_page) + 1;
    subplot(2, 3, page_idx);
    
    % 绘图
    plot(real(x_resampled), imag(x_resampled), '.', 'MarkerSize', 6);
    axis square; grid on;
    title(sprintf('D = %.6f\n(P=%d, Q=%d)', D, P, Q), 'FontSize', 10, 'FontWeight', 'bold');
    
    % 固定坐标轴以便人眼对比
    limit = 0.5; % 根据归一化后的数据范围，通常在 0.5 左右
    xlim([-limit, limit]);
    ylim([-limit, limit]);
    
    % 翻页判定
    if page_idx == plots_per_page || idx == num_steps
        sgtitle(sprintf('Scanning: %d / %d (Page %d/%d)', idx, num_steps, ceil(idx/plots_per_page), num_pages));
        drawnow;
        
        command = input(sprintf('[Page %d] Press Enter for next, "q" to quit: ', ceil(idx/plots_per_page)), 's');
        
        if strcmpi(command, 'q')
            fprintf('用户请求退出。\n');
            break;
        end
        
        if idx < num_steps
            clf(hFig); % 清空当前窗口
        end
    end
end
fprintf('扫描结束。\n');

%% Local Helper Functions
function n = local_num_samples_limit(inFile, reqSamples, startSample0)
    bytesPerComplexSample = 4;
    info = dir(inFile);
    if isempty(info)
        error('找不到文件: %s', inFile);
    end
    Ntotal = floor(info.bytes / bytesPerComplexSample);
    
    if startSample0 >= Ntotal
        n = 0;
    else
        available = Ntotal - startSample0;
        n = min(reqSamples, available);
    end
end

function [x, meta] = local_iq_read(filename, startSample, numSamples, startBytes)
    if nargin < 4, startBytes = 0; end
    if nargin < 3, numSamples = []; end
    
    fid = fopen(filename, 'rb');
    if fid < 0
        error('无法打开文件: %s', filename);
    end
    
    status = fseek(fid, startBytes + startSample * 4, 'bof');
    if status ~= 0
        fclose(fid);
        error('定位文件失败');
    end
    
    if isempty(numSamples)
        raw = fread(fid, Inf, 'int16');
    else
        raw = fread(fid, numSamples * 2, 'int16');
    end
    fclose(fid);
    
    if isempty(raw)
        x = [];
        meta.numSamplesRead = 0;
        return;
    end
    
    i_data = raw(1:2:end);
    q_data = raw(2:2:end);
    len = min(length(i_data), length(q_data));
    x = complex(i_data(1:len), q_data(1:len));
    meta.numSamplesRead = len;
end
