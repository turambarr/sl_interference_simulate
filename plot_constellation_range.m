% 独立脚本：绘制指定区间 IQ 数据的星座图
% 用途：观察特定时间段内的信号调制特性
% 输入：IQ 文件路径，起始样点，结束样点

clear; clc; close all;

%% 1. 参数设置
% inFile = 'test2.iq'; % 请修改您的文件名
inFile = 'sigtest1.iq'; 

% 设置要观察的样本区间 (0-based index)
% 例如：只想看刚才检测到的那个 "protrusion" 区域
regionStart = 14475; 
regionEnd   = 14475+6992; 

normalizeToUnit = true; % 根据文件格式一般选 true (int16 -> float[-1,1])
removeMean = false;      % 去除直流分量

% 降采样设置
enableResample = true;  % 开关：是否启用降采样
% 比例因子遍历：从 minPts 到 maxPts
% 你的 intent: int a=100; int b=6992; for i=0:1:b-a ...
% MATLAB 翻译: 遍历目标点数，每次增加一定步长 (例如 250)

% startPts = 2186;
startPts = 2186;
% endPts = 6992;
% stepPts = 1; 
% targetLengthList = startPts : stepPts : endPts;
targetLengthList = [2186]; % 用户指定单一长度

%% 1.1 读取全量背景数据 (用于宏观展示)
readBackground = true;
bgLoadLimit = 1e5; % 读取前10万点作为背景
x_bg = [];
if readBackground
    % 尝试读取背景数据 (需确保文件存在且 iq_read_int16_le 可用)
    try
        [x_bg_raw, ~] = iq_read_int16_le(inFile, 0, bgLoadLimit);
        x_bg = double(x_bg_raw);
        if normalizeToUnit
            x_bg = x_bg / 32768;
        end
    catch
        warning('无法读取背景数据，将跳过宏观展示');
    end
end

%% 2. 读取数据
numToRead = regionEnd - regionStart + 1;
if numToRead <= 0
    error('结束点必须大于起始点');
end

[x, meta] = iq_read_int16_le(inFile, regionStart, local_num_samples_limit(inFile, numToRead, regionStart));
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

fprintf('已读取区间 [%d, %d], 长度 %d\n', regionStart, regionStart+L-1, L);

% 备份原始数据用于时域绘图
x_orig = x;

%% 2.5 & 3. 批量降采样与绘图
% 模式：不再只画一张，而是遍历 targetLengthList 画多张图
% 每页 10 张图 (2x5 布局)

x_source = x; % 复数源数据
srcLen = length(x_source);

% 1. 如果有背景数据，画一个独立的 Overview 窗口 (Reference)
if exist('x_bg', 'var') && ~isempty(x_bg)
    figure('Name', 'Global Overview', 'Position', [50, 50, 800, 400]);
    t_bg = 0:length(x_bg)-1;
    plot(t_bg, abs(x_bg), 'Color', [0.7 0.7 0.7]); hold on;
    idx_start = regionStart + 1;
    idx_end = min(regionStart + L, length(x_bg));
    if idx_start < idx_end
        plot((idx_start:idx_end)-1, abs(x_bg(idx_start:idx_end)), 'b', 'LineWidth', 1.0);
    end
    xline(regionStart, 'r--'); xline(regionStart + L, 'r--');
    title('Time Domain Overview (Reference)');
    xlabel('Sample Index'); ylabel('Amplitude');
    grid on; axis tight;
end

% 计算统一的坐标轴范围 (基于原始数据最大幅度)
% maxVal = max(abs([real(x_source); imag(x_source)]));
% plotLimit = max(0.01, maxVal * 1.1); % 留 10% 边距
plotLimit = 0.2; % 用户指定固定范围
fprintf('使用固定坐标轴范围: [-%.4f, %.4f]\n', plotLimit, plotLimit);

% 2. 遍历并分页绘图
numPlots = length(targetLengthList);
plotsPerPage = 10;
numPages = ceil(numPlots / plotsPerPage);

fprintf('即将生成 %d 组页面，共 %d 个降采样测试用例...\n', numPages, numPlots);

for p = 1:numPages
    % 创建新 Figure
    figTitle = sprintf('Resample Scan: Page %d/%d (Start: %d pts)', p, numPages, targetLengthList((p-1)*plotsPerPage+1));
    figure('Name', figTitle, 'Position', [100, 100, 1400, 700]);
    
    startIdx = (p-1) * plotsPerPage + 1;
    endIdx = min(p * plotsPerPage, numPlots);
    
    subPlotIdx = 0;
    for k = startIdx : endIdx
        subPlotIdx = subPlotIdx + 1;
        tLen = targetLengthList(k); % 修正：定义 tLen

        % 增加调试信息以回答您的疑问
        fprintf('[Debug] k=%d, tLen (Target)=%d, srcLen (Source)=%d\n', k, round(tLen), srcLen);
        
        % 执行降采样
        if tLen >= srcLen
            x_curr = x_source;
            actualLen = srcLen;
            resampleNote = 'Original';
        else
            % 使用 resample (带抗混叠)
            % 确保 P, Q 为整数，防止 resample 异常
            P = round(tLen);
            Q = round(srcLen);
            
            % 再次检查防止 P, Q 异常
            if P < 1, P = 1; end
            if Q < 1, Q = 1; end
            
            try
                x_curr = resample(x_source, P, Q);
                fprintf('   -> resample(x, %d, %d) executed. Output Len: %d\n', P, Q, length(x_curr));
            catch ME
                warning('Resample failed: %s. Fallback to interp1.', ME.message);
                % Fallback (interp1)
                t_old = linspace(0, 1, srcLen);
                t_new = linspace(0, 1, P);
                x_curr = interp1(t_old, x_source, t_new, 'linear').';
            end
            
            actualLen = length(x_curr);
            resampleNote = 'Resampled';
        end
        
        % 绘图
        if numPlots > 1
             subplot(2, 5, subPlotIdx);
        else
             subplot(1, 1, 1);
        end
        
        % 奇偶着色法 (Odd/Even Coloring)
        % 红色=奇数索引点, 蓝色=偶数索引点
        % 如果是 2 倍过采样，红蓝点群通常会分离
        idx_odd = 1:2:length(x_curr);
        idx_even = 2:2:length(x_curr);
        
        plot(real(x_curr(idx_odd)), imag(x_curr(idx_odd)), 'r.', 'MarkerSize', 8, 'DisplayName', 'Odd');
        hold on;
        plot(real(x_curr(idx_even)), imag(x_curr(idx_even)), 'b.', 'MarkerSize', 8, 'DisplayName', 'Even');
        
        % 样式美化
        axis equal; grid on; box on;
        legend show;
        
        % 使用自动计算的统一范围
        xlim([-plotLimit plotLimit]); ylim([-plotLimit plotLimit]); 
        
        % 标题信息
        ratio = actualLen / srcLen;
        title(sprintf('N=%d (Ratio=%.2f)', actualLen, ratio), 'FontSize', 10);
        
        % 只在左侧加 Label 节省空间
        if mod(subPlotIdx, 5) == 1
            ylabel('Q');
            if subPlotIdx > 5, xlabel('I'); end
        end
    end
    
    % 交互式等待
    if p < numPages
        promptStr = sprintf('>> [交互] 已显示第 %d/%d 页。按 Enter 继续，输入 q 退出: ', p, numPages);
        userIn = input(promptStr, 's');
        
        if strcmpi(userIn, 'q') || strcmpi(userIn, 'exit')
            fprintf('>> 用户终止遍历。\n');
            break;
        end
        
        close(gcf); % 关闭当前页，准备显示下一页
    else
        fprintf('>> [交互] 所有页面显示完毕。\n');
    end
end

% 不再执行下面的单图逻辑了，直接结束脚本部分
return; 

% 下面是旧代码，已被 return 跳过，保留或删除均可 (这里直接覆盖掉旧代码)
%% ===== local helpers =====
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

function [x, meta] = iq_read_int16_le(filename, startSample, numSamples, headerBytes)
    if nargin < 4, headerBytes = 0; end
    if nargin < 3, numSamples = []; end
    
    fid = fopen(filename, 'rb');
    if fid < 0
        error('无法打开文件: %s', filename);
    end
    
    % Seek
    status = fseek(fid, headerBytes + startSample * 4, 'bof');
    if status ~= 0
        fclose(fid);
        error('定位文件失败 (startSample 可能越界)');
    end
    
    % Read
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
    
    % Interleaved I/Q
    i_data = raw(1:2:end);
    q_data = raw(2:2:end);
    
    % Handle odd length (shouldn't happen for valid IQ files but for safety)
    len = min(length(i_data), length(q_data));
    x = complex(i_data(1:len), q_data(1:len));
    
    meta.numSamplesRead = len;
end
