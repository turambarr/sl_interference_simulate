% 独立脚本：计算 sigtest1.iq 和 sigtest2.iq 的互相关
% 用途：分析两个 IQ 信号之间的时间延迟 (Lag) 和相似度
% 注意：如果 sigtest2.iq 不存在，请确保文件已放置在同目录下

clear; clc; close all;

%% 1. 参数设置
file1 = 'sigtest1.iq';
file2 = 'sigtest2.iq';

% 读取最大样本数 (防止内存溢出)
% 如果文件很大，建议限制长度，或者分段处理
maxSamples = 1e6; % 读取前 100 万个点

%% 2. 读取数据
fprintf('正在读取文件 1: %s ...\n', file1);
[x1, meta1] = iq_read_wrapper(file1, maxSamples);

fprintf('正在读取文件 2: %s ...\n', file2);
[x2, meta2] = iq_read_wrapper(file2, maxSamples);

if isempty(x1) || isempty(x2)
    error('读取失败：其中一个文件为空或不存在。');
end

fprintf('文件 1 长度: %d\n', length(x1));
fprintf('文件 2 长度: %d\n', length(x2));

%% 3. 计算互相关 (Cross-Correlation)
% 计算 x1 和 x2 的互相关 (直接使用原始数据，不归一化)
fprintf('正在计算互相关...\n');

[c, lags] = xcorr(x1, x2,'unbiased');

fprintf('\n=== 计算完成 ===\n');
% 不再自动搜索峰值，直接绘图

%% 4. 绘图展示
figure('Position', [100, 100, 1000, 600], 'Name', 'Cross-Correlation (Raw)');

% 实部
subplot(2, 1, 1);
plot(lags, real(c), 'b');
title('Cross-Correlation (Real Part)');
xlabel('Lag'); ylabel('Amplitude');
grid on; axis tight;

% 虚部
subplot(2, 1, 2);
plot(lags, imag(c), 'r');
title('Cross-Correlation (Imaginary Part)');
xlabel('Lag'); ylabel('Amplitude');
grid on; axis tight;

%% 5. 辅助函数：封装读取逻辑
function [x, meta] = iq_read_wrapper(filename, numSamples)
    if ~isfile(filename)
        warning('文件不存在: %s', filename);
        x = [];
        meta.numSamplesRead = 0;
        return;
    end

    % 调用内部读取函数
    [x, meta] = iq_read_int16_le(filename, 0, numSamples);
    
    % 转 double 并归一化
    x = double(x) / 32768;
    x = x - mean(x); % 去直流
end

% 原始读取函数 (从 plot_constellation.m 复制)
function [x, meta] = iq_read_int16_le(filename, startSample, numSamples)
    if nargin < 3, numSamples = []; end
    
    fid = fopen(filename, 'rb');
    if fid < 0
        error('无法打开文件: %s', filename);
    end
    
    status = fseek(fid, startSample * 4, 'bof');
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
