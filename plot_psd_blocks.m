% 独立脚本：分块绘制功率谱密度(PSD)
% 需求：在 [startSample, endSample] 范围内，每50000点画一张PSD，最后全部显示。
%
% 约定：startSample/endSample 为“复采样点”索引（0-based），输入为纯数据（无文件头）
% 数据格式：int16 小端序，I/Q交织

clear; clc;

%% 参数区（按需修改）
inFile = '20250912222305_part1.iq';

Fs = 409.6e6;      % Hz
Fc = 635e6;        % Hz
showAbsoluteFreq = false; % true: Fc+f；false: 基带f

startSample = 800000;          % 0-based（复采样点）
endSample = 1850000;   % 0-based（复采样点，包含端点）
blockLen = 1000000;         % 每块复采样点数（按需求固定50000）

normalizeToUnit = true;   % int16/32768
removeMean = true;        % 去直流

% Welch参数（会自动裁剪到 <= blockLen）
segLen = 8192;
overlapRatio = 0.5;
nfft = 8192;

% 显示布局：每个figure放多少张子图（太大会很挤）
plotsPerFigure = 16;

%% 计算块数
if endSample < startSample
    error('endSample 必须 >= startSample');
end

numAvail = endSample - startSample + 1;
numBlocks = floor(numAvail / blockLen);
if numBlocks <= 0
    error('范围内样点不足一个block：需要至少%d点，当前只有%d点。', blockLen, numAvail);
end

fprintf('总范围点数=%d，blockLen=%d，块数=%d\n', numAvail, blockLen, numBlocks);

% Welch参数安全化
segLenEff = min(segLen, blockLen);
nfftEff = max(nfft, segLenEff);
if nfftEff < segLenEff
    nfftEff = segLenEff;
end

%% 循环分块画PSD
for b = 1:numBlocks
    % 当前块在原始信号中的起点（0-based）
    blkStart = startSample + (b-1)*blockLen;

    % 读取这一块
    [x, meta] = iq_read_int16_le(inFile, blkStart, blockLen);
    if meta.numSamplesRead < blockLen
        warning('第%d块读取不足（%d/%d），停止。', b, meta.numSamplesRead, blockLen);
        break;
    end

    x = double(x);
    if normalizeToUnit
        x = x / 32768;
    end
    if removeMean
        x = x - mean(x);
    end

    % PSD
    [psd, f] = welch_psd_centered(x, Fs, segLenEff, overlapRatio, nfftEff);
    psd_db = 10*log10(psd + eps);

    if showAbsoluteFreq
        f_plot = (f + Fc) / 1e6;
        xlab = 'Frequency (MHz)';
    else
        f_plot = f / 1e6;
        xlab = 'Baseband Frequency (MHz)';
    end

    % 分页显示
    figIndex = floor((b-1) / plotsPerFigure) + 1;
    posInFig = mod((b-1), plotsPerFigure) + 1;

    figure(figIndex);
    set(gcf, 'Name', sprintf('PSD Blocks (page %d)', figIndex));

    if exist('tiledlayout', 'file') == 2
        % 用持久化 tiledlayout（每个figure第一次创建）
        tl = getappdata(gcf, 'tl');
        if isempty(tl) || ~isvalid(tl)
            rows = ceil(sqrt(plotsPerFigure));
            cols = ceil(plotsPerFigure / rows);
            tl = tiledlayout(rows, cols, 'TileSpacing', 'compact', 'Padding', 'compact');
            setappdata(gcf, 'tl', tl);
        end
        nexttile(posInFig);
    else
        rows = ceil(sqrt(plotsPerFigure));
        cols = ceil(plotsPerFigure / rows);
        subplot(rows, cols, posInFig);
    end

    plot(f_plot, psd_db, 'LineWidth', 1);
    grid on;
    xlabel(xlab);
    ylabel('PSD (dB/Hz)');
    title(sprintf('blk %d: [%d..%d]', b, blkStart, blkStart+blockLen-1));
end

% 打开交互
for k = 1:figIndex
    figure(k);
    zoom on;
    pan on;
    datacursormode on;
end


function [Pxx, f] = welch_psd_centered(x, Fs, segLen, overlapRatio, nfft)
% 简单Welch PSD（双边、fftshift后中心化）

x = x(:);
N = numel(x);

if segLen <= 0 || nfft < segLen
    error('segLen需要>0，且nfft需要>=segLen');
end
if overlapRatio < 0 || overlapRatio >= 1
    error('overlapRatio 必须在 [0,1)');
end
if N < segLen
    error('数据长度不足：N=%d, segLen=%d', N, segLen);
end

hop = max(1, round(segLen * (1 - overlapRatio)));

% Hann窗（无工具箱依赖）
n = (0:segLen-1).';
w = 0.5 - 0.5*cos(2*pi*n/(segLen-1));
U = sum(w.^2);

acc = zeros(nfft, 1);
numSeg = 0;

for start = 1:hop:(N - segLen + 1)
    seg = x(start:start+segLen-1);
    seg = seg .* w;
    X = fft(seg, nfft);
    acc = acc + (abs(X).^2);
    numSeg = numSeg + 1;
end

Pxx = (acc / numSeg) / (Fs * U);
Pxx = fftshift(Pxx);
f = ((-nfft/2):(nfft/2-1)).' * (Fs / nfft);
end
