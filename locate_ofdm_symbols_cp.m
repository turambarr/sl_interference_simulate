% 独立脚本：用循环前缀(CP)相关定位OFDM符号边界
% 已知：输入为纯数据（无文件头），int16小端序，I/Q交织
% 参数：子载波数 N=1024，循环前缀长度 Ng=192
%
% 输出：
% - symCpStart0 : 每个OFDM符号CP起点（0-based复采样点索引）
% - symDataStart0: 每个OFDM符号有效数据起点（0-based）= symCpStart0 + Ng

clear; clc;

%% 参数区（按需修改）
inFile = '20250912222305_part1_cut2.iq';

Fs = 409.6e6;   % Hz（仅用于可选显示）
N = 1024;       % 有效符号长度（FFT点数）
Ng = 192;       % CP长度
Lsym = N + Ng;  % 一个OFDM符号总长度

% 要搜索的范围（0-based复采样点索引，包含端点）；留空则用全文件
startSample0 = 0;
endSample0 = []; % 例如：startSample0+Lsym*200-1

% 归一化/去直流（通常有利）
normalizeToUnit = true; % int16/32768
removeMean = true;

% 从全局最大峰开始建立“符号栅格”，每隔Lsym在附近微调
refineRadius = round(Lsym/4); % 每个期望位置±refineRadius内找局部最大
maxSymbolsToReport = 5000;    % 防止打印太多

% 鲁棒性：真实信号里全局最大峰可能是杂峰。
% 这里尝试多个 top 峰作为“锚点”，用生成的 symRel 处 M 值均值打分，选最优锚点。
numAnchorsToTry = 20;
anchorMinSeparation = round(Lsym/2);

% 诊断：如果你想看某个符号的度量值（比如与 run_segment_similarity.m 的 idx 对齐）
debugSymIndex = 140;

%% 读取数据（按文件大小自动读取范围）
bytesPerComplexSample = 4; % int16 I + int16 Q
info = dir(inFile);
if isempty(info)
    error('找不到文件: %s', inFile);
end
Ntotal = floor(info.bytes / bytesPerComplexSample);

if isempty(endSample0)
    endSample0 = Ntotal - 1;
end
if startSample0 < 0 || endSample0 < startSample0
    error('startSample0/endSample0 设置不合法');
end

Nneed = endSample0 - startSample0 + 1;
[x, meta] = iq_read_int16_le(inFile, startSample0, Nneed);
if meta.numSamplesRead < Nneed
    warning('实际读取%d点，小于期望%d点（已到文件尾）。', meta.numSamplesRead, Nneed);
end
x = double(x);
if normalizeToUnit
    x = x / 32768;
end
if removeMean
    x = x - mean(x);
end

%% 计算CP相关度量 M(d)
% P(d) = sum_{k=0..Ng-1} x[d+k] * conj(x[d+k+N])
% R(d) = sum_{k=0..Ng-1} |x[d+k+N]|^2
% M(d) = |P(d)|^2 / (R(d)^2 + eps)

L = numel(x);
if L < (N + Ng + 1)
    error('数据太短：至少需要 N+Ng+1=%d 点，当前只有 %d 点。', N+Ng+1, L);
end

prod = x(1:end-N) .* conj(x(1+N:end));
P = conv(prod, ones(Ng,1), 'valid');
ref = abs(x(1+N:end)).^2;
R = conv(ref, ones(Ng,1), 'valid');

M = (abs(P).^2) ./ (R.^2 + eps);

% 诊断：度量整体“尖锐程度”
Mmax = max(M);
Mmed = median(M);
Mmean = mean(M);
fprintf('M(d)统计：max=%.6g, median=%.6g, mean=%.6g, max/median=%.3g\n', Mmax, Mmed, Mmean, Mmax/(Mmed+eps));

% 诊断：分位数（看尾部是否显著）
Mq95 = prctile(M, 95);
Mq99 = prctile(M, 99);
Mq999 = prctile(M, 99.9);
fprintf('M(d)分位数：p95=%.6g, p99=%.6g, p99.9=%.6g\n', Mq95, Mq99, Mq999);

% M对应的候选CP起点 d（相对读取窗口、0-based）
d0_rel = (0:numel(M)-1).';

%% 选锚点：尝试多个 top 峰
% 选出分离度足够的一组峰位置（避免一坨相邻点）
[MsSorted, ordSorted] = sort(M, 'descend');
anchors = zeros(0,1);
for ii = 1:numel(ordSorted)
    cand = ordSorted(ii) - 1; % 转 0-based rel
    if isempty(anchors) || all(abs(anchors - cand) >= anchorMinSeparation)
        anchors(end+1,1) = cand; %#ok<AGROW>
    end
    if numel(anchors) >= numAnchorsToTry
        break;
    end
end
if isempty(anchors)
    error('无法生成锚点（anchors为空）');
end

bestScore = -Inf;
bestFirstRel = anchors(1);
bestSymRel = [];

minRel = 0;
maxRel = d0_rel(end);

for ai = 1:numel(anchors)
    firstRel = anchors(ai);

    % 生成期望栅格：向前/向后扩展
    kMin = ceil((minRel - firstRel) / Lsym);
    kMax = floor((maxRel - firstRel) / Lsym);
    ks = kMin:kMax;
    expected = firstRel + ks(:) * Lsym;

    % 对每个期望位置在附近找局部最大
    symRelTry = zeros(size(expected));
    for i = 1:numel(expected)
        c = expected(i);
        lo = max(minRel, c - refineRadius);
        hi = min(maxRel, c + refineRadius);
        seg = M(lo+1:hi+1); % +1转MATLAB索引
        [~, im] = max(seg);
        symRelTry(i) = (lo + (im-1));
    end

    symRelTry = unique(symRelTry);
    if isempty(symRelTry)
        continue;
    end

    % 用 symRel 上的 M 值打分：均值越大，说明“整列符号边界”越一致
    score = mean(M(symRelTry + 1));
    if score > bestScore
        bestScore = score;
        bestFirstRel = firstRel;
        bestSymRel = symRelTry;
    end
end

firstRel = bestFirstRel;
symRel = bestSymRel;
fprintf('锚点选择：尝试%d个候选，选中 firstRel=%d（0-based rel），score=%.6g\n', numel(anchors), firstRel, bestScore);

%% 按最优锚点得到的 symRel 输出

% 转成“全文件0-based索引”
symCpStart0 = startSample0 + symRel;
symDataStart0 = symCpStart0 + Ng;

% 诊断：打印指定符号的 M 值（在 symRel 上取样）
if ~isempty(symRel) && debugSymIndex >= 1 && debugSymIndex <= numel(symRel)
    fprintf('debug sym%04d: CP=%d, M(CP)=%.6g\n', debugSymIndex, symCpStart0(debugSymIndex), M(symRel(debugSymIndex)+1));
end

% 诊断：symRel 处的 M 值整体质量
if ~isempty(symRel)
    Msym = M(symRel + 1);
    fprintf('M@symRel统计：min=%.6g, median=%.6g, mean=%.6g, max=%.6g\n', min(Msym), median(Msym), mean(Msym), max(Msym));
end

% 诊断：全局 top-K 峰位置（不依赖 findpeaks 工具箱）
topK = 20;
[Ms, ord] = sort(M, 'descend');
ord = ord(1:min(topK, numel(ord)));
topD0 = startSample0 + (ord - 1); % 0-based
fprintf('M(d) Top-%d（全局）：\n', numel(ord));
for ii = 1:numel(ord)
    fprintf('%2d) d0=%d, M=%.6g\n', ii, topD0(ii), Ms(ii));
end
if numel(ord) >= 2
    gapsTop = diff(sort(topD0));
    fprintf('Top峰间隔统计：min=%d, median=%.2f, max=%d（期望≈%d）\n', min(gapsTop), median(double(gapsTop)), max(gapsTop), Lsym);
end

%% 打印结果（截断）
fprintf('定位到 %d 个OFDM符号（CP起点，0-based）\n', numel(symCpStart0));
showN = min(numel(symCpStart0), maxSymbolsToReport);
for i = 1:showN
    fprintf('sym%04d: CP=%d, DATA=%d\n', i, symCpStart0(i), symDataStart0(i));
end
if numel(symCpStart0) > showN
    fprintf('...（已截断，只显示前%d个）\n', showN);
end

%% 可视化：度量曲线 + 标记符号起点
figure('Name', 'OFDM CP Correlation Metric');
plot(startSample0 + d0_rel, M, 'LineWidth', 1);
hold on;
plot(symCpStart0, M(symRel+1), 'r.', 'MarkerSize', 12);
grid on;
xlabel('Sample Index (0-based)');
ylabel('M(d)');
title(sprintf('CP metric (N=%d, Ng=%d, Lsym=%d)', N, Ng, Lsym));

zoom on;
pan on;
datacursormode on;

% 可选：打印符号间隔统计（看是否稳定为Lsym）
if numel(symCpStart0) >= 2
    gaps = diff(symCpStart0);
    fprintf('符号间隔统计：min=%d, max=%d, median=%.2f（期望=%d）\n', min(gaps), max(gaps), median(double(gaps)), Lsym);
end
