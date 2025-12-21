function out = estimate_ofdm_cp_locations(inFile, N, Ng, startSample0, endSample0, Fs, opts)
%ESTIMATE_OFDM_CP_LOCATIONS 通过CP相关度量定位OFDM符号边界（估计功能）
%
% out = estimate_ofdm_cp_locations(inFile, N, Ng, startSample0, endSample0, Fs, opts)
%
% 约定：
% - 输入为纯IQ数据文件（int16 little-endian，I/Q交织），复采样点索引为 0-based。
% - 输出的 symCpStart0 / symDataStart0 均为 0-based 复采样点索引。
%
% 输入：
% - inFile       : 纯数据IQ文件路径
% - N            : 有效符号长度（FFT点数）
% - Ng           : CP长度
% - startSample0 : 搜索范围起点（0-based，包含）
% - endSample0   : 搜索范围终点（0-based，包含）；可传 [] 表示到文件尾
% - Fs           : 采样率（Hz，可传 []，仅用于 out.diag）
% - opts         : 可选结构体
%   .normalizeToUnit      (default true)
%   .removeMean           (default true)
%   .refineRadius         (default round((N+Ng)/4))
%   .maxSymbolsToReport   (default 5000)
%   .numAnchorsToTry      (default 20)
%   .anchorMinSeparation  (default round((N+Ng)/2))
%   .debugSymIndex        (default [])  % 若设置，打印该符号的 M(CP)
%   .topK                 (default 20)  % 诊断：打印全局 topK
%   .verbose              (default true)
%
% 输出 out：
% - symCpStart0     : CP起点（0-based）
% - symDataStart0   : 数据起点（0-based）= symCpStart0 + Ng
% - diag            : 诊断结构体（M、分位数、top峰等）

if nargin < 7 || isempty(opts)
    opts = struct();
end
if nargin < 6
    Fs = [];
end
if nargin < 5
    endSample0 = [];
end
if nargin < 4 || isempty(startSample0)
    startSample0 = 0;
end

opts = local_defaults(opts, N, Ng);

bytesPerComplexSample = 4;
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

Lsym = N + Ng;

Nneed = endSample0 - startSample0 + 1;
[x, meta] = iq_read_int16_le(inFile, startSample0, Nneed);
if meta.numSamplesRead < Nneed
    warning('实际读取%d点，小于期望%d点（已到文件尾）。', meta.numSamplesRead, Nneed);
end
x = double(x);
if opts.normalizeToUnit
    x = x / 32768;
end
if opts.removeMean
    x = x - mean(x);
end

L = numel(x);
if L < (N + Ng + 1)
    error('数据太短：至少需要 N+Ng+1=%d 点，当前只有 %d 点。', N+Ng+1, L);
end

% CP相关度量 M(d)
prod = x(1:end-N) .* conj(x(1+N:end));
P = conv(prod, ones(Ng,1), 'valid');
ref = abs(x(1+N:end)).^2;
R = conv(ref, ones(Ng,1), 'valid');
M = (abs(P).^2) ./ (R.^2 + eps);

d0_rel = (0:numel(M)-1).';

% 诊断：分布
Mmax = max(M);
Mmed = median(M);
Mmean = mean(M);
Mq95 = prctile(M, 95);
Mq99 = prctile(M, 99);
Mq999 = prctile(M, 99.9);

% 选锚点：多个 top 峰
[MsSorted, ordSorted] = sort(M, 'descend');
anchors = zeros(0,1);
for ii = 1:numel(ordSorted)
    cand = ordSorted(ii) - 1; % 0-based rel
    if isempty(anchors) || all(abs(anchors - cand) >= opts.anchorMinSeparation)
        anchors(end+1,1) = cand; %#ok<AGROW>
    end
    if numel(anchors) >= opts.numAnchorsToTry
        break;
    end
end
if isempty(anchors)
    error('无法生成锚点（anchors为空）');
end

minRel = 0;
maxRel = d0_rel(end);

bestScore = -Inf;
bestFirstRel = anchors(1);
bestSymRel = [];

for ai = 1:numel(anchors)
    firstRel = anchors(ai);

    kMin = ceil((minRel - firstRel) / Lsym);
    kMax = floor((maxRel - firstRel) / Lsym);
    ks = kMin:kMax;
    expected = firstRel + ks(:) * Lsym;

    symRelTry = zeros(size(expected));
    for i = 1:numel(expected)
        c = expected(i);
        lo = max(minRel, c - opts.refineRadius);
        hi = min(maxRel, c + opts.refineRadius);
        seg = M(lo+1:hi+1);
        [~, im] = max(seg);
        symRelTry(i) = (lo + (im-1));
    end

    symRelTry = unique(symRelTry);
    if isempty(symRelTry)
        continue;
    end

    score = mean(M(symRelTry + 1));
    if score > bestScore
        bestScore = score;
        bestFirstRel = firstRel;
        bestSymRel = symRelTry;
    end
end

symRel = bestSymRel;
if isempty(symRel)
    error('未能得到有效的 symRel（可能数据/参数不匹配）。');
end

symCpStart0 = startSample0 + symRel;
symDataStart0 = symCpStart0 + Ng;

% 诊断：symRel上的度量
Msym = M(symRel + 1);

% 诊断：topK
K = min(opts.topK, numel(M));
topOrd = ordSorted(1:K);
topD0 = startSample0 + (topOrd - 1);

if opts.verbose
    fprintf('M(d)统计：max=%.6g, median=%.6g, mean=%.6g, max/median=%.3g\n', Mmax, Mmed, Mmean, Mmax/(Mmed+eps));
    fprintf('M(d)分位数：p95=%.6g, p99=%.6g, p99.9=%.6g\n', Mq95, Mq99, Mq999);
    fprintf('锚点选择：尝试%d个候选，选中 firstRel=%d（0-based rel），score=%.6g\n', numel(anchors), bestFirstRel, bestScore);
    fprintf('定位到 %d 个OFDM符号（CP起点，0-based）\n', numel(symCpStart0));
    showN = min(numel(symCpStart0), opts.maxSymbolsToReport);
    for i = 1:showN
        fprintf('sym%04d: CP=%d, DATA=%d\n', i, symCpStart0(i), symDataStart0(i));
    end
    if numel(symCpStart0) > showN
        fprintf('...（已截断，只显示前%d个）\n', showN);
    end

    fprintf('M@symRel统计：min=%.6g, median=%.6g, mean=%.6g, max=%.6g\n', min(Msym), median(Msym), mean(Msym), max(Msym));

    fprintf('M(d) Top-%d（全局）：\n', K);
    for ii = 1:K
        fprintf('%2d) d0=%d, M=%.6g\n', ii, topD0(ii), MsSorted(ii));
    end
    if K >= 2
        gapsTop = diff(sort(topD0));
        fprintf('Top峰间隔统计：min=%d, median=%.2f, max=%d（期望≈%d）\n', min(gapsTop), median(double(gapsTop)), max(gapsTop), Lsym);
    end

    if ~isempty(opts.debugSymIndex)
        si = opts.debugSymIndex;
        if si >= 1 && si <= numel(symRel)
            fprintf('debug sym%04d: CP=%d, M(CP)=%.6g\n', si, symCpStart0(si), M(symRel(si)+1));
        end
    end

    if numel(symCpStart0) >= 2
        gaps = diff(symCpStart0);
        fprintf('符号间隔统计：min=%d, max=%d, median=%.2f（期望=%d）\n', min(gaps), max(gaps), median(double(gaps)), Lsym);
    end
end

out = struct();
out.symCpStart0 = symCpStart0;
out.symDataStart0 = symDataStart0;
out.diag = struct();
out.diag.Fs = Fs;
out.diag.N = N;
out.diag.Ng = Ng;
out.diag.Lsym = Lsym;
out.diag.startSample0 = startSample0;
out.diag.endSample0 = endSample0;
out.diag.M = M;
out.diag.d0_rel = d0_rel;
out.diag.Mstats = struct('max', Mmax, 'median', Mmed, 'mean', Mmean, 'p95', Mq95, 'p99', Mq99, 'p999', Mq999, 'maxOverMedian', Mmax/(Mmed+eps));
out.diag.anchorsRel = anchors;
out.diag.bestFirstRel = bestFirstRel;
out.diag.bestScore = bestScore;
out.diag.symRel = symRel;
out.diag.Msym = Msym;
out.diag.topD0 = topD0;
out.diag.topM = MsSorted(1:K);
end

function opts = local_defaults(opts, N, Ng)
if ~isfield(opts, 'normalizeToUnit');     opts.normalizeToUnit = true; end
if ~isfield(opts, 'removeMean');          opts.removeMean = true; end
if ~isfield(opts, 'refineRadius');        opts.refineRadius = round((N+Ng)/4); end
if ~isfield(opts, 'maxSymbolsToReport');  opts.maxSymbolsToReport = 5000; end
if ~isfield(opts, 'numAnchorsToTry');     opts.numAnchorsToTry = 20; end
if ~isfield(opts, 'anchorMinSeparation'); opts.anchorMinSeparation = round((N+Ng)/2); end
if ~isfield(opts, 'debugSymIndex');       opts.debugSymIndex = []; end
if ~isfield(opts, 'topK');                opts.topK = 20; end
if ~isfield(opts, 'verbose');             opts.verbose = true; end
end
