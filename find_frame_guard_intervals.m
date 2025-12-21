function result = find_frame_guard_intervals(inFile, startSample0, endSample0, opts)
%FIND_FRAME_GUARD_INTERVALS 寻找帧与帧之间的保护间隔（Guard Interval）
%
% result = find_frame_guard_intervals(inFile, startSample0, endSample0, opts)
%
% 适用场景：真实信号中“成帧突发（burst）”之间存在一段明显低能量的空隙（保护间隔）。
% 本函数通过短时能量包络 + 自适应阈值，检测 burst 与 gap。
%
% 约定：
% - 输入文件为纯IQ数据（int16 little-endian，I/Q交织），复采样点索引为 0-based。
% - 输出的起止点均为 0-based 复采样点索引。
%
% 输入：
% - inFile       : 纯数据IQ文件路径
% - startSample0 : 搜索起点（0-based，包含）
% - endSample0   : 搜索终点（0-based，包含）；可传 [] 表示到文件尾
% - opts         : 可选结构体
%   .normalizeToUnit (default true)    int16/32768
%   .removeMean      (default true)    去直流（有利于能量阈值稳定）
%   .winLen          (default 4096)    能量窗口长度（复采样点数）
%   .hop             (default 1024)    能量步进（复采样点数）
%   .smoothLen       (default 9)       对能量序列做移动平均平滑（以“bin”为单位）
%   .thrMethod       (default 'mad')   'mad' 或 'quantile'
%   .thrK            (default 8)       thr = median(E)+k*mad(E)
%   .thrQuantile     (default 0.90)    quantile 方法阈值分位数
%   .minBurstLen     (default 20000)   最短帧长度（复采样点）
%   .minGuardLen     (default 20000)   最短保护间隔长度（复采样点）
%   .chunkSamples    (default 2e6)     分块读取长度（复采样点）
%   .makePlot        (default false)   是否绘制能量与检测结果
%   .verbose         (default true)    打印摘要
%
% 输出 result：
% - framesStart0 / framesEnd0 : 检测到的帧区间（burst）
% - guardsStart0 / guardsEnd0 : 帧间保护间隔（gap）
% - energyBins : 每个能量bin对应的中心/末端样点索引（0-based）
% - energy     : 平滑后的能量序列
% - threshold  : 使用的阈值
%
% 注意：能量窗会引入 ~winLen 的时间模糊，边界是近似定位；如果你需要更精确边界，可在结果附近做二次精细化。

if nargin < 4 || isempty(opts)
    opts = struct();
end
if nargin < 3
    error('需要参数: inFile, startSample0, endSample0(可空), opts(可空)');
end

opts = local_defaults(opts);

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
    error('startSample0/endSample0 不合法');
end

% 为了能量计算，确保 chunk 至少 winLen
opts.chunkSamples = max(opts.chunkSamples, opts.winLen + opts.hop);

% 全局能量bin序列（逐块拼接）
energy = zeros(0,1);
energyBins0 = zeros(0,1); % 每个bin对应的“窗末端样点索引”(0-based)

rangeLen = endSample0 - startSample0 + 1;
numReadSoFar = 0;
prevTail = complex(zeros(0,1));
prevTailStart0 = startSample0;

while numReadSoFar < rangeLen
    thisStart0 = startSample0 + numReadSoFar;
    thisN = min(opts.chunkSamples, rangeLen - numReadSoFar);

    [x, meta] = iq_read_int16_le(inFile, thisStart0, thisN);
    if meta.numSamplesRead == 0
        break;
    end
    x = double(x);
    if opts.normalizeToUnit
        x = x / 32768;
    end
    if opts.removeMean
        x = x - mean(x);
    end

    % 拼接上一块尾巴，保证跨块窗口连续
    if ~isempty(prevTail)
        xCat = [prevTail; x];
        catStart0 = prevTailStart0;
    else
        xCat = x;
        catStart0 = thisStart0;
    end

    % 计算短时能量（用 filter 实现滑动和，避免巨型 conv）
    e = abs(xCat).^2;
    acc = filter(ones(opts.winLen,1), 1, e);
    valid = acc(opts.winLen:end) / opts.winLen; % 对应窗末端从 (winLen-1) 到 end

    % 以 hop 抽取能量bin
    Ebins = valid(1:opts.hop:end);

    % 每个 bin 的“窗末端”全局索引（0-based）
    % valid(1) 对应窗末端 sample = catStart0 + (opts.winLen-1)
    k = (0:numel(Ebins)-1).';
    bins0 = catStart0 + (opts.winLen-1) + k*opts.hop;

    energy = [energy; Ebins]; %#ok<AGROW>
    energyBins0 = [energyBins0; bins0]; %#ok<AGROW>

    % 留下尾巴给下一块：至少 winLen-1 点
    keep = min(numel(xCat), opts.winLen-1);
    prevTail = xCat(end-keep+1:end);
    prevTailStart0 = catStart0 + (numel(xCat) - keep);

    numReadSoFar = numReadSoFar + meta.numSamplesRead;
    if meta.numSamplesRead < thisN
        break;
    end
end

if isempty(energy)
    error('未能生成能量序列（energy为空）。');
end

% 平滑能量（按bin）
energySmooth = local_smooth(energy, opts.smoothLen);

% 阈值
thr = local_threshold(energySmooth, opts);

% 活动区（burst）
isActive = energySmooth > thr;

% 连通段提取 + 形态学清理（按最小长度）
[burstRuns, gapRuns] = local_runs(isActive);

minBurstBins = ceil(opts.minBurstLen / opts.hop);
minGuardBins = ceil(opts.minGuardLen / opts.hop);

% 丢弃太短 burst
burstRuns = burstRuns((burstRuns(:,2) - burstRuns(:,1) + 1) >= minBurstBins, :);

% 重新生成 gap（基于 burstRuns）
if isempty(burstRuns)
    framesStart0 = zeros(0,1);
    framesEnd0 = zeros(0,1);
    guardsStart0 = zeros(0,1);
    guardsEnd0 = zeros(0,1);
else
    % burstRuns 已按出现顺序
    framesStart0 = local_bins_to_start0(energyBins0, burstRuns(:,1), opts);
    framesEnd0 = local_bins_to_end0(energyBins0, burstRuns(:,2));

    % 保护间隔 = 相邻两帧之间的 gap
    gS = framesEnd0(1:end-1) + 1;
    gE = framesStart0(2:end) - 1;
    keepG = (gE - gS + 1) >= opts.minGuardLen;
    guardsStart0 = gS(keepG);
    guardsEnd0 = gE(keepG);
end

% 裁剪到用户范围
framesStart0 = max(framesStart0, startSample0);
framesEnd0 = min(framesEnd0, endSample0);

guardsStart0 = max(guardsStart0, startSample0);
guardsEnd0 = min(guardsEnd0, endSample0);

result = struct();
result.inFile = inFile;
result.startSample0 = startSample0;
result.endSample0 = endSample0;
result.opts = opts;
result.energyBins0 = energyBins0;
result.energy = energySmooth;
result.threshold = thr;
result.framesStart0 = framesStart0;
result.framesEnd0 = framesEnd0;
result.guardsStart0 = guardsStart0;
result.guardsEnd0 = guardsEnd0;

if opts.verbose
    fprintf('=== Guard Interval Detection ===\n');
    fprintf('range=[%d..%d], bins=%d, winLen=%d, hop=%d\n', startSample0, endSample0, numel(energyBins0), opts.winLen, opts.hop);
    fprintf('thr=%.6g, frames=%d, guards=%d\n', thr, numel(framesStart0), numel(guardsStart0));
    fprintf('energy: min=%.6g, median=%.6g, max=%.6g, activeBins=%d/%d\n', ...
        min(energySmooth), median(energySmooth), max(energySmooth), sum(energySmooth > thr), numel(energySmooth));
    if ~isempty(guardsStart0)
        gLen = guardsEnd0 - guardsStart0 + 1;
        fprintf('guardLen: min=%d, median=%.1f, max=%d (samples)\n', min(gLen), median(double(gLen)), max(gLen));
    end
end

if opts.makePlot
    figure('Name', 'Frame/Guard Detection (Energy)');
    plot(energyBins0, 10*log10(energySmooth + eps), 'b'); hold on;
    yThr = 10*log10(thr + eps);
    yline(yThr, 'r--', 'thr');
    grid on;
    xlabel('Sample Index (0-based)');
    ylabel('Energy (dB)');
    title('Short-time Energy + Threshold');

    % 标注帧与保护间隔
    yl = ylim;
    for i = 1:numel(framesStart0)
        x1 = framesStart0(i); x2 = framesEnd0(i);
        patch([x1 x2 x2 x1], [yl(1) yl(1) yl(2) yl(2)], [0.8 0.9 1.0], 'FaceAlpha', 0.15, 'EdgeColor', 'none');
    end
    for i = 1:numel(guardsStart0)
        x1 = guardsStart0(i); x2 = guardsEnd0(i);
        patch([x1 x2 x2 x1], [yl(1) yl(1) yl(2) yl(2)], [0.9 0.8 0.8], 'FaceAlpha', 0.12, 'EdgeColor', 'none');
    end
    uistack(findobj(gca,'Type','line'),'top');

    zoom on; pan on; datacursormode on;
end

end

function opts = local_defaults(opts)
if ~isfield(opts, 'normalizeToUnit'); opts.normalizeToUnit = true; end
if ~isfield(opts, 'removeMean');      opts.removeMean = true; end
if ~isfield(opts, 'winLen');          opts.winLen = 4096; end
if ~isfield(opts, 'hop');             opts.hop = 1024; end
if ~isfield(opts, 'smoothLen');       opts.smoothLen = 9; end
if ~isfield(opts, 'thrMethod');       opts.thrMethod = 'mad'; end
if ~isfield(opts, 'thrK');            opts.thrK = 8; end
if ~isfield(opts, 'thrQuantile');     opts.thrQuantile = 0.90; end
if ~isfield(opts, 'minBurstLen');     opts.minBurstLen = 20000; end
if ~isfield(opts, 'minGuardLen');     opts.minGuardLen = 20000; end
if ~isfield(opts, 'chunkSamples');    opts.chunkSamples = 2e6; end
if ~isfield(opts, 'makePlot');        opts.makePlot = false; end
if ~isfield(opts, 'verbose');         opts.verbose = true; end

% basic checks
if opts.winLen <= 0 || opts.hop <= 0
    error('winLen/hop 必须为正');
end
if opts.hop > opts.winLen
    % 允许，但提示用户能量时间分辨率更粗
end
if opts.smoothLen < 1
    opts.smoothLen = 1;
end
end

function y = local_smooth(x, smoothLen)
if smoothLen <= 1
    y = x;
    return;
end
x = x(:);
w = ones(smoothLen,1) / smoothLen;
y = filter(w, 1, x);
% 对齐：前 smoothLen-1 点滤波尚未“充满”，简单用原值替代
if numel(y) >= smoothLen
    y(1:smoothLen-1) = x(1:smoothLen-1);
end
end

function thr = local_threshold(E, opts)
E = E(:);
switch lower(opts.thrMethod)
    case 'mad'
        med = median(E);
        % 与注释一致：MAD = median(|E - median(E)|)
        % 不依赖 Statistics Toolbox 的 mad()
        m = median(abs(E - med));
        thr = med + opts.thrK * m;
    case 'quantile'
        q = opts.thrQuantile;
        q = max(0, min(1, q));
        thr = local_quantile(E, q);
    otherwise
        error('thrMethod 仅支持 mad / quantile');
end

% 避免阈值过小
thr = max(thr, eps);
end

function v = local_quantile(x, q)
% 线性插值分位数（不依赖 Statistics Toolbox）
x = x(:);
if isempty(x)
    v = NaN;
    return;
end
x = sort(x);
n = numel(x);
if q <= 0
    v = x(1);
    return;
end
if q >= 1
    v = x(end);
    return;
end
pos = 1 + q * (n - 1); % 1-based
lo = floor(pos);
hi = ceil(pos);
if lo == hi
    v = x(lo);
else
    frac = pos - lo;
    v = x(lo) * (1 - frac) + x(hi) * frac;
end
end

function [runs1, runs0] = local_runs(mask)
% 返回 mask==1 的连通段 runs1，以及 mask==0 的连通段 runs0
mask = logical(mask(:));
if isempty(mask)
    runs1 = zeros(0,2);
    runs0 = zeros(0,2);
    return;
end

d = diff([false; mask; false]);
starts1 = find(d == 1);
ends1 = find(d == -1) - 1;
runs1 = [starts1 ends1];

mask0 = ~mask;
d0 = diff([false; mask0; false]);
starts0 = find(d0 == 1);
ends0 = find(d0 == -1) - 1;
runs0 = [starts0 ends0];
end

function s0 = local_bins_to_start0(energyBins0, binStartIdx, opts)
% bin 对应窗末端索引，换算为窗起点近似索引：end - (winLen-1)
energyBins0 = energyBins0(:);
binStartIdx = binStartIdx(:);
end0 = energyBins0(binStartIdx);
s0 = end0 - (opts.winLen - 1);
end

function e0 = local_bins_to_end0(energyBins0, binEndIdx)
energyBins0 = energyBins0(:);
binEndIdx = binEndIdx(:);
e0 = energyBins0(binEndIdx);
end
