function out = estimate_main(action, varargin)
%ESTIMATE_MAIN 估计功能主入口（不包含验证类功能）
%
% out = estimate_main(action, Name,Value,...)
%
% action（选择要执行的估计功能）：
% - 'ofdm_cp'         : OFDM 符号边界定位（基于 CP 相关度量）
% - 'guard_intervals' : 帧间保护间隔（burst gap）检测
%
% 通用约定：
% - 输入 IQ 文件为纯数据（int16 little-endian，I/Q 交织），索引为 0-based 复采样点。
% - 本入口只“做估计/检测”，不会把验证函数（例如切片自证）装进来。
%
% 示例：
%   % 1) OFDM CP 定位
%   out = estimate_main('ofdm_cp', ...
%       'inFile','20250912222305_part1_cut2.iq', 'N',1024, 'Ng',192, ...
%       'startSample0',0, 'endSample0',[], 'Fs',409.6e6, 'verbose',true);
%
%   % 2) 帧间保护间隔
%   out = estimate_main('guard_intervals', ...
%       'inFile','20250912222305_part1_cut2.iq', 'startSample0',0, 'endSample0',[], ...
%       'winLen',4096, 'hop',1024, 'minGuardLen',20000, 'makePlot',true);

if nargin < 1
    error('需要 action');
end

switch lower(string(action))
    case "ofdm_cp"
        out = local_ofdm_cp(varargin{:});

    case "guard_intervals"
        out = local_guard_intervals(varargin{:});

    otherwise
        error('未知 action=%s。支持：ofdm_cp, guard_intervals', string(action));
end

end

function out = local_ofdm_cp(varargin)
% 参数
p = inputParser;
p.addParameter('inFile', '', @(s)ischar(s) || isstring(s));
p.addParameter('N', [], @(x)isnumeric(x) && isscalar(x) && x>0);
p.addParameter('Ng', [], @(x)isnumeric(x) && isscalar(x) && x>0);
p.addParameter('startSample0', 0, @(x)isnumeric(x) && isscalar(x) && x>=0);
p.addParameter('endSample0', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x>=0));
p.addParameter('Fs', [], @(x)isnumeric(x) && isscalar(x) && (isempty(x) || x>0));

% opts 映射
p.addParameter('normalizeToUnit', true, @(x)islogical(x) || isnumeric(x));
p.addParameter('removeMean', true, @(x)islogical(x) || isnumeric(x));
p.addParameter('refineRadius', [], @(x)isnumeric(x) && isscalar(x) && x>=0);
p.addParameter('maxSymbolsToReport', 5000, @(x)isnumeric(x) && isscalar(x) && x>=0);
p.addParameter('numAnchorsToTry', 20, @(x)isnumeric(x) && isscalar(x) && x>=1);
p.addParameter('anchorMinSeparation', [], @(x)isnumeric(x) && isscalar(x) && x>=0);
p.addParameter('debugSymIndex', [], @(x)isnumeric(x) && isscalar(x) && x>=1);
p.addParameter('topK', 20, @(x)isnumeric(x) && isscalar(x) && x>=1);
p.addParameter('verbose', true, @(x)islogical(x) || isnumeric(x));
p.parse(varargin{:});

inFile = string(p.Results.inFile);
N = p.Results.N;
Ng = p.Results.Ng;
if strlength(inFile) == 0 || isempty(N) || isempty(Ng)
    error('ofdm_cp 需要参数 inFile, N, Ng');
end

opts = struct();
opts.normalizeToUnit = logical(p.Results.normalizeToUnit);
opts.removeMean = logical(p.Results.removeMean);
if ~isempty(p.Results.refineRadius); opts.refineRadius = p.Results.refineRadius; end
opts.maxSymbolsToReport = p.Results.maxSymbolsToReport;
opts.numAnchorsToTry = p.Results.numAnchorsToTry;
if ~isempty(p.Results.anchorMinSeparation); opts.anchorMinSeparation = p.Results.anchorMinSeparation; end
if ~isempty(p.Results.debugSymIndex); opts.debugSymIndex = p.Results.debugSymIndex; end
opts.topK = p.Results.topK;
opts.verbose = logical(p.Results.verbose);

out = estimate_ofdm_cp_locations(char(inFile), N, Ng, p.Results.startSample0, p.Results.endSample0, p.Results.Fs, opts);
end

function out = local_guard_intervals(varargin)
% 参数
p = inputParser;
p.addParameter('inFile', '', @(s)ischar(s) || isstring(s));
p.addParameter('startSample0', 0, @(x)isnumeric(x) && isscalar(x) && x>=0);
p.addParameter('endSample0', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x>=0));

% opts 映射（与 find_frame_guard_intervals.m 一致）
p.addParameter('normalizeToUnit', true, @(x)islogical(x) || isnumeric(x));
p.addParameter('removeMean', true, @(x)islogical(x) || isnumeric(x));
p.addParameter('winLen', 4096, @(x)isnumeric(x) && isscalar(x) && x>0);
p.addParameter('hop', 1024, @(x)isnumeric(x) && isscalar(x) && x>0);
p.addParameter('smoothLen', 9, @(x)isnumeric(x) && isscalar(x) && x>=1);
p.addParameter('thrMethod', 'mad', @(s)ischar(s) || isstring(s));
p.addParameter('thrK', 8, @(x)isnumeric(x) && isscalar(x) && x>=0);
p.addParameter('thrQuantile', 0.90, @(x)isnumeric(x) && isscalar(x) && x>=0 && x<=1);
p.addParameter('minBurstLen', 20000, @(x)isnumeric(x) && isscalar(x) && x>=0);
p.addParameter('minGuardLen', 20000, @(x)isnumeric(x) && isscalar(x) && x>=0);
p.addParameter('chunkSamples', 2e6, @(x)isnumeric(x) && isscalar(x) && x>0);
p.addParameter('makePlot', false, @(x)islogical(x) || isnumeric(x));
p.addParameter('verbose', true, @(x)islogical(x) || isnumeric(x));

p.parse(varargin{:});

inFile = string(p.Results.inFile);
if strlength(inFile) == 0
    error('guard_intervals 需要参数 inFile');
end

opts = struct();
opts.normalizeToUnit = logical(p.Results.normalizeToUnit);
opts.removeMean = logical(p.Results.removeMean);
opts.winLen = p.Results.winLen;
opts.hop = p.Results.hop;
opts.smoothLen = p.Results.smoothLen;
opts.thrMethod = char(p.Results.thrMethod);
opts.thrK = p.Results.thrK;
opts.thrQuantile = p.Results.thrQuantile;
opts.minBurstLen = p.Results.minBurstLen;
opts.minGuardLen = p.Results.minGuardLen;
opts.chunkSamples = p.Results.chunkSamples;
opts.makePlot = logical(p.Results.makePlot);
opts.verbose = logical(p.Results.verbose);

out = find_frame_guard_intervals(char(inFile), p.Results.startSample0, p.Results.endSample0, opts);
end
