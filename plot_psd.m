% 独立脚本：绘制IQ信号功率谱密度(PSD)
% 已知：Fs=409.6MHz，带宽=320MHz，中心频率=635MHz
% 数据格式：int16 小端序，I/Q交织（纯数据文件，无文件头）

clear; clc;

%% 参数区（按需修改）
inFile = 'test1.iq';

Fs = 409.6e6;      % 采样率 Hz
Fc = 635e6;        % 中心频率 Hz
showAbsoluteFreq = true; % true: 横轴显示 Fc+f；false: 显示基带频率 f

startSample = 0;   % 从第几个“复采样点”开始（0-based）
Nread = 2e6;       % 读取多少复采样点用于PSD（越大越平滑，但更慢）

normalizeToUnit = true; % true: int16 / 32768 归一化到近似[-1,1]
removeMean = false;      % true: 去直流

% Welch参数
segLen = 65536;     % 每段长度（复采样点数）
overlapRatio = 0.5; % 重叠比例 [0,1)
nfft = 1024;       % FFT点数（>=segLen）

%% 读取数据
[x, meta] = iq_read_int16_le(inFile, startSample, Nread);
N = meta.numSamplesRead;
if N == 0
    error('未读取到数据。');
end

x = double(x);
if normalizeToUnit
    x = x / 32768;
end
if removeMean
    x = x - mean(x);
end

%% Welch PSD（自实现，避免工具箱依赖）
% 若数据太短，自动调整 Welch 参数，避免 segLen>N 直接报错
[segLenEff, nfftEff] = sanitize_welch_params(N, segLen, nfft);
if segLenEff ~= segLen
    fprintf('NOTE: segLen=%d > N=%d, auto set segLen=%d\n', segLen, N, segLenEff);
end
if nfftEff ~= nfft
    fprintf('NOTE: nfft=%d < segLen=%d or too small, auto set nfft=%d\n', nfft, segLenEff, nfftEff);
end

[psd, f] = welch_psd_centered(x, Fs, segLenEff, overlapRatio, nfftEff);
psd_db = 10*log10(psd + eps);

if showAbsoluteFreq
    f_plot = (f + Fc) / 1e6; % MHz
    xlab = 'Frequency (MHz)';
else
    f_plot = f / 1e6; % MHz
    xlab = 'Baseband Frequency (MHz)';
end

%% 绘图
figure('Name', sprintf('PSD: %s', inFile));
plot(f_plot, psd_db, 'LineWidth', 1);
grid on;
xlabel(xlab);
ylabel('PSD (dB/Hz)');
title(sprintf('Welch PSD (start=%d, N=%d, seg=%d, nfft=%d)', startSample, N, segLenEff, nfftEff));

% 交互放大/查看
zoom on;
pan on;


function [Pxx, f] = welch_psd_centered(x, Fs, segLen, overlapRatio, nfft)
%WELCH_PSD_CENTERED 简单Welch PSD（双边、fftshift后中心化）
% x: 复数序列
% Fs: 采样率
% segLen: 分段长度
% overlapRatio: 重叠比例
% nfft: FFT点数

x = x(:);
N = numel(x);

if segLen < 2
    error('segLen 必须 >= 2（当前 segLen=%d）', segLen);
end

if segLen <= 0 || nfft < segLen
    error('segLen需要>0，且nfft需要>=segLen');
end
if overlapRatio < 0 || overlapRatio >= 1
    error('overlapRatio 必须在 [0,1)');
end

hop = max(1, round(segLen * (1 - overlapRatio)));

% Hann窗（工具箱无依赖）
n = (0:segLen-1).';
w = 0.5 - 0.5*cos(2*pi*n/(segLen-1));
U = sum(w.^2); % 用于PSD归一化

acc = zeros(nfft, 1);
numSeg = 0;

for start = 1:hop:(N - segLen + 1)
    seg = x(start:start+segLen-1);
    seg = seg .* w;

    X = fft(seg, nfft);
    acc = acc + (abs(X).^2);
    numSeg = numSeg + 1;
end

if numSeg == 0
    error('数据长度不足以形成一个分段：N=%d, segLen=%d', N, segLen);
end

% 平均功率谱，再换算为功率谱密度(每Hz)
% 对应关系：Pxx = E{|X|^2} / (Fs * U)
Pxx = (acc / numSeg) / (Fs * U);

% 频率轴（中心化）
Pxx = fftshift(Pxx);
f = ((-nfft/2):(nfft/2-1)).' * (Fs / nfft);
end

function [segLenEff, nfftEff] = sanitize_welch_params(N, segLen, nfft)
%SANITIZE_WELCH_PARAMS 让 Welch 参数适配当前数据长度 N
% - 当 segLen > N 时，将 segLen 降到 2^floor(log2(N))（至少 256；若 N<256 则取 N）
% - 保证 nfft >= segLen，且 nfft 为 2 的幂（便于 FFT）
segLenEff = segLen;
nfftEff = nfft;

if N < segLenEff
    if N >= 256
        segLenEff = 2^floor(log2(N));
    else
        segLenEff = max(2, N);
    end
end

if nfftEff < segLenEff
    nfftEff = 2^nextpow2(segLenEff);
end

% 也把 nfft 拉到 2 的幂（如果用户填了非 2 幂）
if nfftEff ~= 2^nextpow2(nfftEff)
    nfftEff = 2^nextpow2(nfftEff);
end
end
