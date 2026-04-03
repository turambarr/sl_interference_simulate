% 独立脚本：直接绘制IQ信号频谱（单次FFT）
% 已知：Fs=409.6MHz，带宽=320MHz，中心频率=635MHz
% 数据格式：int16 小端序，I/Q交织（纯数据文件，无文件头）

clear; clc;

%% 参数区（按需修改）
inFile = 'sigtest6.iq';

Fs = 409.6e6;      % 采样率 Hz
Fc = -63.5e6;        % 中心频率 Hz
showAbsoluteFreq = true; % true: 横轴显示 Fc+f；false: 显示基带频率 f

startSample = 0;   % 从第几个“复采样点”开始（0-based）
Nread = 5e6;       % 读取多少复采样点用于PSD（越大越平滑，但更慢）

% 按“起止点”选取信号段（0-based，含 endPoint）
useStartEndRange = false;  % true: 使用 startPoint/endPoint；false: 使用 startSample/Nread
startPoint = 11562;            % 选段起点（复采样点）
endPoint = 14472 ;       % 选段终点（复采样点，含该点）

normalizeToUnit = true; % true: int16 / 32768 归一化到近似[-1,1]
removeMean = false;      % true: 去直流

% 直接FFT参数
nfft = 2^18;            % FFT点数（建议2的幂；可设[]表示自动=2^nextpow2(N)）
useHannWindow = false;  % true: 加Hann窗；false: 矩形窗（更“直接”）

%% 选段参数处理
if useStartEndRange
    if endPoint < startPoint
        error('endPoint 必须 >= startPoint（当前 startPoint=%d, endPoint=%d）', startPoint, endPoint);
    end
    startSampleEff = round(startPoint);
    NreadEff = round(endPoint - startPoint + 1);
else
    startSampleEff = round(startSample);
    NreadEff = round(Nread);
end

if startSampleEff < 0
    error('startSample/startPoint 必须 >= 0（当前=%d）', startSampleEff);
end
if NreadEff <= 0
    error('Nread 或 (endPoint-startPoint+1) 必须 > 0（当前=%d）', NreadEff);
end

%% 文件长度与读区间检查（避免 fseek 失败）
fInfo = dir(inFile);
if isempty(fInfo)
    error('文件不存在: %s', inFile);
end
bytesPerComplexSample = 4; % int16 I + int16 Q
totalSamples = floor(fInfo.bytes / bytesPerComplexSample);
if totalSamples <= 0
    error('文件长度异常，无法读取复采样点: %s', inFile);
end

if startSampleEff >= totalSamples
    error(['起始点超出文件范围：start=%d, total=%d。\n' ...
           '请减小 startPoint/startSample（最大允许起点=%d）。'], ...
          startSampleEff, totalSamples, totalSamples-1);
end

maxReadable = totalSamples - startSampleEff;
if NreadEff > maxReadable
    fprintf(['NOTE: 请求读取长度超出文件末尾，自动截断 Nread: %d -> %d ' ...
             '(start=%d, total=%d)\n'], ...
            NreadEff, maxReadable, startSampleEff, totalSamples);
    NreadEff = maxReadable;
end

%% 读取数据
[x, meta] = iq_read_int16_le(inFile, startSampleEff, NreadEff);
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

%% 直接FFT频谱
if isempty(nfft) || nfft <= 0
    nfftEff = 2^nextpow2(N);
else
    nfftEff = 2^nextpow2(round(nfft));
end

x = x(:);
if useHannWindow
    n = (0:N-1).';
    w = 0.5 - 0.5*cos(2*pi*n/(N-1));
else
    w = ones(N,1);
end

X = fftshift(fft(x .* w, nfftEff));
spec_db = 20*log10(abs(X) + eps);
spec_db = spec_db - max(spec_db); % 归一化到峰值0 dB，便于观察

f = ((-nfftEff/2):(nfftEff/2-1)).' * (Fs / nfftEff);

if showAbsoluteFreq
    f_plot = (f + Fc) / 1e6; % MHz
    xlab = 'Frequency (MHz)';
else
    f_plot = f / 1e6; % MHz
    xlab = 'Baseband Frequency (MHz)';
end

%% 绘图
figure('Name', sprintf('Spectrum: %s', inFile));
plot(f_plot, spec_db, 'LineWidth', 1);
grid on;
xlabel(xlab);
ylabel('Magnitude (dB, normalized)');
title(sprintf('Direct FFT Spectrum (start=%d, end=%d, N=%d, nfft=%d)', ...
    startSampleEff, startSampleEff + N - 1, N, nfftEff));

% 交互放大/查看
zoom on;
pan on;
