% 独立脚本：对 test2.iq 做“整段自相关”（一次性对全信号计算）
% 要求：仅自相关；横轴=延迟(lag)，纵轴=自相关值（不做任何滑窗/热力图/其它分析）
%
% 线性自相关（非负 lag）：
%   r[k] = sum_{n=0}^{L-1-k} x[n+k] * conj(x[n])
% 归一化：
%   rho[k] = r[k] / r[0]

clear; clc;

%% 参数区（按需修改）
inFile = 'test2.iq';

startSample0 = 0;   % 0-based
numSamples = [];    % [] 表示读到文件尾（文件很大时建议给一个有限值）

normalizeToUnit = true; % int16/32768
removeMean = true;      % 去直流

% 只画 0..maxLag；设为 [] 表示画到 L-1
maxLag = [];

plotDb = false;          % true: 20log10(|rho|)；false: |rho|

%% 读取
[x, meta] = iq_read_int16_le(inFile, startSample0, local_num_samples(inFile, numSamples, startSample0));
L = meta.numSamplesRead;
if L == 0
    error('未读取到数据。');
end
x = double(x);
if normalizeToUnit
    x = x / 32768;
end
if removeMean
    x = x - mean(x);
end

%% 整段自相关（非负 lag）
% 用 FFT 计算线性自相关：ifft( FFT(x).*conj(FFT(x)) ) 在足够零填充下，前 L 项等价
if isempty(maxLag)
    maxLag = L - 1;
else
    maxLag = min(maxLag, L - 1);
end

nfft = 2^nextpow2(2*L - 1);
X = fft(x(:), nfft);
r = ifft(X .* conj(X));
r = r(1:maxLag+1);

r0 = r(1);
rho = r / (r0 + eps);
lags = (0:maxLag).';

fprintf('=== Autocorrelation (test2, full) ===\n');
fprintf('file=%s, startSample0=%d, L=%d\n', inFile, startSample0, L);
fprintf('removeMean=%d, normalizeToUnit=%d\n', removeMean, normalizeToUnit);
fprintf('maxLag=%d\n', maxLag);

figure('Name', 'Autocorr |rho| (test2, full)');
if plotDb
    plot(lags, 20*log10(abs(rho) + eps), 'LineWidth', 1);
    ylabel('|\rho(lag)| (dB, normalized to lag=0)');
else
    plot(lags, abs(rho), 'LineWidth', 1);
    ylabel('|\rho(lag)|');
end
grid on;
xlabel('Lag (samples)');
title(sprintf('Autocorrelation (full), file=%s, L=%d', inFile, L));

zoom on; pan on; datacursormode on;

%% ===== local helpers =====
function n = local_num_samples(inFile, numSamples, startSample0)
if ~isempty(numSamples)
    n = numSamples;
    return;
end
bytesPerComplexSample = 4;
info = dir(inFile);
if isempty(info)
    error('找不到文件: %s', inFile);
end
Ntotal = floor(info.bytes / bytesPerComplexSample);
if startSample0 >= Ntotal
    n = 0;
else
    n = Ntotal - startSample0;
end
end
