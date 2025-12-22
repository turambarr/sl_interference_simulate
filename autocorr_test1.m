% 独立脚本：对 test1.iq 做自相关（FFT实现，无工具箱）
% 目的：观察在 lag=N（以及 lag=Lsym=N+Ng）处是否存在明显相关峰
%
% 输入约定：纯数据IQ（int16 little-endian，I/Q交织，无文件头），0-based 复采样点索引

clear; clc;

%% 参数区（按需修改）
inFile = 'test1.iq';

% 如果只想取 test1 的一部分，可在这里设置（0-based，针对 test1 文件内部）
startSample0 = 0;
numSamples = []; % [] 表示读到文件尾

% OFDM 参考参数（用于在图上标注 lag）
N = 1024;
Ng = 192;
Lsym = N + Ng;

normalizeToUnit = true; % int16/32768
removeMean = true;      % 自相关前去均值（去直流）

% 绘图范围（lag 轴）
maxLag = 4000; % 只画 0..maxLag；设为 [] 表示画到 L-1

% 只看峰：打印前 topK 个相关峰（排除 lag=0）
topKPeaks = 20;

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

if isempty(maxLag)
    maxLag = L - 1;
else
    maxLag = min(maxLag, L - 1);
end

%% 自相关（非负 lag）
% 线性自相关 r[k] = sum_{n=0}^{L-1-k} x[n+k] * conj(x[n])
% 用 FFT 计算：ifft( FFT(x).*conj(FFT(x)) ) 在足够零填充下，前 L 项等价于线性自相关
nfft = 2^nextpow2(2*L - 1);
X = fft(x, nfft);
r = ifft(X .* conj(X));
r = r(1:maxLag+1);

% 归一化：rho[k] = r[k]/r[0]
r0 = r(1);
rho = r / (r0 + eps);

lags = (0:maxLag).';

%% 打印关键 lag
fprintf('=== Autocorrelation (test1) ===\n');
fprintf('file=%s, startSample0=%d, L=%d\n', inFile, startSample0, L);
fprintf('removeMean=%d, normalizeToUnit=%d\n', removeMean, normalizeToUnit);
fprintf('r0=%.6g + j%.6g, |r0|=%.6g\n', real(r0), imag(r0), abs(r0));

local_print_lag(rho, N, 'N');
local_print_lag(rho, Lsym, 'Lsym');

local_print_top_peaks(rho, topKPeaks);

%% 绘图
figure('Name', sprintf('Autocorr |rho| (test1)'));
plot(lags, 20*log10(abs(rho) + eps), 'LineWidth', 1);
grid on;
xlabel('Lag (samples)');
ylabel('|\rho(lag)| (dB, normalized to lag=0)');
title(sprintf('Autocorrelation (L=%d), N=%d, Ng=%d, Lsym=%d', L, N, Ng, Lsym));

hold on;
yl = ylim;
if N <= maxLag
    xline(N, 'r--', 'N');
end
if Lsym <= maxLag
    xline(Lsym, 'm--', 'Lsym');
end
ylim(yl);

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

function local_print_lag(rho, lag, name)
if lag <= 0
    return;
end
if lag + 1 > numel(rho)
    fprintf('%s lag=%d 超出当前 maxLag\n', name, lag);
    return;
end
v = rho(lag+1);
fprintf('%s lag=%d: |rho|=%.6g, dB=%.2f\n', name, lag, abs(v), 20*log10(abs(v) + eps));
end

function local_print_top_peaks(rho, topK)
absRho = abs(rho);
if numel(absRho) <= 1
    return;
end
absRhoNoZero = absRho(2:end); % lag=1..end
[vals, ord] = sort(absRhoNoZero, 'descend');
K = min(topK, numel(vals));
fprintf('Top-%d peaks (exclude lag=0):\n', K);
for i = 1:K
    lag = ord(i); % 因为 absRhoNoZero(1) 对应 lag=1
    v = vals(i);
    fprintf('%2d) lag=%d: |rho|=%.6g, dB=%.2f\n', i, lag, v, 20*log10(v + eps));
end
end
