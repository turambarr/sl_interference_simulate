% 独立脚本：对 test2.iq 做“滑动自相关结构检测”
% 目标：利用 Schmidl & Cox 算法原理寻找 PSS 类似的重复结构
% 参考：pss_test.m

clear; clc;

%% 参数区（按需修改）
inFile = 'test2.iq';

startSample0 = 0;   % 0-based
numSamples = [];    % [] 表示读到文件尾（如果文件很大，建议改成一个有限值）

normalizeToUnit = true; % int16/32768
removeMean = true;      % 去直流

W = 219;                % 窗长（固定 128）

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

if L < W
    error('数据太短：L=%d < W=%d', L, W);
end

%% 3. 核心算法: 滑动自相关 (Sliding Autocorrelation)
% 目标: 寻找类似 pss_test.m 中的重复结构
% 参数：积分窗口 W=128, 延迟 D=128

D = 219; % 延迟长度 (Delay)

if length(x) <= D
    error('信号长度不足以进行 D=%d 的延迟处理', D);
end

% 构造延迟流: r(n+D) 与 r(n)
rx_delayed = x(1+D:end);
rx_base    = x(1:end-D);

% 1. 共轭乘积: x*(n) * x(n+D)
conj_prod = conj(rx_base) .* rx_delayed;

% 2. 滑动求和 (使用 filter 实现)
% 使用 ones(1, W) 矩形窗进行移动求和
% 注意：x(1:end-D) 的长度为 L-D
P_metric = filter(ones(1, W), 1, conj_prod);

% 3. 能量计算
% 计算窗口内的能量用于归一化
rx_power = abs(rx_base).^2;
R_energy = filter(ones(1, W), 1, rx_power);

% 4. 计算最终度量 M(n)
% 修改说明：
% 根据 PSS 结构，Block0 是反相的 (-B)，Block1 是 B。
% 1. 在 Block0 和 Block1 交界处，rx_base 对应 Block0(-B)，rx_delayed 对应 Block1(B)。
%    此时 conj_prod = (-B)* . B = -|B|^2，求和也是负实数。
%    Real(P_metric) 应该出现显著负峰值。
% 2. 在 Block1 和 Block2 交界处，rx_base 对应 Block1(B)，rx_delayed 对应 Block2(B)。
%    Real(P_metric) 应该出现显著正峰值。
% 因此我们保留 P_metric 的实部（归一化后），这样可以看到正负极性。

% 仅归一化 P_metric (保留复数信息，稍后取实部)
M_complex = P_metric ./ (R_energy + 1e-10);

% 取实部作为特征曲线 (理论上 PSS 区域相位为 0 或 pi，即全是实数)
M_real = real(M_complex);

% 补齐长度以便与原信号对齐绘图 (补齐 D)
M_real_full = [M_real; zeros(D, 1)];

%% 4. 绘图展示
figure('Position', [100, 100, 1200, 600], 'Name', 'Sliding Autocorrelation Structure Search');

t_axis = (0:length(x)-1) + startSample0;

% 子图 1: 信号幅度
subplot(2,1,1);
plot(t_axis, abs(x), 'Color', [0.6 0.6 0.6]);
ylabel('幅度 |x[n]|');
title(['信号时域波形 - ' inFile ' (N=' num2str(length(x)) ')']);
grid on; axis tight;

% 子图 2: 滑动自相关实部 Real[M(n)]
subplot(2,1,2);
plot(t_axis, M_real_full, 'LineWidth', 1.5, 'Color', 'b');
hold on;
yline(0, 'k-', 'Zero'); % 零线
yline(0.5, 'r:', 'Pos Thresh'); 
yline(-0.5, 'r:', 'Neg Thresh');

% 在图中标注反相特性
title(sprintf('滑动自相关实部 (W=%d, D=%d) - 负峰值指示反相块起始', W, D));
ylabel('归一化相关值 (Real)');
xlabel('样本索引 (Samples)');
grid on; axis tight;
ylim([-1.1, 1.1]); % 范围改为 [-1, 1]

linkaxes(findall(gcf,'type','axes'), 'x');
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
