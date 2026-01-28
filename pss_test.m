clc; clear; close all;

%% 1. 参数配置
N = 1024;               % 符号长度 (FFT Size)
oversample_K = 4;       % 过采样因子 (Oversampling Factor)
L_block = 128;          % PSS 分块长度 (N/8)
SNR_dB = 15;            % 信噪比
num_data_symbols = 4;   % 模拟的数据符号数量 (Data Payload)

%% 2. 信号生成

% --- (1) 生成 PSS (基于您提供的逻辑) ---
% 原始 Hex 序列
Qpss_hex = 'C1B5D191024D3DC3F8EC52FAA16F3958';

% Hex 转 Bin (手动实现以确保兼容性)
Qpss_bin = [];
for i = 1 : length(Qpss_hex)
    hex_val = hex2dec(Qpss_hex(i));
    bin_str = dec2bin(hex_val, 4);
    Qpss_bin = [Qpss_bin, bin_str - '0']; % 转换为 0/1 数组
end

% 双极性映射 + 翻转 (b = flip(2*x - 1))
b = flip(2 * Qpss_bin - 1);

% 生成基础块 (Differential Encoding)
pss_block_base = zeros(1, L_block);
for k = 1 : L_block
    b_sum = sum(b(1:k));
    pss_block_base(k) = exp(1j * pi * (-0.25 - 0.5 * b_sum));
end

% 构建时域 PSS:
% 结构: [B, B, -B, B, -B, B, -B, B]
% 1. (B, B)   -> +1
% 2. (B, -B)  -> -1 (负峰起始)
% 3. (-B, B)  -> -1
% 4. (B, -B)  -> -1
% 5. (-B, B)  -> -1
% ...
pattern_signs = [1, 1, -1, 1, -1, 1, -1, 1];
pss_signal = [];
for s = pattern_signs
    pss_signal = [pss_signal, s * pss_block_base];
end

% 新增需求: 拼一个 -B 的前32个点在序列最后面 (Cyclic Suffix)
suffix_len = 32;
pss_signal = [pss_signal, -pss_block_base(1:suffix_len)];

% --- (2) 生成 SSS (简化为随机 OFDM 符号) ---
sss_signal = gen_ofdm_symbol(N);

% --- (3) 生成 Data Payload (OFDM 数据符号) ---
data_payload = [];
for i = 1 : num_data_symbols
    data_payload = [data_payload, gen_ofdm_symbol(N)];
end

%% 3. 组装帧结构
% 结构: [前置噪声] + [PSS] + [SSS] + [Data x 4] + [后置噪声]
noise_len = 600;
noise_pre = (randn(1, noise_len) + 1j*randn(1, noise_len)) * 0.1;
noise_post = (randn(1, noise_len) + 1j*randn(1, noise_len)) * 0.1;

tx_signal = [noise_pre, pss_signal, sss_signal, data_payload, noise_post];

% (User requested removal of oversampling/downsampling logic)


%% 4. 添加信道噪声 (AWGN)
sig_power = mean(abs(tx_signal).^2);
noise_power = sig_power / (10^(SNR_dB/10));
noise_scale = sqrt(noise_power/2);
noise_awgn = (randn(size(tx_signal)) + 1j*randn(size(tx_signal))) * noise_scale;

rx_signal = tx_signal + noise_awgn;

%% 5. 接收机算法实现

% === 算法 A: 滑动自相关 (Sliding Autocorrelation) ===
% 原理: P(d) = sum( r(m)' * r(m+L) )
% 窗口 W = L_block, 延迟 D = L_block

W = L_block; % 积分窗口长度 (128 * K)
D = L_block; % 延迟长度 (128 * K)

% 构造延迟流: r(n+D)
% 有效长度: length(rx) - D
rx_delayed = rx_signal(1+D:end);
rx_base    = rx_signal(1:end-D);

% 共轭相乘
conj_prod = conj(rx_base) .* rx_delayed;

% 滑动求和 (使用 filter 实现)
% 'ones(1,W)' 是矩形窗，filter 做移动求和
P_metric = filter(ones(1, W), 1, conj_prod);

% 能量归一化 (Energy Normalization)
rx_power = abs(rx_base).^2;
R_energy = filter(ones(1, W), 1, rx_power);

% 计算最终度量 M(n) - 修改为保留实部以观察反相特性
% M_n = abs(P_metric) ./ (R_energy + 1e-10); % 原代码
M_complex = P_metric ./ (R_energy + 1e-10);
M_n = real(M_complex);

% 调整 M_n 的长度以便绘图 (补齐尾部延迟)
M_n = [M_n, zeros(1, D)]; 

% === 算法 B: 互相关 (Cross Correlation / Matched Filter) ===
% 使用本地已知 PSS 进行相关
[acor, lag] = xcorr(rx_signal, pss_signal);
% 取模并归一化
cross_corr = abs(acor);
cross_corr = cross_corr / max(cross_corr);

% 截取有效部分 (只看正延迟部分，对齐时间轴)
% xcorr 的中心点是 lag=0，对应 PSS 完全对齐
center_idx = find(lag == 0);
% 简单的对齐映射用于绘图
start_plot_idx = center_idx - noise_len; 
% 为了画图方便，我们直接截取和 rx_signal 长度一致的部分
% (实际工程中需要根据 lag 换算)
valid_len = length(rx_signal);
cross_corr_plot = cross_corr(center_idx : min(end, center_idx + valid_len - 1));
% 补齐长度
if length(cross_corr_plot) < valid_len
    cross_corr_plot = [cross_corr_plot, zeros(1, valid_len - length(cross_corr_plot))];
end

%% 6. 绘图展示
figure('Position', [100, 100, 1000, 800]);

% 子图 1: 帧结构 (幅度)
subplot(3,1,1);
plot(abs(rx_signal), 'Color', [0.5 0.5 0.5]); hold on;
xline(noise_len, 'r--', 'LineWidth', 1.5); text(noise_len, max(abs(rx_signal)), ' PSS Start', 'Color', 'r');
xline(noise_len+N, 'g--', 'LineWidth', 1.5); text(noise_len+N, max(abs(rx_signal)), ' SSS Start', 'Color', 'g');
xline(noise_len+2*N, 'b--', 'LineWidth', 1.5); text(noise_len+2*N, max(abs(rx_signal)), ' Data Start', 'Color', 'b');
title(['接收信号时域波形 (包含 ' num2str(num_data_symbols) ' 个 OFDM 数据符号)']);
ylabel('幅度'); xlim([1 length(rx_signal)]);
legend('Rx Envelope', 'Frame Boundaries');

% 子图 2: 滑动自相关 (Schmidl & Cox) - 实部
subplot(3,1,2);
plot(M_n, 'LineWidth', 1.5, 'Color', 'b'); hold on;
yline(0, 'k-', 'Zero'); % 零线
yline(-0.5, 'r:', 'Neg Peak'); % 负峰参考
% 高亮 PSS 区域 (因为Y轴现在是[-1,1]，需要重新获取yl或手动指定)
ylim([-1.1 1.1]);
yl = ylim;
patch([noise_len noise_len+N noise_len+N noise_len], [yl(1) yl(1) yl(2) yl(2)], 'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
title('滑动自相关 M(n)实部 - (注意：负峰起始 + 正值平台)');
ylabel('归一化度量(Real)'); xlim([1 length(rx_signal)]);
grid on;

% 子图 3: 互相关 (Matched Filter)
subplot(3,1,3);
plot(cross_corr_plot, 'LineWidth', 1.5, 'Color', [0.5 0 0.5]);
title('互相关 (Matched Filter) - 精确检测峰值');
ylabel('归一化相关值'); xlim([1 length(rx_signal)]);
grid on;
xline(noise_len, 'r--', 'True Start');

%% 7. PSS 星座图分析
% 提取接收信号中的 PSS 全长 (1024点)
pss_rx_full = rx_signal(noise_len + 1 : noise_len + length(pss_signal));

% 提取 PSS 的第一个 Block (长度 128*K)
% 也就是所谓的 "取 PSS 的 1/8" (因为 1024/128 = 8)
L_unit = L_block;
pss_rx_block1 = pss_rx_full(1 : L_unit);

figure('Position', [1120, 100, 600, 800], 'Name', 'PSS Constellation');

% 上图：完整 PSS (1024*K 点)
subplot(2,1,1);
plot(real(pss_rx_full), imag(pss_rx_full), 'b.', 'MarkerSize', 6);
hold on; xline(0,'k'); yline(0,'k'); grid on; axis equal;
title(sprintf('完整 PSS (%d 点) 星座图 (SNR=%d dB)', length(pss_rx_full), SNR_dB));
xlabel('I'); ylabel('Q');

% 下图：单截段 PSS (128*K 点)
subplot(2,1,2);
plot(real(pss_rx_block1), imag(pss_rx_block1), 'r.', 'MarkerSize', 8);
hold on; xline(0,'k'); yline(0,'k'); grid on; axis equal;
title(sprintf('PSS 第一个块 (%d 点) 星座图', length(pss_rx_block1)));
xlabel('I'); ylabel('Q');

% 简单的 EVM 估计
evm_val = std(abs(pss_rx_full) - mean(abs(pss_rx_full)));
text(min(xlim)*0.8, min(ylim)*0.8, sprintf('Full Amp Std: %.4f', evm_val));

%% 辅助函数: 生成随机 OFDM 符号
function sym = gen_ofdm_symbol(N)
    % 1. 参数
    num_subcarriers = N - 200; % 留出保护带 (Guard Band)
    
    % 2. 生成随机 QPSK 数据
    bits = randi([0 3], 1, num_subcarriers);
    qpsk_data = exp(1j * (bits * pi/2 + pi/4));
    
    % 3. 子载波映射 (Map to Subcarriers)
    % 频率结构: [0(DC), 正频率..., 0...0(Guard), 负频率...]
    spectrum = zeros(1, N);
    half = floor(num_subcarriers / 2);
    
    % 正频率部分 (Index 2 : half+1) -> MATLAB下标从1开始，DC是1
    spectrum(2 : half+1) = qpsk_data(1:half);
    
    % 负频率部分 (Index N-half+1 : N)
    spectrum(N-half+1 : N) = qpsk_data(half+1:end);
    
    % 4. IFFT
    sym = ifft(spectrum, N) * sqrt(N); % 能量归一化
end