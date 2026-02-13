% constellation_interp.m
% 使用插值法 (interp1) 进行非整数倍降采样并绘制星座图
% 目的：替代 resample 函数，更精确地控制采样时刻（相位），解决采样点对不准的问题

clear; clc; close all;

%% 1. 参数设置
inFile = 'sigtest8.iq'; 
startSample = 19934-874*5;   % 读取起始点
readLen     = 874*8;           % 读取长度 (适当读长一点，保证足够插值)

D = 6.3975;                 % 降采样倍率 (每个符号占用的原始样点数)
interpolation_method = 'spline'; % 插值方法: 'linear', 'spline', 'pchip'

fs = 409.6e6;          % 采样率

%% 2. 读取原始数据
fprintf('读取文件: %s (Start=%d, Len=%d)\n', inFile, startSample, readLen);
[x_raw, ~] = iq_read_int16_le(inFile, startSample, readLen);
x_raw = double(x_raw);
x_raw = x_raw - mean(x_raw); % 去直流

% 归一化幅度
x_raw = x_raw / mean(abs(x_raw));

%% 2.1 频偏预处理 (已禁用，改用后续 Costas 环)
fprintf('跳过预先频偏补偿，将在插值后使用 Costas 环进行同步...\n');
f_offset = 0; % 设为 0，让 Costas 环去追踪

%% 3. 插值核心逻辑
% 原始时间轴: 0, 1, 2, ..., N-1
t_raw = 0:(length(x_raw)-1);

% 我们不再搜索最佳的采样相位，而是使用固定的 Offset
% 如果需要微调，请修改这里的 offset_ratio (0 ~ 1 之间)
offset_ratio = 0; 
current_offset = offset_ratio * D;

fprintf('使用固定 Offset: %.4f (Ratio=%.2f)\n', current_offset, offset_ratio);

%% 3.5 抗混叠滤波器 (Anti-Aliasing Filter)
% 在降采样前，对 x_raw (近似 Baseband 但采样率高) 进行低通滤波
% Fs=409.6MHz, Target=60MHz => Nyquist=30MHz
fc_cutoff = 32e6; 
fprintf('设计并应用抗混叠滤波器 (Fc=%.1f MHz)...\n', fc_cutoff/1e6);

[b_lpf, a_lpf] = butter(5, fc_cutoff / (fs/2), 'low');

% 分别对 I/Q 路进行滤波，防止某些环境下的 filtfilt 不支持复数
fprintf('分别针对 I 路和 Q 路进行零相位滤波...\n');
x_I = real(x_raw);
x_Q = imag(x_raw);

x_I_filt = filtfilt(b_lpf, a_lpf, x_I);
x_Q_filt = filtfilt(b_lpf, a_lpf, x_Q);

x_filtered = complex(x_I_filt, x_Q_filt);

% 调试信息：检查滤波前后 Q 路能量是否异常
pwr_I_pre = mean(x_I.^2); pwr_Q_pre = mean(x_Q.^2);
pwr_I_post = mean(x_I_filt.^2); pwr_Q_post = mean(x_Q_filt.^2);
fprintf('   Pre-Filter: P_I=%.4f, P_Q=%.4f, Ratio(Q/I)=%.4f\n', pwr_I_pre, pwr_Q_pre, pwr_Q_pre/pwr_I_pre);
fprintf('   Post-Filter: P_I=%.4f, P_Q=%.4f, Ratio(Q/I)=%.4f\n', pwr_I_post, pwr_Q_post, pwr_Q_post/pwr_I_post);

% 【关键】滤波后能量会衰减，必须重新归一化，否则星座图会变很小
x_filtered = x_filtered / mean(abs(x_filtered));

%% 4. 执行单次插值与绘图
figure('Position', [100, 100, 800, 800], 'Name', 'Constellation (Fixed Offset)');

% --- 构造新的采样时间点 ---
% 从 current_offset 开始，每次步进 D，直到不超过原始数据长度
t_new = current_offset : D : (length(x_raw)-1);

% --- 执行插值 (interp1) ---
% 使用滤波后的信号 x_filtered 进行插值
x_resampled = interp1(t_raw, x_filtered, t_new, interpolation_method);

% --- 绘图 (星座图) ---
plot(x_resampled, '.', 'Color', 'b', 'MarkerSize', 8);
axis square; grid on;
title(sprintf('简单插值降采样星座图 (Fixed Offset)\nD=%.4f', D));
xlabel('In-Phase (I)'); ylabel('Quadrature (Q)');
xlim([-2 2]); ylim([-2 2]);

fprintf('\n=== 结果说明 ===\n');
fprintf('点数: %d\n', length(x_resampled));

%% 6. 绘制差分星座图 (Differential Constellation)
% 针对 DPSK/SDPSK 信号，信息包含在相邻符号的相位差中
% z_diff[k] = z[k] * conj(z[k-1])
% 这可以消除绝对相位模糊，甚至抵消部分低频相位噪声

fprintf('正在生成差分星座图...\n');
x_diff = x_resampled(2:end) .* conj(x_resampled(1:end-1));

figure('Position', [500, 100, 500, 500], 'Name', 'Differential Constellation');
plot(real(x_diff), imag(x_diff), '.', 'Color', [0 .6 0], 'MarkerSize', 8);
axis square; grid on;
title('差分星座图 (Differential)');
xlabel('Real( z[k] * z*[k-1] )'); 
ylabel('Imag( z[k] * z*[k-1] )');
xlim([-2.5 2.5]); ylim([-2.5 2.5]);
xline(0, 'k:'); yline(0, 'k:');

%% 7. SDPSK 解调与 Hex 输出
fprintf('\n正在尝试 SDPSK 解调...\n');

% 1. 频偏校正 (基于差分信号)
% 理想 SDPSK 差分相位为 +pi/2 或 -pi/2
% 平方后应全为 -1 (相位 pi)
% 偏差即为 2*CFO_phase
% 利用 mean(x_diff.^2) 估计平均相位偏差
mean_sq = mean(x_diff.^2);
phase_bias = angle(mean_sq) / 2;
% 如果 phase_bias 接近 pi/2 或 -pi/2，可能是被 flip 了，这里假设 CFO 较小
% 将星座图旋转回 +/- pi/2 附近 (即虚轴)
x_diff_corr = x_diff * exp(-1j * (phase_bias - pi)); % 这里的 -pi 是因为理想平方是 -1 (angle=pi)，我们要把它转到 angle=pi 的位置?
% 修正逻辑:
% 理想 x_diff 就在 +/- 1j. x_diff^2 就在 -1.
% 实际 x_diff^2 的角度是 phase_sq.
% 误差是 phase_sq - pi.
% 因此我们需要把 x_diff 旋转 -(phase_sq - pi)/2.
rot_angle = -(angle(mean_sq) - pi) / 2;
x_diff_corr = x_diff * exp(1j * rot_angle);

% 2. 判决 (Soft Decision: Angle / Hard Decision: Sign of Imag)
% 映射: angle > 0 -> bit 1 (+pi/2), angle < 0 -> bit 0 (-pi/2)
% 注意: 之前的生成逻辑是 b = 2*bit - 1. b=1 -> +pi/2. b=-1 -> -pi/2.
% 所以 +pi/2 对应 bit 1, -pi/2 对应 bit 0.
demod_bits_all = zeros(1, length(x_diff_corr));
demod_bits_all(imag(x_diff_corr) > 0) = 1;

% 3. 寻找并提取 Hex
% 我们知道 Block 长度为 128
% 原始生成逻辑是 b = flip(2*Qpss_bin - 1).
% 这意味着时间上先发出的数据对应 Qpss_bin 的末尾。
% 若要还原 Qpss_hex，我们需要把解调出的比特流 flip 回去。

% 尝试滑动窗口寻找最佳匹配? 或者直接按 128 分块打印
% 由于读了 readLen = 874*8，大概有 1000+ 个符号
% PSS 结构是 [B, B, -B, B, -B, B, -B, B]
% 我们尝试打印每 128 个符号对应的 Hex
L_blk = 127; % m-sequence length often 127, but here pss_test used 128
% pss_test.m: L_block = 128. Qpss_hex has 32 chars -> 128 bits.
% So we expect 128 bits chunks.

fprintf('--- 解调数据流 (每 128 bit 一行) ---\n');
num_chunks = floor(length(demod_bits_all) / 128);

for k = 1:num_chunks
    idx_start = (k-1)*128 + 1;
    bits_chunk = demod_bits_all(idx_start : idx_start+127);
    
    % 尝试1: 直接转换 (Assume Forward Order)
    hex_str_fwd = bits2hex_str(bits_chunk);
    
    % 尝试2: 翻转后转换 (Assume Reverse Order / Flip)
    bits_chunk_flipped = flip(bits_chunk);
    hex_str_rev = bits2hex_str(bits_chunk_flipped);
    
    % 尝试3: 取反 (Invert bits) - 对应 -B 块
    hex_str_inv = bits2hex_str(~bits_chunk);
    hex_str_rev_inv = bits2hex_str(~bits_chunk_flipped);
    
    fprintf('Chunk %d: Raw=%s | Flip=%s\n', k, hex_str_fwd, hex_str_rev);
end

% 辅助函数：Bit流转Hex字符串
function h_str = bits2hex_str(b_vec)
    % b_vec: 1xN array of 0/1
    % Pad to multiple of 4
    len = length(b_vec);
    rem4 = mod(len, 4);
    if rem4 > 0
        b_vec = [b_vec, zeros(1, 4-rem4)];
    end
    
    h_str = '';
    for i = 1:4:length(b_vec)
        nibble = b_vec(i:i+3);
        val = nibble(1)*8 + nibble(2)*4 + nibble(3)*2 + nibble(4)*1;
        h_str = [h_str, dec2hex(val)];
    end
end




