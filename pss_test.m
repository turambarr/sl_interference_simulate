clc; clear; close all;

%% 1. 参数配置
N = 1024;               % 符号长度 (FFT Size)
%% oversample_K = 4;    % (已移除：不再通过固定倍数推算)
Fs_high = 409.6e6;      % 高采样率 (Hz)
Fs_base = 60e6;         % 基带采样率 (带宽设置为 60MHz)
freq_offset = 1000;     % 频偏 (Hz)

L_block = 128;          % PSS 分块长度 (N/8)
SNR_dB = 10;            % 信噪比
num_data_symbols = 4;   % 模拟的数据符号数量 (Data Payload)

%% 2. 信号生成

% --- (1) 生成 PSS (修改为 SDPSK) ---
% 原始 Hex 序列
Qpss_hex = 'C1B5D191024D3DC3F8EC52FAA16F3958';

% Hex 转 Bin
Qpss_bin = [];
for i = 1 : length(Qpss_hex)
    hex_val = hex2dec(Qpss_hex(i));
    bin_str = dec2bin(hex_val, 4);
    Qpss_bin = [Qpss_bin, bin_str - '0']; % 转换为 0/1 数组
end

% SDPSK 调制逻辑
% 映射: 0/1 -> -1/1 (或者反转)
% 差分相位: dphi = b * (pi/2)
b = flip(2 * Qpss_bin - 1); 

dphi = b * (pi/2);
phase_seq = cumsum(dphi);
pss_block_base = exp(1j * phase_seq);

% 构建时域 PSS:
% 结构: [B, B, -B, B, -B, B, -B, B]
pattern_signs = [1, 1, -1, 1, -1, 1, -1, 1];
pss_signal = [];
for s = pattern_signs
    pss_signal = [pss_signal, s * pss_block_base];
end

% PSS Cyclic Suffix (-B 的前32个点)
suffix_len = 32;
pss_signal = [pss_signal, -pss_block_base(1:suffix_len)];

% --- (2) 生成 SSS (OFDM + 4QAM Fixed) ---
fprintf('生成 SSS 信号 (OFDM + 4QAM Fixed Pattern)...\n');
% 构造固定频域数据: 0, 1, 2, 3, 0, 1 ... 便于解调验证
sss_data_len = N - 200; % 824 个有效子载波
sss_fixed_msg = mod(0:sss_data_len-1, 4); % 生成固定序列

% 4QAM 映射 (Phase: pi/4, 3pi/4, -3pi/4, -pi/4)
sss_qam_syms = exp(1j * (sss_fixed_msg * pi/2 + pi/4));

% OFDM 调制 (映射到子载波)
sss_spectrum = zeros(1, N);
half_sc = floor(sss_data_len / 2);
% 正频率部分 (Index 2 : half+1)
sss_spectrum(2 : half_sc+1) = sss_qam_syms(1:half_sc);
% 负频率部分 (Index N-half+1 : N)
sss_spectrum(N-half_sc+1 : N) = sss_qam_syms(half_sc+1:end);

% IFFT 变换 (无 CP)
sss_signal = ifft(sss_spectrum, N) * sqrt(N);

% 保存参考数据用于接收端验证
sss_ref_data = sss_fixed_msg;
sss_ref_syms = sss_qam_syms;

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

tx_signal_base = [noise_pre, pss_signal, sss_signal, data_payload, noise_post];

%% 4. 信道模拟 (上采样 + 噪声 -> 接收前端)
fprintf('正在进行上采样和信道处理...\n');
% 定义基带时间轴 (用于参考)
time_base = (0:length(tx_signal_base)-1) / Fs_base;

% A. 上采样 (Resample: 60MHz -> 409.6MHz)
[P, Q] = rat(Fs_high / Fs_base); 
tx_signal_high = resample(tx_signal_base, P, Q);

% 重新计算时间轴 (因为 resample 后长度变了)
time_high = (0:length(tx_signal_high)-1) / Fs_high;

% B. 信道 (添加频偏)
% 原代码: rx_signal_high = tx_signal_high; (无频偏)
% 新代码: 添加 3 MHz 频偏
cfo_val = 3e6; % 3 MHz
t_vec = (0:length(tx_signal_high)-1) / Fs_high;
phase_rot = exp(1j * 2 * pi * cfo_val * t_vec);
rx_signal_high = tx_signal_high .* phase_rot; 

% D. 添加信道噪声 (AWGN) 到高信号上
sig_power = mean(abs(rx_signal_high).^2);
noise_power = sig_power / (10^(SNR_dB/10));
noise_scale = sqrt(noise_power/2);
noise_awgn = (randn(size(rx_signal_high)) + 1j*randn(size(rx_signal_high))) * noise_scale;

rx_signal_high = rx_signal_high + noise_awgn;

% 为了兼容 第5部分 滑动自相关代码 (原本它是拿基带信号算的)
% 我们暂时粗暴降采样一下给 第5、6 部分用来画图/找峰值
% 注意：rx_signal 现在是带噪声的基带信号
rx_signal = resample(rx_signal_high, Q, P); % 降回来

fprintf('信号处理完成。\n');

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

%% 7. 高精度定时恢复 (Gardner + Farrow)
% 此部分替代单纯的降采样，使用闭环算法寻找最佳采样点

% --- A. 准备数据: 截取 PSS 附近的高速采样数据 ---
% 这里的 rx_signal_high 是 409.6 MHz 的信号 (含噪声+频偏)
% 我们根据 Section 5 (滑动自相关) 得到的粗略位置来截取
% 粗略位置 (60MHz): noise_len + 1 (因为是仿真已知)
% 实际工程中应使用 argmax(M_n) 来确定 start_idx_est

% 使用已知位置演示:
known_start_time = time_base(noise_len + 1);
% 为了让 Gardner 有足够的锁定时间 (虽然 PSS 很短), 我们往前多取一点
start_offset_time = 20 * (1/Fs_base); % 往前 20 个符号
t_extract_start = known_start_time - start_offset_time;
if t_extract_start < 0; t_extract_start = 0; end

% 往后取足够长 (PSS + SSS + Some Data = 300 symbols)
extract_duration = 300 * (1/Fs_base); 
t_extract_end   = t_extract_start + extract_duration;

% 转为 High Rate 索引
idx_start_high = find(time_high >= t_extract_start, 1);
idx_end_high   = find(time_high <= t_extract_end, 1, 'last');
if isempty(idx_end_high); idx_end_high = length(rx_signal_high); end

rx_segment_high = rx_signal_high(idx_start_high : idx_end_high);

% --- B. RRC 匹配滤波 (在 409.6MHz 进行) ---
fprintf('Step 7.B: RRC 匹配滤波 (409.6MHz)...\n');
rolloff = 0.5; span = 10; sps_raw = Fs_high / Fs_base;
sps_high_res = 100;
h_proto = rcosdesign(rolloff, span, sps_high_res, 'sqrt');
t_proto = linspace(-span/2, span/2, length(h_proto));
t_query = -span/2 : (1/sps_raw) : span/2; 
h_rrc   = interp1(t_proto, h_proto, t_query, 'spline');
h_rrc   = h_rrc / sqrt(sum(h_rrc.^2));

rx_filtered_high = conv(rx_segment_high, h_rrc, 'same');

% --- C. 粗降采样至 ~2 SPS ---
% Gardner 需要约 2 SPS 的输入。
% 我们从 409.6降到 120 (SPS=2)。
target_sps = 2;
fs_intermediate = target_sps * Fs_base; % 120 MHz
[P2, Q2] = rat(fs_intermediate / Fs_high);
rx_mid = resample(rx_filtered_high, P2, Q2);

% --- AGC / 幅度能一化 (关键修正) ---
% Gardner 误差检测器的增益直接依赖于信号幅度。
% 必须将信号归一化到单位幅度附近，否则计算出的 Kp/Ki 将失效。
rx_mid = rx_mid - mean(rx_mid);       % 去直流
rx_mid = rx_mid / mean(abs(rx_mid));  % 归一化幅度

% --- D. Farrow + Gardner Loop ---
fprintf('Step 7.D: 执行 Gardner 闭环定时恢复...\n');
% 参数
fs_sym_loop = Fs_base; % 60 MHz
nominal_step = fs_intermediate / (2 * fs_sym_loop); % ~1
Bn_T = 0.02; % 带宽稍大一点以快速锁定
zeta = 0.707;
Kp_ted = 2.7; K0_nco = -1;
Kp = (4 * zeta * Bn_T) / (Kp_ted * K0_nco);
Ki = (4 * Bn_T^2)      / (Kp_ted * K0_nco);

cnt_nco = 0; idx_input = 3; 
sample_buffer = zeros(1, 3);
strobe_cnt = 0;
loop_filter_state = 0;
y_out_all = [];
timing_errs = [];

len_mid = length(rx_mid);
while idx_input < len_mid - 2
    cnt_nco = cnt_nco - 1;
    if cnt_nco <= 0
        mu = cnt_nco + 1; % 
        % Farrow Interpolation
        p0 = rx_mid(idx_input - 2); p1 = rx_mid(idx_input - 1);
        p2 = rx_mid(idx_input);     p3 = rx_mid(idx_input + 1);
        c3 = -1/6 * p0 + 1/2 * p1 - 1/2 * p2 + 1/6 * p3;
        c2 =  1/2 * p0 - p1 + 1/2 * p2;
        c1 = -1/3 * p0 - 1/2 * p1 + p2 - 1/6 * p3;
        c0 = p1;
        y_now = ((c3 * mu + c2) * mu + c1) * mu + c0;
        
        y_out_all = [y_out_all, y_now];
        sample_buffer = [sample_buffer(2:3), y_now];
        strobe_cnt = strobe_cnt + 1;
        
        if mod(strobe_cnt, 2) == 0
            % Gardner Error: (y(n) - y(n-1)) * conj(y(n-1/2))
            % buffer: [Sym(n-1), Mid, Sym(n)]
            err = real( (sample_buffer(3) - sample_buffer(1)) * conj(sample_buffer(2)) );
            timing_errs = [timing_errs, err];
            
            % Loop Filter
            pll_out = err * Kp + loop_filter_state;
            loop_filter_state = loop_filter_state + err * Ki;
            
            cnt_nco = cnt_nco + (nominal_step + pll_out);
        else
             cnt_nco = cnt_nco + (nominal_step + loop_filter_state);
        end
    end
    idx_input = idx_input + 1;
end

% 抽取最佳采样点 (SPS=1)
y_symbols = y_out_all(2:2:end);

%% 8. Costas 环 (去除残留频偏与相位)
% Gardner 只管定时，不管相位/频偏。
fprintf('Step 8: Costas 环载波同步...\n');
x_loop_in = y_symbols;
N_sym = length(x_loop_in);
x_synced = zeros(size(x_loop_in));

phase_est = -pi/4; % 初始相位猜测
freq_est = 0;
% Costas 参数
alpha_c = 0.05; beta_c = 0.002;

for n = 1:N_sym
    val = x_loop_in(n);
    z = val * exp(-1j * phase_est);
    x_synced(n) = z;
    
    % QPSK/SDPSK 鉴相
    I = real(z); Q = imag(z);
    err_c = Q * sign(I) - I * sign(Q);
    
    freq_est = freq_est + beta_c * err_c;
    phase_est = phase_est + freq_est + alpha_c * err_c;
end

% 旋转 -45 度以便观察
x_final = x_synced * exp(-1j * pi/4);

%% 9. 结果展示
% figure('Position', [200, 200, 1000, 400]);
% subplot(1,3,1); plot(timing_errs); title('Gardner Timing Error'); grid on;
% subplot(1,3,2); plot(real(x_loop_in), imag(x_loop_in), '.'); title('Gardner Output (No Costas)'); axis square;
% subplot(1,3,3); plot(real(x_final), imag(x_final), 'b.'); 
% hold on; 
% % 绘制差分
% x_diff = x_final(2:end) .* conj(x_final(1:end-1));
% plot(real(x_diff), imag(x_diff), 'g.'); 
% title('Final Output (Blue) & Diff (Green)'); 
% legend('Sync', 'Diff'); axis square; grid on; xlim([-2 2]); ylim([-2 2]);

%% 9.1 差分眼图 (SPS=2)
% 使用 Gardner 输出的 2倍过采样信号查看眼图
% figure('Position', [250, 250, 800, 600], 'Name', 'Differential Eye Diagram (SPS=2)');

% % 计算 SPS=2 的差分
% % d[k] = y[k] * conj(y[k-2]) (因为 2 samples = 1 symbol)
% y_all_diff = y_out_all(3:end) .* conj(y_out_all(1:end-2));

% % 折叠显示 (2 Symbols width = 4 samples)
% % 由于 SPS Only=2, 眼图可能略显折线感
% trace_width = 4; % 2 Symbols
% n_traces = floor(length(y_all_diff) / 2); % Step size = 2 (1 Symbol)
% eye_I = zeros(n_traces, trace_width);
% eye_Q = zeros(n_traces, trace_width);

% for k = 1:n_traces
%     % 滑动窗口
%     idx_start = (k-1)*2 + 1;
%     idx_end   = idx_start + trace_width - 1;
%     if idx_end > length(y_all_diff); break; end
    
%     seg = y_all_diff(idx_start : idx_end);
%     eye_I(k, :) = real(seg);
%     eye_Q(k, :) = imag(seg);
% end

% subplot(2,1,1);
% plot(eye_I', 'b', 'Color', [0 0 1 0.05]); 
% title('Diff Eye Diagram (I/Real) - SPS=2'); grid on;
% subplot(2,1,2);
% plot(eye_Q', 'r', 'Color', [1 0 0 0.05]);
% title('Diff Eye Diagram (Q/Imag) - SPS=2'); grid on;

% 验证误码
% dphi_rx = angle(x_diff);
% b_est = sign(dphi_rx);

% 我们只验证 PSS 部分
% 需要找到 x_final 中 PSS 的起始点
% 我们是基于 known_start_time 往前 20个 符号截取的
% 加上 loop 延迟，大概在 20 附近
% offset_est = 20; 
% valid_len = 127;
% if length(b_est) >= offset_est + valid_len
%     b_est_pss = b_est(offset_est : offset_est + valid_len - 1);
%     b_ref_segment = b(2:valid_len+1);
    
%     errors = sum(b_est_pss ~= b_ref_segment);
%     fprintf('PSS 解调误码: %d / %d\n', errors, valid_len);
% end

%% 10. SSS 解调验证 (Ideal Case)
fprintf('\nStep 10: SSS 解调与验证 (Ideal Timing)...\n');
% 从基带接收信号中直接截取 SSS 用于验证
% 理想起始点: noise_len + length(pss_signal) + 1
rss_idx_start = noise_len + length(pss_signal) + 1;
rss_idx_end   = rss_idx_start + N - 1;

if rss_idx_end <= length(rx_signal)
    rss_rx = rx_signal(rss_idx_start : rss_idx_end);
    
    % 1. FFT
    rss_rx_freq = fft(rss_rx, N) / sqrt(N);
    
    % 2. 提取有效子载波
    % 正频率
    rx_sc_pos = rss_rx_freq(2 : half_sc+1);
    % 负频率
    rx_sc_neg = rss_rx_freq(N-half_sc+1 : N);
    
    rx_sc_total = [rx_sc_pos, rx_sc_neg];
    
    % 3. 绘制星座图
    figure('Position', [300, 300, 500, 500], 'Name', 'SSS Constellation (Ideal Timing)');
    plot(real(rx_sc_total), imag(rx_sc_total), 'b.'); hold on;
    % 绘制参考中心点
    ref_pts = [1+1j, -1+1j, -1-1j, 1-1j] / sqrt(2); % QPSK
    plot(real(ref_pts), imag(ref_pts), 'rx', 'MarkerSize', 10, 'LineWidth', 2);
    title('SSS Constellation (After FFT)'); axis square; grid on;
    xlabel('I'); ylabel('Q'); xlim([-2 2]); ylim([-2 2]);
    
    % 4. 简单的硬判决与误码率
    rx_phase = angle(rx_sc_total);
    % 映射回 0,1,2,3: (phase - pi/4) / (pi/2) -> rounds to int
    rx_int = mod(round((rx_phase - pi/4) / (pi/2)), 4);
    
    bit_errs = sum(rx_int ~= sss_ref_data);
    fprintf('   SSS Symbol Errors (Raw): %d / %d\n', bit_errs, length(sss_ref_data));
else
    fprintf('   [Warning] Received signal too short to verify SSS.\n');
end
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