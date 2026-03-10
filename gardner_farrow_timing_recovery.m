% gardner_farrow_timing_recovery.m
% 修改的降采样与同步方案：
% 1. 先进行 3倍 整数降采样 (resample)
% 2. 使用 Farrow 插值器 + Gardner 误差检测器 (TED) 进行闭环定时恢复
% 目标：从非整数倍采样率中恢复出最佳符号时刻 (Strobe)

clear; clc; close all;

%% 1. 参数设置
inFile = 'sigtest8.iq'; 

% 自动获取文件全部长度进行处理
fInfo = dir(inFile);
fileBytes = fInfo.bytes;
totalSamples = floor(fileBytes / 4); % int16 I/Q = 4 bytes

% 定义此时我们关心的“信号范围” (Original ROI)
% 这是用户之前 manually set 的范围
roi_start_in = 15564;
roi_len_in   = 874*8;

startSample = 0;             % 从头读
readLen     = totalSamples;  % 读全部

% 基础采样率参数
fs_source   = 409.6e6;       % 源采样率
fs_symbol   = 60e6;          % 目标符号速率 (Symbol Rate)
freq_shift_hz = 63e6;        % 频谱归基带向左搬移 63MHz

% 整数降采样参数
int_decim   = 3;             % 第一步降采样倍数
fs_intermediate = fs_source / int_decim; % 中间采样率 (~136.53 MHz)

% 环路参数
Bn_T        = 0.01;          % 环路带宽归一化参数 (可以调小以获得更稳的锁定，但收敛慢)
zeta        = 0.707;         % 阻尼系数
sps_out     = 1;             % 最终输出每符号点数 (我们通常只需要1个最佳点)
sps_loop    = 2;             % 环路工作在 2倍符号率 (Gardner 需要中间点)

%% 2. 读取原始数据
fprintf('读取文件: %s (Start=%d, Len=%d)\n', inFile, startSample, readLen);
[x_raw, ~] = iq_read_int16_le(inFile, startSample, readLen);
x_raw = double(x_raw);
x_raw = x_raw - mean(x_raw); 
x_raw = x_raw / mean(abs(x_raw)); % 归一化

%% 2.x 频谱下变频 (DDC 归基带)
fprintf('Shifting spectrum left by %.1f MHz at 409.6MHz...\n', freq_shift_hz/1e6);
t_vec = (0:length(x_raw)-1) / fs_source;
x_raw = x_raw .* exp(-1j * 2 * pi * freq_shift_hz * t_vec).';

%% 2.1 时域波形概览 (新增)
figure('Position', [50, 50, 1000, 400], 'Name', 'Time Domain Overview');
plot(abs(x_raw), 'Color', [0.6 0.6 0.6]); hold on;
% 标记 ROI
xline(roi_start_in, 'g-', 'LineWidth', 1.5);
xline(roi_start_in + roi_len_in, 'r-', 'LineWidth', 1.5);
title(['原始信号幅度图 (绿色=Start, 红色=End) - ' inFile]);
xlabel('Sample Index (409.6 MHz)'); ylabel('Normalized Amplitude');
xlim([1 length(x_raw)]);
grid on;

%% 3. 重采样至 540MHz (9 SPS)
fprintf('Step 3: Resample to 540MHz (SPS=9)...\n');

% 目标采样率
fs_target = 540e6; 
% 计算重采样比例 P/Q
[P, Q] = rat(fs_target / fs_source);

fprintf('   Resampling ratio: %d/%d (%.4f)\n', P, Q, P/Q);

% 执行重采样 (MATLAB resample 自带抗混叠滤波)
x_high = resample(x_raw, P, Q);

% --- SPS 设置 ---
sps     = 9; % 540MHz / 60MHz = 9

% --- RRC 滤波器 (移除：信号特性不匹配导致的恶化) ---
% 用户反馈加了 RRC 后效果变差，说明该信号可能不是标准的脉冲成型信号，
%或者是恒包络信号 (如 GMSK/SDPSK)，RRC 会引入不必要的 ISI。
fprintf('   Skipping RRC Filtering (Reverted). Using raw resampled data (SPS=%d)...\n', sps);
x_filtered = x_high;

% beta_rrc = 0.5; 
% span     = 10; 
% h_rrc    = rcosdesign(beta_rrc, span, sps);
% x_filtered = conv(x_high, h_rrc, 'same');

%% 4. Farrow + Gardner 闭环定时恢复 (Input: 540MHz, Output: 2*SymRate)
fprintf('Step 4: Gardner Timing Recovery Loop (Input SPS=9, Loop SPS=2)...\n');

% --- 关键修正：鲁棒的归一化 (Robust Normalization) ---
% Gardner TED 的增益与输入幅度平方成正比。
% 为了防止 Burst 信号导致有效段幅度过大（从而导致环路震荡），
% 我们检查信号的峰值与均值之比，智能选择归一化因子。

% 1. 计算统计量
amp_vals = abs(x_filtered);
rms_val  = sqrt(mean(amp_vals.^2));
% 取 99% 分位点作为"近似峰值"，抗噪点干扰
sorted_amp = sort(amp_vals);
peak_ref   = sorted_amp(floor(length(sorted_amp) * 0.99));

% 2. 决策归一化因子
% 如果是连续信号，Peak 也就比 RMS 大一点点 (比如 QPSK 是 1.0 vs 1.0, 也就是 PAPR 低)
% 如果是 Burst，绝大部分是噪声(0)，RMS 会极低，导致 Peak/RMS 巨大。
if peak_ref > 3.0 * rms_val
    fprintf('   [Normalization] Detect Burst Signal (Peak >> RMS). Using Peak ref: %.4f\n', peak_ref);
    norm_factor = peak_ref; 
else
    fprintf('   [Normalization] Detect Continuous Signal. Using RMS ref: %.4f\n', rms_val);
    norm_factor = rms_val;
end

% 3. 执行归一化
x_filtered = x_filtered / norm_factor;

% --- 参数设置 ---
% 输入采样率: 540 MHz
% 目标符号率: 60 MHz
% 环路工作率: 120 MHz (SPS=2)
% 标称步进 (Nominal Step): 540 / 120 = 4.5
% 强制归一化输入信号功率
% 这一步非常关键！否则误差检测器的增益会乱套
pwr = mean(abs(x_filtered).^2); 
x_filtered = x_filtered / sqrt(pwr);
step_nominal = fs_target / (2 * fs_symbol); % 4.5

% 环路滤波器参数
Bn_T   = 0.01;   % 环路带宽
zeta   = 0.707;  % 阻尼系数
Kp_ted = 2.7;    % TED 增益
K0_nco = -1;     % NCO 增益

Kp = (4 * zeta * Bn_T) / (Kp_ted * K0_nco);
Ki = (4 * Bn_T^2)      / (Kp_ted * K0_nco);

fprintf('   Nominal Step: %.4f\n', step_nominal);
fprintf('   Loop Gains: Kp=%.6f, Ki=%.6f\n', Kp, Ki);

% --- 初始化状态 ---
cnt_nco     = 0;
idx_input   = 4; % 从第4个点开始，保证Farrow有足够历史 (i-2, i-1, i, i+1)
strobe_cnt  = 0;
sample_buf  = zeros(1, 3); % [y(n-2), y(n-1), y(n)] -> [PrevSym, Mid, CurrSym]
loop_state  = 0;
w_control   = 0; % 保持当前的控制量 (P+I)，避免奇偶抖动

err_history = []; % 记录误差
y_out_all   = []; % 记录所有输出 (SPS=2)

% 输入数据
x_in = x_filtered; 
len_in = length(x_in);

% 预分配内存 (估计输出长度)
est_out_len = ceil(len_in / step_nominal);
y_out_all = zeros(1, est_out_len);
out_idx = 0;

% --- 环路循环 ---
% 这里的可视化代码已移除，以保持清洁
while idx_input <= len_in - 2
    
    cnt_nco = cnt_nco - 1;
    
    if cnt_nco <= 0
        % 产生采样点 (Strobe)
        mu = cnt_nco + 1; % 小数延迟 (0 <= mu < 1)
        
        % Farrow 3rd Order Interpolation
        % 基准点: idx_input
        % 指针: p0=i-2, p1=i-1, p2=i, p3=i+1
        % 注意: idx_input 随着外层循环增加
        
        p0 = x_in(idx_input - 2);
        p1 = x_in(idx_input - 1);
        p2 = x_in(idx_input);
        p3 = x_in(idx_input + 1);
        
        c3 = -1/6 * p0 + 1/2 * p1 - 1/2 * p2 + 1/6 * p3;
        c2 =  1/2 * p0 - p1 + 1/2 * p2;
        c1 = -1/3 * p0 - 1/2 * p1 + p2 - 1/6 * p3;
        c0 = p1;
        
        % 计算插值输出
        y_now = ((c3 * mu + c2) * mu + c1) * mu + c0;
        
        % 保存输出
        out_idx = out_idx + 1;
        y_out_all(out_idx) = y_now;
        
        % [Visualization Removed]
        % debug_strobe_idx...
        
        % 更新 Gardner Buffer
        sample_buf = [sample_buf(2:3), y_now];
        strobe_cnt = strobe_cnt + 1;
        
        % 误差检测 (每2个采样点计算一次, SPS=2)
        % buffer: [Sym(k-1), Mid(k-0.5), Sym(k)]
        if mod(strobe_cnt, 2) == 0
            % Gardner Equation:
            % e = Re{ (y(k) - y(k-1)) * conj(y(k-0.5)) }
            % buffer(3) = Sym(k)
            % buffer(1) = Sym(k-1)
            % buffer(2) = Mid(k-0.5)
            
            ted_err = real( (sample_buf(3) - sample_buf(1)) * conj(sample_buf(2)) );
            err_history(end+1) = ted_err;
            
            % 环路滤波
            prop_term = ted_err * Kp;
            loop_state = loop_state + ted_err * Ki;
            
            % 更新控制量 P+I
            w_control = prop_term + loop_state;
        else
            % 奇数时刻保持控制量，避免抖动
        end
        
        % 更新 NCO
        % Step = Nominal + Control
        current_step = step_nominal + w_control;
        cnt_nco = cnt_nco + current_step;
        
    end
    
    idx_input = idx_input + 1;
end

% 截取实际输出
y_out_all = y_out_all(1:out_idx);
fprintf('   Timing Recovery Done. Output Samples: %d\n', length(y_out_all));


%% 4.1 可视化：[已移除]
% figure('Position', [100, 100, 1000, 400], 'Name', 'Gardner Sampling Instants');...

%% 5. 抽取最佳符号流 (Downsample to SPS=1)
% Gardner 稳定后，我们需要选取最佳采样点。
% SPS=2 的流：Sym, Mid, Sym, Mid...
% 也就是奇数或偶数索引。
% 为避开收敛过程，丢弃前 500 个符号
stable_start = 1000; % samples
if stable_start > length(y_out_all); stable_start = 1; end

% 这里的相位模糊度：究竟 1:2:end 是符号，还是 2:2:end 是符号？
% 可以看能量？通常 Mid 点的能量比 Sym 点小 (由于过零点多)。
pwr_odd  = mean(abs(y_out_all(stable_start:2:end)).^2);
pwr_even = mean(abs(y_out_all(stable_start+1:2:end)).^2);

if pwr_odd > pwr_even
    y_symbols_all = y_out_all(1:2:end);
    fprintf('   Selecting ODD samples as Symbols (Power: %.2f vs %.2f)\n', pwr_odd, pwr_even);
else
    y_symbols_all = y_out_all(2:2:end);
    fprintf('   Selecting EVEN samples as Symbols (Power: %.2f vs %.2f)\n', pwr_even, pwr_odd);
end

%% 6. ROI 截取与绘图
% 映射 ROI (Input Index -> Symbol Index)
% 1. Input(409.6M) -> 540M
% 我们其实可以不管中间过程，直接按点数比例映射
total_input_samples = length(x_raw);
total_output_symbols = length(y_symbols_all);

ratio_in_out = total_output_symbols / total_input_samples; % 应该接近 60/409.6
roi_start_sym = floor(roi_start_in * ratio_in_out);
roi_len_sym   = floor(roi_len_in * ratio_in_out);

% 修正起始点 (加上一些余量)
roi_start_sym = roi_start_sym + 10; % 稍微往后一点避开边界误差
if roi_start_sym < 1; roi_start_sym = 1; end

roi_end_sym = roi_start_sym + roi_len_sym;
if roi_end_sym > length(y_symbols_all); roi_end_sym = length(y_symbols_all); end

y_roi = y_symbols_all(roi_start_sym : roi_end_sym);
fprintf('   ROI Mapped: Input[%d] -> Symbol[%d] (Len=%d)\n', roi_start_in, roi_start_sym, roi_len_sym);

%% 7. 频偏估计与补偿 (M-th Power Spectrum Estimation)
% 针对 SDPSK / pi/2-BPSK (M=4)
fprintf('Step 7: Frequency Offset Compensation (M=4 Power)...\n');

M = 4;
% 使用 ROI 数据进行估计 (信号质量最好)
scan_data = y_roi; 

% 1. M次方 (去除调制信息)
z_M = scan_data .^ M;

% 2. FFT 频谱搜索
Nfft = 2^nextpow2(length(z_M)*4); % 补零提高频率分辨率
Z_spec = fftshift(fft(z_M, Nfft));
f_axis = (-Nfft/2 : Nfft/2 - 1) / Nfft; % 归一化频率 (-0.5 ~ 0.5)

% 3. 搜索峰值
[~, max_idx] = max(abs(Z_spec));
f_est_M = f_axis(max_idx); % M * Delta_f

% 4. 计算频偏
delta_f_norm = f_est_M / M; % cycle/symbol
fprintf('   Estimated CFO: %.6f (cycles/symbol)\n', delta_f_norm);

% 5. 补偿
% 对整个 ROI 进行补偿
idx_vec = 0:(length(y_roi)-1);
phase_correction = exp(-1j * 2 * pi * delta_f_norm * idx_vec);
y_roi_corr = y_roi .* phase_correction;

% --- 绘图 (更新版) ---
figure('Position', [200, 200, 1200, 500], 'Name', 'Gardner Results + CFO Comp');

% 1. 误差收敛曲线
subplot(1, 3, 1);
plot(err_history);
title('Timing Error (Unnormalized gain?)');
grid on; xlabel('Symbol Index'); ylabel('Error');
xlim([0 length(err_history)]);

% 2. 补偿后的星座图 (SPS=1)
subplot(1, 3, 2);
plot(real(y_roi_corr), imag(y_roi_corr), '.', 'Color', 'b');
hold on;
% 画一个单位圆做参考
theta = 0:0.01:2*pi;
plot(cos(theta), sin(theta), 'k:', 'LineWidth', 0.5);
axis square; grid on;
title(['ROI Constellation (Compensated)\nCFO=' num2str(delta_f_norm, '%.5f')]);
xlabel('I'); ylabel('Q');
xlim([-2 2]); ylim([-2 2]);

% 3. 补偿后的差分星座图
subplot(1, 3, 3);
% 差分: y[n] * conj(y[n-1])
y_diff_corr = y_roi_corr(2:end) .* conj(y_roi_corr(1:end-1));

plot(real(y_diff_corr), imag(y_diff_corr), '.', 'Color', [0 0.6 0]);
axis square; grid on;
title('ROI Differential Constellation (Final)');
xlabel('I'); ylabel('Q');
xlim([-2.5 2.5]); ylim([-2.5 2.5]);

%% 8. Costas Loop (Fine Phase Sync)
fprintf('Step 8: Costas Loop (Fine Phase Sync)...\n');

% 归一化幅度，确保环路增益稳定
% 这一步只是整个序列的均值归一化
y_costas_in = y_roi_corr / mean(abs(y_roi_corr)); 

% Loop Parameters (大幅降低增益，提高稳定性)
alpha = 0.005;  % 原 0.05 -> 0.005 (降低10倍)
beta  = 0.0001; % 原 0.0005 -> 0.0001 (降低5倍)
phase_out = 0;
freq_track  = 0;

y_costas = zeros(size(y_costas_in));
costas_err_hist = zeros(size(y_costas_in));

% 调试：记录频率跟踪轨迹
freq_hist = zeros(size(y_costas_in));

for k = 1:length(y_costas_in)
    % 1. Apply Loop Phase
    z = y_costas_in(k) * exp(-1j * phase_out);
    y_costas(k) = z;
    
    % 2. Phase Error Detector (Robust)
    % 改进：先进行幅度归一化，再计算误差，消除幅度噪声对鉴相器的影响
    z_norm = z / (abs(z) + 1e-10); % 防止除零
    
    % QPSK/SDPSK -> 4th Power PED
    % err = imag(z_norm^4) / 4
    % 范围限制在 [-0.25, 0.25] 之间
    err = 0.25 * imag(z_norm^4); 
    
    costas_err_hist(k) = err;
    
    % 3. Loop Filter
    freq_track  = freq_track + beta * err;
    phase_out   = phase_out + alpha * err + freq_track;
    
    freq_hist(k) = freq_track;
end

fprintf('   Costas Loop Locked. Final Freq Correction: %.6f\n', freq_track);

% --- 汇总与输出估计参数 ---
fprintf('\n================================================\n');
fprintf('   ESTIMATED SYNCHRONIZATION PARAMETERS\n');
fprintf('================================================\n');

% 1. SRO (时钟频偏)
% 理论步进: step_nominal
% 实际稳态步进: step_nominal + loop_state (积分项即为稳态误差)
% loop_state 是累积的，代表了为追踪 SRO 需要的额外步进
% Gardner 循环结束后 loop_state 保存着最后的状态
loop_integral_val = loop_state; 
avg_step_steady = step_nominal + loop_integral_val;

% 计算 ppm: (Actual - Nominal) / Nominal * 1e6
sro_ppm = (loop_integral_val / step_nominal) * 1e6;

fprintf('1. Sampling Rate Offset (SRO):\n');
fprintf('   Nominal Step    : %.6f samples\n', step_nominal);
fprintf('   Loop Integrator : %.6f samples (Steady State)\n', loop_integral_val);
fprintf('   Steady Step     : %.6f samples\n', avg_step_steady);
fprintf('   Clock Offset    : %.2f ppm\n', sro_ppm);

% 2. CFO (载波频偏)
% 总频偏 = 粗估 (M=4 FFT) + 精估 (Costas Loop Integrator)
% delta_f_norm 是粗估的归一化频偏 (cycles/symbol)
% freq_track 是 Costas 环的归一化频偏 (rad/symbol)，需转为 cycles/symbol
% Costas loop update: phase += ... + freq_track
% Compensation: exp(-j * phase)
% 这意味着 freq_track 正在追踪残留的频偏 (rad/s)

delta_f_fine_cycles = freq_track / (2*pi); % rad -> cycle
total_cfo_norm      = delta_f_norm + delta_f_fine_cycles; % cycles/symbol

% 转换为 Hz
% 1 cycle/symbol = fs_symbol Hz
cfo_hz = total_cfo_norm * fs_symbol;

fprintf('2. Carrier Frequency Offset (CFO):\n');
fprintf('   Coarse Est (M=4): %.6f (cycles/symbol)\n', delta_f_norm);
fprintf('   Fine Est (Loop) : %.6f (cycles/symbol)\n', delta_f_fine_cycles);
fprintf('   Total Norm Freq : %.6f (cycles/symbol)\n', total_cfo_norm);
fprintf('   Total CFO       : %.2f Hz (at SymRate %d Msps)\n', cfo_hz, fs_symbol/1e6);
fprintf('================================================\n');

% --- 最终绘图 (含 Costas) ---
figure('Position', [250, 250, 1200, 500], 'Name', 'Final Results (Gardner + CFO + Costas)');

% 1. Costas 误差收敛
subplot(1, 3, 1);
plot(costas_err_hist);
hold on;
plot(freq_hist * 100, 'r', 'LineWidth', 1.5); % 把频率轨迹也画上去（放大显示）
legend('Phase Err', 'Freq Track (x100)');
title('Costas Loop Status');
grid on; xlabel('Symbol Index'); 
ylim([-0.3 0.3]); % 限制显示范围，误差被限幅了

% 2. 最终星座图 (Costas Output)
subplot(1, 3, 2);
plot(real(y_costas), imag(y_costas), '.', 'Color', 'b');
hold on;
theta = 0:0.01:2*pi;
plot(cos(theta), sin(theta), 'k:', 'LineWidth', 0.5);
axis square; grid on;
title(['Final Constellation (Costas)\nResid Freq=' num2str(freq_track, '%.5f')]);
xlabel('I'); ylabel('Q');
xlim([-2 2]); ylim([-2 2]);

% 3. 最终差分星座图
subplot(1, 3, 3);
y_diff_final = y_costas(2:end) .* conj(y_costas(1:end-1));
plot(real(y_diff_final), imag(y_diff_final), '.', 'Color', [0 0.6 0]);
axis square; grid on;
title('Final Differential Constellation');
xlabel('I'); ylabel('Q');
xlim([-2.5 2.5]); ylim([-2.5 2.5]);


%% 9. [Final Step] Demodulation & m-sequence Extraction
fprintf('\nStep 9: SDPSK Demodulation and m-sequence Extraction...\n');

% 1. 获取同步后的符号流
% 假设 y_costas 是已经锁定相位的符号 (1 Sample/Symbol)
% 截取前 1024 个点 (8 blocks * 128 symbols)
total_sym_needed = 8 * 128;
if length(y_costas) < total_sym_needed
    warning('Not enough symbols tracked! Expected 1024, got %d', length(y_costas));
    y_sync = y_costas;
else
    y_sync = y_costas(1:total_sym_needed);
end

% 2. 差分检波 (Differential Demodulation)
% SDPSK 信息包含在相邻符号的相位差中
% d[k] = y[k] * conj(y[k-1])
% 对于 pi/2-BPSK, 差分相位应该是 +pi/2 (Bit 0) 或 -pi/2 (Bit 1)
% 即 d[k] 应该接近 +j 或 -j
% 注意：如果是 m 序列结构，且块之间有反相 (++-+-+-+)，
% 差分检波会自动消除块反相的影响 (因为 (-A)*conj(-B) = A*conj(B))
% 所以 differential demod 后的比特流应该是 8 个完全相同的 m 序列 (除了边界点)

% 为了保持长度一致，我们假定第一个符号通过某种方式参考，或者直接丢弃第一个差分
% 这里为了对齐，我们直接计算 full differential
y_diff_vec = y_sync(2:end) .* conj(y_sync(1:end-1));

% 判决: Imaginary part > 0 -> 0, < 0 -> 1
% (这取决于具体的映射规则，暂定 Im>0 为 0, Im<0 为 1)
raw_soft_bits = -imag(y_diff_vec); % 正值为1，负值为0 (逻辑取反适应习惯)

% 3. 规整化与合并 (Reshape & Combine)
% --- [关键修改] 使用差分检波数据 ---
% 优势：自动消除 Block 的正负反相影响 ((-A)*(-B)' = AB')
% 我们直接使用 raw_soft_bits (Step 2 计算出的)
% 注意：raw_soft_bits 长度为 1023 (因为丢了第一个)，需要补齐或截断
% 为了整齐，我们丢弃每一块的第一个差分点？
% 或者假设 raw_soft_bits 对应的是 1..1023
% 没关系，先 reshape 看看
target_len = 128 * 8;
if length(raw_soft_bits) < target_len
    raw_soft_bits(end+1 : target_len) = 0; % 补零
else
    raw_soft_bits = raw_soft_bits(1:target_len);
end

% 不再使用相干解调结果
% soft_syms = real(y_bpsk); 

% 4. 结构合并
% 此时不需要再根据 ++-+-+-+ 手动翻转了，差分天然去除了符号翻转
try
    mat_syms = reshape(raw_soft_bits, 128, 8);
    
    % [已移除] 手动翻转代码
    % mat_syms(:, 3) = -mat_syms(:, 3); ...
    
    % --- 验证：逐块输出 ---
    % 将每一块独立解调并打印，方便肉眼比对一致性
    fprintf('\n--- Individual Block Verification (Sign Corrected) ---\n');
    for col = 1:8
         col_soft = mat_syms(:, col);
         col_bits = col_soft > 0;
         col_hex = '';
         for i = 1:4:128
            chunk = col_bits(i:i+3);
            val = chunk(1)*8 + chunk(2)*4 + chunk(3)*2 + chunk(4)*1;
            col_hex = [col_hex, dec2hex(val)];
         end
         
         % 简单的质量评估 (软比特的平均幅度，越大越确信)
         avg_mag = mean(abs(col_soft));
         fprintf('Block %d (Conf=%.2f): %s\n', col, avg_mag, col_hex);
    end
    fprintf('------------------------------------------------\n');

    % 合并 (求平均) - 提升 SNR
    m_seq_soft = mean(mat_syms, 2);
    
    % 5. 判决与 Hex 输出
    % 判决: >0 -> ? <0 -> ?
    % 这一步有 180度 模糊度 (即全 0 变 全 1)。只能输出一种，另一种是其反码。
    % 假设 >0 为 1
    m_bits = m_seq_soft > 0;
    
    % 转 Hex
    % 128 bits = 128/4 = 32 Hex digits
    hex_str = '';
    for i = 1:4:128
        chunk = m_bits(i:i+3);
        % Convert [b3 b2 b1 b0] to int
        val = chunk(1)*8 + chunk(2)*4 + chunk(3)*2 + chunk(4)*1;
        hex_str = [hex_str, dec2hex(val)];
    end
    
    fprintf('\n------------------------------------------------\n');
    fprintf('Extracted m-sequence (Hex):\n');
    fprintf('%s\n', hex_str);
    fprintf('(Note: If incorrect, try inverting all bits: ~Hex)\n');
    fprintf('------------------------------------------------\n');
    
catch e
    fprintf('Error in structural decoding: %s\n', e.message);
end

return;
