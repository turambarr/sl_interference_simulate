% gardner_farrow_timing_recovery.m
% 彻底修改的降采样与同步方案：
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
roi_start_in = 19934 - 874*5;
roi_len_in   = 874*8;

startSample = 0;             % 从头读
readLen     = totalSamples;  % 读全部

% 基础采样率参数
fs_source   = 409.6e6;       % 源采样率
fs_symbol   = 60e6;          % 目标符号速率 (Symbol Rate)

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

% 跳过 RRC 滤波，直接使用重采样后的信号
fprintf('   Skipping RRC Filtering, using raw resampled data (SPS=%d)...\n', sps);
x_filtered = x_high;

%% 4. Farrow + Gardner 闭环定时恢复 (Input: 540MHz, Output: 2*SymRate)
fprintf('Step 4: Gardner Timing Recovery Loop (Input SPS=9, Loop SPS=2)...\n');

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
y_costas_in = y_roi_corr / mean(abs(y_roi_corr)); 

% Loop Parameters (调整增益)
alpha = 0.05; 
beta  = 0.0005;
phase_out = 0;
freq_track  = 0;
y_costas = zeros(size(y_costas_in));
costas_err_hist = zeros(size(y_costas_in));

for k = 1:length(y_costas_in)
    % Apply Loop Phase
    z = y_costas_in(k) * exp(-1j * phase_out);
    y_costas(k) = z;
    
    % Phase Error Detector (QPSK/SDPSK -> 4th Power)
    % error = imag(z^4) pushes to axes (0, 90, 180, 270)
    % Scaling factor 0.25 to reduce magnitude
    err = 0.25 * imag(z^4); 
    costas_err_hist(k) = err;
    
    % Loop Filter
    freq_track  = freq_track + beta * err;
    phase_out   = phase_out + alpha * err + freq_track;
end

% --- 最终绘图 (含 Costas) ---
figure('Position', [250, 250, 1200, 500], 'Name', 'Final Results (Gardner + CFO + Costas)');

% 1. Costas 误差收敛
subplot(1, 3, 1);
plot(costas_err_hist);
title('Costas Phase Error');
grid on; xlabel('Symbol Index'); ylabel('Error');
xlim([0 length(costas_err_hist)]);
ylim([-0.5 0.5]);

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

return;
