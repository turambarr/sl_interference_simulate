% gardner_sps_demo.m
% 使用 Gardner 环路进行符号定时同步 (Timing Recovery)
% 目的：自适应跟踪采样率误差 (SROS)，消除时钟漂移，输出最佳采样点的星座图

clear; clc; close all;

%% 1. 参数设置
filename = 'sigtest8.iq'; 
% 用户指定的 PSS 起始偏移
start_offset = 19934 - 874 * 5;   
read_len = 874 * 12;      % 读取长度

% 系统速率参数
Fs_input = 409.6e6;       % 输入采样率
Rs_symbol = 60e6;         % 符号速率
SPS_nominal = Fs_input / Rs_symbol; % 标称每符号采样数 (~6.8266)

%% 2. 读取数据
if ~isfile(filename)
    error('文件 %s 不存在', filename);
end
[x_raw, ~] = iq_read_int16_le(filename, start_offset, read_len);
x_raw = double(x_raw);
x_raw = x_raw - mean(x_raw);
x_raw = x_raw / mean(abs(x_raw)); % 归一化

%% 3. Gardner 环路初始化
% 目标：输出 2 SPS (Gardner 算法需要每符号 2 个点：中间点和最佳点)
target_osr = 2; 

% 初始步进 (输入样点/输出样点)
% 我们希望输出间隔是 0.5 个符号，即 0.5 * SPS_nominal
step_nominal = SPS_nominal / target_osr; 
step_curr = step_nominal;

% 环路滤波器参数
% 阻尼系数 zeta, 带宽 BL
zeta = 0.707;
BL_normalized = 0.005; % 归一化环路带宽 (相对于符号速率)

% 计算 PI 控制器增益
Kp = 2 * zeta * BL_normalized;
Ki = BL_normalized^2;

% 状态变量
idx_curr = 20; % 从第 20 个点开始，留出插值余量
mu = 0;        % 分数延迟 (0 <= mu < 1)
strobe = 0;    % 抽样时刻指示

rx_syms_buffer = zeros(1, ceil(read_len / step_nominal)+100);
rx_idx = 0;

% 误差记录
err_log = [];
step_log = [];

% 环路积分器
NCO_acc = 0;

fprintf('开始 Gardner 环路同步 (Nominal Step = %.4f)...\n', step_nominal);

% 预分配三个寄存器用于 Gardner TED 
% r0: y[k] (当前), r1: y[k-1] (半符号前/Mid), r2: y[k-2] (一符号前/Prev)
% 注意：实际上 Gardner 是在 2 SPS 下工作
% 偶数点是符号点，奇数点是中间点，或者反之。
% 我们在每个输出点都计算，利用 strobe 降采样。
% 但标准实现通常是：一直是 2 SPS 输出，只在“符号时刻”更新误差。

% 简化实现：Interpolated Timing Recovery Loop
% 参考: F. M. Gardner, "A BPSK/QPSK Timing-Error Detector for DSP Receivers"

out_trace = [];
timing_err = 0;

while idx_curr < length(x_raw) - 10
    
    % --- 1. 插值 (Cubic Lagrange Farrow 结构简化 或 Linear) ---
    % 线性插值：y = (1-mu)*x[i] + mu*x[i+1]
    % i = floor(idx_curr)
    
    int_idx = floor(idx_curr);
    frac = idx_curr - int_idx; % 当前的真实 mu 是基于 idx 的
    % 注意：通常 mu 是由 NCO 控制的下一次采样时刻的小数部分
    % 这里我们用 idx_curr 直接表示浮点索引
    
    % 立方插值 (4-point) 精度更好
    p_vec = x_raw(int_idx-1 : int_idx+2);
    % 利用 Farrow 结构计算插值 (Cubic)
    % y = sum(c_k * mu^k)
    % 这里使用 MATLAB 简单的样条或 interp1 并不适合逐点循环
    % 手写 Lagrange 3rd order 或 Catmull-Rom
    % 简化：线性插值 (对于高过采样率通常足够，这里 ~6.8 SPS 足够)
    sample_now = (1-frac) * x_raw(int_idx) + frac * x_raw(int_idx+1);
    
    % 存储输出
    rx_idx = rx_idx + 1;
    rx_syms_buffer(rx_idx) = sample_now;
    
    % --- 2. 误差检测 (Gardner TED) ---
    % 需要 3 个点: y[k-2] (Prev Sym), y[k-1] (Mid), y[k] (Curr Sym)
    % 只要 rx_idx 是偶数，我们就假设完成了一个符号周期 (2 samples per symbol)
    
    if rx_idx >= 3 && mod(rx_idx, 2) == 0
        y_curr = rx_syms_buffer(rx_idx);     % y(n)
        y_mid  = rx_syms_buffer(rx_idx-1);   % y(n - 1/2)
        y_prev = rx_syms_buffer(rx_idx-2);   % y(n - 1)
        
        % Gardner Error:
        % e = real( (y_prev - y_curr) * conj(y_mid) )
        % 适用于 QPSK / BPSK 等
        
        ted_err = real( (y_prev - y_curr) * conj(y_mid) );
        
        % --- 3. 环路滤波 (Loop Filter) ---
        NCO_acc = NCO_acc + Ki * ted_err;
        step_adj = Kp * ted_err + NCO_acc;
        
        % 更新步进 (负反馈: 误差大则减慢或加快)
        step_curr = step_nominal - step_adj;
        
        err_log(end+1) = ted_err;
        step_log(end+1) = step_curr;
    else
        % 非符号时刻，保持步进
        % step_curr 保持不变
    end
    
    % --- 4. 更新索引 ---
    idx_curr = idx_curr + step_curr;
    
end

% 截取有效输出
rx_syms = rx_syms_buffer(1:rx_idx);

% 降采样取符号点 (1, 3, 5...) 或 (2, 4, 6...) 取决于锁定相位
% 通常取偶数点作为最佳采样点 (因为那是做 TED 计算的点)
syms_out = rx_syms(2:2:end); 

% 计算最终收敛的 SPS
final_step = mean(step_log(end-100:end));
estimated_sps = final_step * target_osr;

fprintf('Gardner 环路完成。\n');
fprintf('收敛 SPS: %.5f (Nominal: %.5f)\n', estimated_sps, SPS_nominal);

%% 4. 后处理：Costas 环 (去频偏)
% Gardner 只管定时，不管频偏。定时对准后，仍需 Costas 去掉相位旋转。

fprintf('启动 Costas 环去频偏...\n');
x_loop_in = syms_out;
N_sym = length(x_loop_in);
x_synced = zeros(size(x_loop_in));

BL_T_costas = 0.05; 
w_n_c = BL_T_costas / (0.707 + 1/(4*0.707));
alpha_c = 2 * 0.707 * w_n_c;
beta_c  = w_n_c^2;

phase_est = -pi/4; % 初始相位
freq_est = 0;

for n = 1:N_sym
    z_in = x_loop_in(n);
    z_out = z_in * exp(-1j * phase_est);
    x_synced(n) = z_out;
    
    I = real(z_out); Q = imag(z_out);
    % QPSK PED
    err = Q * sign(I) - I * sign(Q);
    
    freq_est = freq_est + beta_c * err;
    phase_est = phase_est + freq_est + alpha_c * err;
end

% 旋转回坐标轴显示
x_final = x_synced * exp(-1j * pi/4);

%% 5. 绘图
figure('Position', [100, 100, 1000, 600], 'Name', 'Gardner Timing Recovery');

subplot(2, 2, 1);
plot(step_log * target_osr);
title('SPS 跟踪曲线 (Estimated SPS)');
ylabel('Samples Per Symbol'); xlabel('Symbol Index');
yline(SPS_nominal, 'r--');
grid on;

subplot(2, 2, 2);
plot(err_log);
title('Gardner TED 误差');
grid on;

subplot(2, 2, 3);
plot(real(x_final), imag(x_final), 'k.', 'MarkerSize', 6);
axis square; grid on; xline(0); yline(0);
title('同步后星座图 (Timing + Carrier)');
xlabel('I'); ylabel('Q');
xlim([-2 2]); ylim([-2 2]);

subplot(2, 2, 4);
% 绘制前 200 个点的眼图 (实部)
% 需要重塑为眼图形式，这里简单画一下过采样序列
output_trace_real = real(rx_syms(1:min(end, 400)));
plot(output_trace_real, '-o', 'MarkerSize', 2);
title('输出波形片段 (2 sps)');
grid on;
