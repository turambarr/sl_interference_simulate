% estimate_sss_cp_60mhz.m
% 先降采样到60MHz，再在60MHz基带寻找精准 CP 长度与起点
clear; clc; close all;

inFile = 'sigtest1.iq';

% ===== 与 sss_demodulation.m 对齐的读数参数 =====
% 稍微往前多读一点，防止重采样边缘效应吃掉开头的 CP
read_start_sample = 16386+874*6-200;
read_length = 10000;

% 与 sss_demodulation.m 对齐的同步参数
fs_source = 409.6e6;
fs_target = 60e6;
sro_ppm = 0;
cfo_hz = 17556;

fprintf('Loading file: %s from %d, len %d...\n', inFile, read_start_sample, read_length);
[x_raw, ~] = iq_read_int16_le(inFile, read_start_sample, read_length);
x_raw = double(x_raw);
x_raw = x_raw - mean(x_raw);
x_raw = x_raw / (mean(abs(x_raw)) + 1e-12);

% --- 新增：频谱下变频 (DDC 归基带) ---
freq_shift_hz = 63e6; % 频谱归基带向左搬移 63MHz
fprintf('Shifting spectrum left by %.1f MHz at 409.6MHz...\n', freq_shift_hz/1e6);
t_vec = (0:length(x_raw)-1) / fs_source;
x_raw = x_raw .* exp(-1j * 2 * pi * freq_shift_hz * t_vec).';

% --- 新增：CFO 补偿（与 SSS 解调流程一致） ---
fprintf('Applying separately CFO Correction: %.2f Hz...\n', cfo_hz);
x_raw = x_raw .* exp(-1j * 2 * pi * cfo_hz * t_vec).';

% --- 新增：使用零相移(Zero-Phase)低通滤波器抗混叠 ---
% 与 SSS 解调流程保持一致：35MHz 截止
fprintf('Applying Zero-Phase Anti-aliasing LPF...\n');
Wn = 35e6 / (fs_source / 2);
% 设计一个适中阶数 (例如 30 阶) 的 FIR 滤波器
lpf_order = 30;
b_lpf = fir1(lpf_order, Wn);
% filtfilt 执行正反双向滤波：彻底消除群时延(相位偏移)，完美保护 CP 在时域上的绝对位置！
x_raw_filtered = filtfilt(b_lpf, 1, x_raw);

% 1. 重采样到 60MHz (采用 Farrow 分数阶延迟插值消除 resample 带来的边界失真)
fs_eff = fs_source * (1 + sro_ppm/1e6);

T_in = 1 / fs_eff;
T_out = 1 / fs_target;
% 避免两端由于缺少足量插值基准点而越界，舍弃最后微小的尾部
t_out = 0 : T_out : (length(x_raw_filtered)-3)*T_in; 

% 映射到输入序列的虚拟索引 (1-based)
idx_frac = t_out / T_in + 1; 
idx_base = floor(idx_frac);
mu = idx_frac - idx_base;

% 确保 base 索引不会越界 (Farrow Cubic 需要 base-1 到 base+2 共4个点)
valid_mask = (idx_base >= 2) & (idx_base <= length(x_raw_filtered)-2);
idx_base = idx_base(valid_mask);
mu = mu(valid_mask);

% Farrow 立方插值滤波器系数 (Cubic Lagrange)
h0 = -(mu - 1) .* (mu - 2) .* mu / 6;
h1 =  (mu - 1) .* (mu + 1) .* (mu - 2) / 2;
h2 = -(mu + 1) .* mu .* (mu - 2) / 2;
h3 =  (mu + 1) .* (mu - 1) .* mu / 6;

% 卷积求和（无边界失真效应），此处的源数据已替换为抗混叠滤波后的结果
x_60 = h0 .* x_raw_filtered(idx_base - 1) + ...
       h1 .* x_raw_filtered(idx_base) + ...
       h2 .* x_raw_filtered(idx_base + 1) + ...
       h3 .* x_raw_filtered(idx_base + 2);
x_60 = x_60(:); % 转为列向量

N = length(x_60);
% 60MHz下的有效数据体长度是确切的 1024 点
delay_D = 1024; 

% 2. 逐点共轭乘积 + 双边能量（Schmidl 风格归一化）
x_base = x_60(1 : N - delay_D);
x_del  = x_60(delay_D + 1 : end);
corr_prod = x_base .* conj(x_del);
pow_base = abs(x_base).^2;
pow_del  = abs(x_del).^2;

% 3. 遍历确切的基带 CP 长度 (通常在 10~150 之间，例如常见的 72, 144 等)
W_test = 10 : 1 : 150;
peak_vals = zeros(size(W_test));
for i = 1:length(W_test)
    W = W_test(i);
    P = filter(ones(W, 1), 1, corr_prod);
    E1 = filter(ones(W, 1), 1, pow_base);
    E2 = filter(ones(W, 1), 1, pow_del);

    % 归一化度量
    M = abs(P) ./ (sqrt(E1 .* E2) + 1e-10);

    % 关键修复：避免“低能量分母导致 M 虚高”，只在有效能量区统计
    e_ref = median(E1(E1 > 0));
    if isempty(e_ref) || ~isfinite(e_ref)
        e_ref = 1e-12;
    end
    e_th = max(1e-12, 0.2 * e_ref);
    valid_e = (E1 > e_th) & (E2 > e_th);

    M_valid = M(valid_e);
    if isempty(M_valid)
        peak_vals(i) = 0;
    else
        % 用高分位数替代 max，减少单点毛刺把整条曲线抹平
        peak_vals(i) = prctile(M_valid, 99);
    end
end

% 4. 用“分段线性最小误差”找拐点（比 Kneedle 更稳）
x = W_test(:); y = peak_vals(:);

% 先轻微平滑，减少局部抖动
y_fit = smoothdata(y, 'movmean', 5);

minSeg = 8; % 每段最少点数
bestErr = inf;
bestIdx = minSeg + 1;
for k = (minSeg+1):(length(x)-minSeg)
    x1 = x(1:k);   y1 = y_fit(1:k);
    x2 = x(k:end); y2 = y_fit(k:end);

    p1 = polyfit(x1, y1, 1);
    p2 = polyfit(x2, y2, 1);

    e1 = y1 - polyval(p1, x1);
    e2 = y2 - polyval(p2, x2);
    sse = sum(e1.^2) + sum(e2.^2);

    if sse < bestErr
        bestErr = sse;
        bestIdx = k;
    end
end

cp_idx = bestIdx;
exact_cp_60 = x(cp_idx);

% 5. 利用归一化相关度量寻找更稳定的起点
P_opt = filter(ones(exact_cp_60, 1), 1, corr_prod);
E1_opt = filter(ones(exact_cp_60, 1), 1, pow_base);
E2_opt = filter(ones(exact_cp_60, 1), 1, pow_del);
M_opt = abs(P_opt) ./ (sqrt(E1_opt .* E2_opt) + 1e-10);

% 同样做能量门限，避免低能量区误导峰值位置
e_ref_opt = median(E1_opt(E1_opt > 0));
if isempty(e_ref_opt) || ~isfinite(e_ref_opt)
    e_ref_opt = 1e-12;
end
e_th_opt = max(1e-12, 0.2 * e_ref_opt);
M_opt((E1_opt <= e_th_opt) | (E2_opt <= e_th_opt)) = 0;

% 局部平滑，去除毛刺
y_smooth = smoothdata(M_opt, 'movmean', 5); 
[~, max_idx] = max(y_smooth);
cp_start_offset_60 = max_idx - exact_cp_60 + 1;

% 6. 核心：计算可以直接填进 sss_demodulation.m 的 startSample
% SSS 有效数据的起点(60MHz下) = CP的起点 + CP的长度
sss_data_start_60 = cp_start_offset_60 + exact_cp_60;
% 将偏移量按比例换算回 409.6MHz 的原生采样点位置
offset_4096 = round((sss_data_start_60 - 1) * (409.6 / 60));
global_startSample = read_start_sample + offset_4096;

fprintf('\n=== 60MHz 基带 CP 分析完成 ===\n');
fprintf('60MHz 下精确 CP 长度 : %d 个点\n', exact_cp_60);
fprintf('CP 相对截取区(60MHz)起点偏移 : %d 个点\n', cp_start_offset_60);
fprintf('\n>>> 请将以下值填入 sss_demodulation.m 的读数参数:\n');
fprintf('read_start_sample = %d;\n', global_startSample);
fprintf('==============================\n');

figure;
plot(x, y, 'b-', 'LineWidth', 2); hold on; grid on;
plot([x(1), x(end)], [y(1), y(end)], 'k--');
plot(exact_cp_60, y(cp_idx), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
xline(exact_cp_60, 'r--'); yline(y(cp_idx), 'r--');
title('60MHz基带下相关得分曲线及分段线性拐点');