% 终极融合版：星链上行 PSS 物理层完美验收接收机
% 包含：Farrow+Gardner定时恢复 -> Costas相干解调 -> 物理相轨迹重构 -> 精准对齐验收
clear; clc; close all;
fprintf('======================================================\n');
fprintf('     👑 终极验收：基带环路恢复 + 物理模板完美对线 👑\n');
fprintf('======================================================\n');

%% 1. 参数设置与读取
inFile = 'sigtest9.iq'; 
roi_start_in = 15564 + 874*8 - 1000;
roi_len_in   = 874*8;

fInfo = dir(inFile);
totalSamples = floor(fInfo.bytes / 4);
[x_raw, ~] = iq_read_int16_le(inFile, 0, totalSamples);
x_raw = double(x_raw);
x_raw = x_raw - mean(x_raw); 
x_raw = x_raw / mean(abs(x_raw)); 

fs_source = 409.6e6; 
fs_symbol = 60e6; 
fs_target = 540e6; 

%% 2. 重采样至 540MHz
fprintf('Step 1: 重采样至 540MHz 与鲁棒归一化...\n');
[P, Q] = rat(fs_target / fs_source);
x_filtered = resample(x_raw, P, Q);

amp_vals = abs(x_filtered);
rms_val  = sqrt(mean(amp_vals.^2));
sorted_amp = sort(amp_vals);
peak_ref   = sorted_amp(floor(length(sorted_amp) * 0.99));
if peak_ref > 3.0 * rms_val
    norm_factor = peak_ref; 
else
    norm_factor = rms_val;
end
x_filtered = x_filtered / norm_factor;
x_filtered = x_filtered / sqrt(mean(abs(x_filtered).^2)); 

%% 3. Farrow + Gardner 闭环定时恢复
fprintf('Step 2: 运行 Gardner 定时恢复环路寻找最佳采样点 (Strobe)...\n');
step_nominal = fs_target / (2 * fs_symbol); 
Bn_T = 0.01; zeta = 0.707; Kp_ted = 2.7; K0_nco = -1;     
Kp = (4 * zeta * Bn_T) / (Kp_ted * K0_nco);
Ki = (4 * Bn_T^2)      / (Kp_ted * K0_nco);

len_in = length(x_filtered);
y_out_all = zeros(1, ceil(len_in / step_nominal));
cnt_nco = 0; idx_input = 4; strobe_cnt = 0; out_idx = 0;
sample_buf = zeros(1, 3); loop_state = 0; w_control = 0;

while idx_input <= len_in - 2
    cnt_nco = cnt_nco - 1;
    if cnt_nco <= 0
        mu = cnt_nco + 1;
        p0 = x_filtered(idx_input - 2); p1 = x_filtered(idx_input - 1);
        p2 = x_filtered(idx_input);     p3 = x_filtered(idx_input + 1);
        
        c3 = -1/6*p0 + 1/2*p1 - 1/2*p2 + 1/6*p3; c2 = 1/2*p0 - p1 + 1/2*p2;
        c1 = -1/3*p0 - 1/2*p1 + p2 - 1/6*p3;     c0 = p1;
        
        y_now = ((c3 * mu + c2) * mu + c1) * mu + c0;
        out_idx = out_idx + 1; y_out_all(out_idx) = y_now;
        sample_buf = [sample_buf(2:3), y_now];
        strobe_cnt = strobe_cnt + 1;
        
        if mod(strobe_cnt, 2) == 0
            ted_err = real( (sample_buf(3) - sample_buf(1)) * conj(sample_buf(2)) );
            w_control = (ted_err * Kp) + loop_state + (ted_err * Ki);
            loop_state = loop_state + (ted_err * Ki);
        end
        cnt_nco = cnt_nco + step_nominal + w_control;
    end
    idx_input = idx_input + 1;
end
y_out_all = y_out_all(1:out_idx);

%% 4. 抽取最佳符号流与极其宽容的 ROI 截取
stable_start = 1000; 
pwr_odd  = mean(abs(y_out_all(stable_start:2:end)).^2);
pwr_even = mean(abs(y_out_all(stable_start+1:2:end)).^2);
if pwr_odd > pwr_even
    y_symbols_all = y_out_all(1:2:end);
else
    y_symbols_all = y_out_all(2:2:end);
end

idx_540M = roi_start_in * (fs_target / fs_source);
idx_120M = idx_540M / step_nominal;
true_sym_start = floor((idx_120M - stable_start) / 2);

% 极其宽容的截取裕量，保证 1024 个符号完整落入池中
roi_start_sym = max(1, true_sym_start - 500); 
roi_end_sym = min(length(y_symbols_all), true_sym_start + 1024 + 2000);
y_roi = y_symbols_all(roi_start_sym : roi_end_sym);

%% 5. M=4 频偏估计与 Costas 锁相
fprintf('Step 3: M=4 频偏粗补与 Costas 载波跟踪...\n');
M = 4;
z_M = y_roi .^ M;
Nfft = 2^nextpow2(length(z_M)*4); 
Z_spec = fftshift(fft(z_M, Nfft));
f_axis = (-Nfft/2 : Nfft/2 - 1) / Nfft; 
[~, max_idx] = max(abs(Z_spec));
delta_f_norm = f_axis(max_idx) / M; 

idx_vec = 0:(length(y_roi)-1);
y_roi_corr = y_roi .* exp(-1j * 2 * pi * delta_f_norm * idx_vec);
y_costas_in = y_roi_corr / mean(abs(y_roi_corr)); 

alpha = 0.005; beta = 0.0001; phase_out = 0; freq_track = 0;
y_costas = zeros(size(y_costas_in));
for k = 1:length(y_costas_in)
    z = y_costas_in(k) * exp(-1j * phase_out);
    y_costas(k) = z;
    z_norm = z / (abs(z) + 1e-10); 
    err = 0.25 * imag(z_norm^4); 
    freq_track  = freq_track + beta * err;
    phase_out   = phase_out + alpha * err + freq_track;
end
fprintf('  >> 环路收敛！获得 1 SPS 锁相符号流 (y_costas)。\n');

%% ==================== 以下为验收对比核心模块 ==================== %%

%% 6. 构建本地纯净 1024 符号 (坚守 8 Block 全正极性)
fprintf('Step 4: 构建本地物理层理想模板 (8 块, 全正极性)...\n');
hexStr = '64B5E9E680A4438BA8FC2676C60DF195'; 
bits = zeros(1, 128);
for k = 1:length(hexStr)
    val = hex2dec(hexStr(k));
    bits((k-1)*4+1 : k*4) = bitget(val, 4:-1:1);
end

phase_diff = zeros(1, 128);
phase_diff(bits == 1) = -pi/2; 
phase_diff(bits == 0) =  pi/2; 
phases = cumsum(phase_diff);
local_syms_single = exp(1j * phases);

% 极性全为正
block_signs = [1, 1, 1, 1, 1, 1, 1, 1];
local_syms_1024 = kron(block_signs, local_syms_single);
local_syms_1024 = local_syms_1024(:); 
N_syms = 1024;

%% 7. 互相关精确找齐与【强制回退提取】
fprintf('Step 5: 互相关精确捕获与强制回退...\n');
[R_sym, lags_sym] = xcorr(y_costas(:), local_syms_1024);
[~, max_idx] = max(abs(R_sym));

% 核心破局点：xcorr 会避开第一块的瞬态畸变对齐到第二块。强制向左回退 128 个符号找回真实起点。
sym_start_idx = lags_sym(max_idx) + 1 - 128; 

fprintf('  >> xcorr 算法捕获点: %d\n', lags_sym(max_idx) + 1);
fprintf('  >> 强制回退 128 符号，真实 PSS 起点定格: %d\n', sym_start_idx);

if sym_start_idx < 1 || sym_start_idx + N_syms - 1 > length(y_costas)
    error('截取失败：向左回退超出了数组边界！');
end
x_rx_pss = y_costas(sym_start_idx : sym_start_idx + N_syms - 1);
x_rx_pss = x_rx_pss(:);

%% 8. 功率归一化与【四重相位模糊】自动盲对齐
fprintf('Step 6: 四象限相位模糊盲消除...\n');
x_rx_pss = x_rx_pss / sqrt(mean(abs(x_rx_pss).^2));
local_syms_1024 = local_syms_1024 / sqrt(mean(abs(local_syms_1024).^2));

base_phase = angle(mean(x_rx_pss .* conj(local_syms_1024)));
trial_phases = [-pi, -pi/2, 0, pi/2];
mse_results = zeros(size(trial_phases));

for iP = 1:length(trial_phases)
    temp_aligned = local_syms_1024 * exp(1j * (base_phase + trial_phases(iP)));
    mse_results(iP) = mean(abs(x_rx_pss - temp_aligned).^2);
end
[~, best_idx] = min(mse_results);
final_offset = base_phase + trial_phases(best_idx);
local_syms_aligned = local_syms_1024 * exp(1j * final_offset);

%% 9. 【终极杀手锏】消除残留频偏 (剔除瞬态区精准拟合)
fprintf('Step 7: 稳态区精准残留频偏拉平 (Detrending)...\n');
raw_phase_diff = unwrap(angle(x_rx_pss)) - unwrap(angle(local_syms_aligned));

% 核心修复：只用环路收敛后的稳定数据（索引 129~1024）来计算频偏斜率，防止被 Block 1 带偏
stable_idx = 129:N_syms; 
p_fit = polyfit(stable_idx, raw_phase_diff(stable_idx)', 1); 

% 将算出的真实稳态斜率外推应用到全部 1024 个点
fine_phase_ramp = polyval(p_fit, 1:N_syms)'; 

% 补偿残留频偏
x_rx_pss_perfect = x_rx_pss .* exp(-1j * fine_phase_ramp);
error_vec_perfect = x_rx_pss_perfect - local_syms_aligned;

% 计算全段 EVM (含畸变的 Block 1) 与 稳态 EVM (Block 2~8)
evm_pct_total = sqrt(mean(abs(error_vec_perfect).^2)) * 100;
error_vec_stable = error_vec_perfect(stable_idx);
evm_pct_stable = sqrt(mean(abs(error_vec_stable).^2)) * 100;

%% 10. 逐点精细对比与符号错误统计
fprintf('\n======================================================\n');
fprintf('     🔍 最终验收：1024 物理符号逐点解剖 🔍\n');
fprintf('======================================================\n');

phase_error_rad = angle(x_rx_pss_perfect .* conj(local_syms_aligned));
phase_error_deg = rad2deg(phase_error_rad);

% 以 45 度为硬判决门限
symbol_errors_total = abs(phase_error_rad) > (pi/4);
symbol_errors_stable = symbol_errors_total(stable_idx);

err_total_count = sum(symbol_errors_total);
err_stable_count = sum(symbol_errors_stable);

fprintf('  >> 全段 EVM (含收敛期): %.2f%%\n', evm_pct_total);
fprintf('  >> 稳态 EVM (Block 2~8): %.2f%%  <-- 【工业级验收标准】\n\n', evm_pct_stable);

fprintf('  >> 总错误数 (1~1024点): %d\n', err_total_count);
fprintf('  >> 稳态错误数 (129~1024点): %d\n', err_stable_count);

if err_stable_count == 0
    fprintf('  >> 🏆 完美！基带稳态区域实现 0 误码恢复！\n');
end
fprintf('======================================================\n');

%% 11. 终极绘图
figure('Name', '终极验收：基带瞬态与稳态相轨迹解剖', 'Position', [150, 200, 1200, 500]);

% 1. 稳态星座图
subplot(1, 3, 1);
plot(real(x_rx_pss_perfect(stable_idx)), imag(x_rx_pss_perfect(stable_idx)), 'b.', 'MarkerSize', 8); hold on;
plot(real(local_syms_aligned(stable_idx)), imag(local_syms_aligned(stable_idx)), 'r+', 'MarkerSize', 12, 'LineWidth', 2);
title(sprintf('稳态锁定星座图 (Block 2~8)\n(稳态 EVM: %.2f%%)', evm_pct_stable));
legend('补偿后接收符号', '理想物理符号', 'Location', 'best');
xlabel('I'); ylabel('Q');
axis square; grid on; xlim([-1.5 1.5]); ylim([-1.5 1.5]);

% 2. 逐点相位误差全景图 (展示收敛全过程)
subplot(1, 3, [2, 3]);
plot(1:N_syms, phase_error_deg, 'b-', 'LineWidth', 1.2); hold on;

yline(45, 'r--', '判决门限 (+45度)', 'LineWidth', 1.5, 'LabelHorizontalAlignment', 'left');
yline(-45, 'r--', '判决门限 (-45度)', 'LineWidth', 1.5, 'LabelHorizontalAlignment', 'left');
xline(128.5, 'k-', '环路锁定分界线 (Block 1 结束)', 'LineWidth', 2, 'LabelOrientation', 'horizontal');

% 标红错点
err_indices = find(symbol_errors_total);
if ~isempty(err_indices)
    plot(err_indices, phase_error_deg(err_indices), 'rx', 'MarkerSize', 10, 'LineWidth', 2);
end

title(sprintf('1024 个符号的逐点相位误差 (找回 Block 1 + 精准频偏补偿)\n稳态 EVM: %.2f%%, 稳态错误数: %d', evm_pct_stable, err_stable_count));
xlabel('物理符号索引 (1 to 1024)'); ylabel('绝对相位误差 (度)');
grid on; xlim([1 N_syms]); ylim([-90 90]);

%% 12. 终极探秘：开环相干积分侦破全段 8 Block 真实极性
fprintf('\n======================================================\n');
fprintf('     🕵️‍♂️ 终极探秘：开环相干积分侦破 8 Block 极性 🕵️‍♂️\n');
fprintf('======================================================\n');

polarities_inferred = zeros(1, 8);
polarities_inferred(1) = 1; % 假设第 1 块为 '+' 作为相对基准
polarities_str = '+';

corr_array = zeros(1, 8);
phase_deg_array = zeros(1, 8);
delta_phase_array = zeros(1, 8);

fprintf('区块 | 原始相干向量相位 | 相邻相位差 (Delta) | 推断极性判定\n');
fprintf('----------------------------------------------------------\n');

for k = 1:8
    % 从开环信号中切出第 k 个 Block
    start_idx = sym_start_idx + (k-1)*128;
    end_idx = start_idx + 127;
    
    % 容错防越界
    if end_idx > length(y_costas_in)
        warning('数据长度不足以提取第 %d 块', k);
        break;
    end
    
    block_raw = y_costas_in(start_idx : end_idx);

    % 相干积分 (内积)：将 128 个点浓缩成一个宏观的复数向量
    corr_val = sum(block_raw(:) .* conj(local_syms_single(:)));
    corr_array(k) = corr_val;
    
    % 计算绝对相位
    phase_deg_array(k) = rad2deg(angle(corr_val));
    
    if k == 1
        delta_phase_array(k) = 0;
        fprintf(' #%d  |    %8.2f 度   |        N/A        |   + (基准)\n', k, phase_deg_array(k));
    else
        % 利用共轭乘法计算相邻两块的相位差 (自动处理跨越 360 度的卷绕)
        delta_rad = angle(corr_array(k) * conj(corr_array(k-1)));
        delta_deg = rad2deg(delta_rad);
        delta_phase_array(k) = delta_deg;
        
        % 极性判定逻辑：
        % 如果相邻相位差在 90 度以内，说明只是残留频偏导致的轻微旋转，极性没变。
        % 如果相位差超过 90 度 (通常接近 180 度)，说明极性绝对翻转了。
        if abs(delta_deg) > 90
            polarities_inferred(k) = -polarities_inferred(k-1);
            pol_char = '- (翻转)';
        else
            polarities_inferred(k) = polarities_inferred(k-1);
            pol_char = '+ (相同)';
        end
        
        % 拼接字符串用于最终展示
        if polarities_inferred(k) == 1
            polarities_str = [polarities_str, ' +'];
        else
            polarities_str = [polarities_str, ' -'];
        end
        
        fprintf(' #%d  |    %8.2f 度   |    %8.2f 度   |   %s\n', k, phase_deg_array(k), delta_deg, pol_char);
    end
end

fprintf('----------------------------------------------------------\n');
fprintf('  >> 💡 物理层真实宏观极性推断结果: [ %s ]\n', polarities_str);
fprintf('======================================================\n');

%% 13. 汇报专用：不同极性假设的逐点相位误差终极对比 (汇报铁证)
fprintf('\n======================================================\n');
fprintf('     📊 第13步：汇报专用 - 极性假设对比铁证图 📊\n');
fprintf('======================================================\n');

% 提取已经完美去除频偏的接收信号
rx_signal_clean = x_rx_pss_perfect; 
steady_idx = 129:1024; % 稳态区索引 (Block 2 ~ 8)

% -----------------------------------------------------------
% 假设 A: 物理层真实推断 [- + + + + + + +]
% -----------------------------------------------------------
signs_A = [-1, 1, 1, 1, 1, 1, 1, 1];
template_A = kron(signs_A, local_syms_single);
template_A = template_A(:);

% 给假设 A 最公平的对齐机会 (用稳态区求平均相位差)
offset_A = angle(mean(rx_signal_clean(steady_idx) .* conj(template_A(steady_idx))));
template_A_aligned = template_A * exp(1j * offset_A);

% 计算假设 A 的相位误差
err_rad_A = angle(rx_signal_clean .* conj(template_A_aligned));
err_deg_A = rad2deg(err_rad_A);
err_count_A = sum(abs(err_rad_A(steady_idx)) > pi/4);

% -----------------------------------------------------------
% 假设 B: 初期猜测 [+ + - + - + - +]
% -----------------------------------------------------------
signs_B = [1, 1, -1, 1, -1, 1, -1, 1];
template_B = kron(signs_B, local_syms_single);
template_B = template_B(:);

% 给假设 B 同等公平的对齐机会
offset_B = angle(mean(rx_signal_clean(steady_idx) .* conj(template_B(steady_idx))));
template_B_aligned = template_B * exp(1j * offset_B);

% 计算假设 B 的相位误差
err_rad_B = angle(rx_signal_clean .* conj(template_B_aligned));
err_deg_B = rad2deg(err_rad_B);
err_count_B = sum(abs(err_rad_B(steady_idx)) > pi/4);

% -----------------------------------------------------------
% 绘制终极对比大图
% -----------------------------------------------------------
figure('Name', '汇报专用：极性假设对比铁证', 'Position', [100, 100, 1200, 800]);

% --- 子图 1: 假设 A ---
subplot(2, 1, 1);
plot(1:1024, err_deg_A, 'b-', 'LineWidth', 1.2); hold on;
yline(45, 'r--', '判决门限 (+45度)', 'LineWidth', 1.5);
yline(-45, 'r--', '判决门限 (-45度)', 'LineWidth', 1.5);
% 画出 Block 的分界线
for b = 1:9
    xline((b-1)*128 + 1, 'k:', 'LineWidth', 1);
end
% 标红稳态区的错点
err_idx_A = steady_idx(abs(err_rad_A(steady_idx)) > pi/4);
if ~isempty(err_idx_A)
    plot(err_idx_A, err_deg_A(err_idx_A), 'rx', 'MarkerSize', 8, 'LineWidth', 2);
end
title(sprintf('假设 A (真实物理推断): [ - + + + + + + + ]   |   稳态区错误数: %d', err_count_A), 'FontSize', 14);
ylabel('绝对相位误差 (度)', 'FontSize', 12);
grid on; xlim([1 1024]); ylim([-200 200]);
set(gca, 'YTick', -180:45:180);

% --- 子图 2: 假设 B ---
subplot(2, 1, 2);
plot(1:1024, err_deg_B, 'Color', [0.85 0.33 0.1], 'LineWidth', 1.2); hold on;
yline(45, 'r--', '判决门限 (+45度)', 'LineWidth', 1.5);
yline(-45, 'r--', '判决门限 (-45度)', 'LineWidth', 1.5);
for b = 1:9
    xline((b-1)*128 + 1, 'k:', 'LineWidth', 1);
end
% 标红稳态区的错点
err_idx_B = steady_idx(abs(err_rad_B(steady_idx)) > pi/4);
if ~isempty(err_idx_B)
    plot(err_idx_B, err_deg_B(err_idx_B), 'rx', 'MarkerSize', 8, 'LineWidth', 2);
end
title(sprintf('假设 B (初期极性猜想): [ + + - + - + - + ]   |   稳态区错误数: %d', err_count_B), 'FontSize', 14);
xlabel('物理符号索引 (1 to 1024)', 'FontSize', 12); 
ylabel('绝对相位误差 (度)', 'FontSize', 12);
grid on; xlim([1 1024]); ylim([-200 200]);
set(gca, 'YTick', -180:45:180);

fprintf('  >> 图表已生成！请查看两者的稳态区对比。\n');
fprintf('======================================================\n');

%% 14. 帧结构逆向：瞬态免疫法精准探测 PSS 前置 CP 长度
fprintf('\n======================================================\n');
fprintf('     📏 第14步：精准测量 PSS 前置 CP 长度 📏\n');
fprintf('======================================================\n');

% 向前探查的最大符号数 (假设 CP 不会超过 256)
max_lookback = 256; 
% 防止索引越界 (如果抓包刚好从信号开头开始，可能前面没有足够的数据)
lookback_syms = min(max_lookback, sym_start_idx - 1); 

if lookback_syms < 10
    fprintf('  >> ⚠️ 警告: PSS 起点前方的数据太少 (仅 %d 符号)，无法有效探测 CP。\n', lookback_syms);
else
    % 1. 提取 PSS 前方的待测数据 (使用开环数据 y_costas_in，免疫锁相环瞬态畸变)
    rx_front = y_costas_in(sym_start_idx - lookback_syms : sym_start_idx - 1);
    rx_front = rx_front(:);
    
    % 2. 提取 PSS 尾部的对应数据 (因为 CP 是尾部的拷贝，它们相距 N_syms = 1024)
    rx_tail = y_costas_in(sym_start_idx + N_syms - lookback_syms : sym_start_idx + N_syms - 1);
    rx_tail = rx_tail(:);
    
    % 3. 计算延迟共轭乘积的幅度 (能量)
    % 物理意义：无视相位旋转，只看波形是否高度相似。是 CP 则能量高，是噪声则能量低。
    cp_correlation_mag = abs(rx_front .* conj(rx_tail));
    
    % 为了图表好看，用尾部能量做个归一化
    tail_power = abs(rx_tail).^2;
    cp_metric_norm = cp_correlation_mag ./ (mean(tail_power) + 1e-6);

    % 4. 提取前端的原始信号幅度 (辅助判断突发信号真实的物理起点)
    front_amplitude = abs(rx_front);
    
    % 5. 可视化绘图
    idx_xaxis = (1 - lookback_syms) : 0; % 构造 X 轴：0 为 PSS 的边界
    
    figure('Name', 'PSS 前置 CP 长度精准探测', 'Position', [150, 150, 1000, 600]);
    
    % --- 子图 1: 信号绝对幅度 (看突发信号从哪里起振) ---
    subplot(2, 1, 1);
    plot(idx_xaxis, front_amplitude, 'Color', [0.6 0.6 0.6], 'LineWidth', 1.2); hold on;
    xline(0, 'r-', 'PSS 第一块起点 (-)', 'LineWidth', 2, 'LabelHorizontalAlignment', 'left');
    title('PSS 边界前方的原始信号幅度 (Envelope)');
    ylabel('幅度'); grid on; axis tight;
    xlim([-lookback_syms, 5]);
    
    % --- 子图 2: 延迟共轭相关度 (看 CP 的真实长度) ---
    subplot(2, 1, 2);
    plot(idx_xaxis, cp_metric_norm, 'b.-', 'LineWidth', 1.5, 'MarkerSize', 8); hold on;
    xline(0, 'r-', 'PSS 第一块起点 (-)', 'LineWidth', 2, 'LabelHorizontalAlignment', 'left');
    yline(0.5, 'k--', 'CP 判定门限 (参考)');
    title('PSS 前方数据与 PSS 尾部数据的复相关度 (距离=1024)');
    xlabel('相对于 PSS 起点的符号索引 (0 为 PSS 边界)');
    ylabel('归一化相关度'); grid on; axis tight;
    xlim([-lookback_syms, 5]); ylim([0, max(1.2, max(cp_metric_norm)*1.1)]);
    
    fprintf('  >> CP 探测图表已生成！\n');
    fprintf('  >> 请看下方的【复相关度】子图：从右向左看（从 0 开始往负数走）。\n');
    fprintf('  >> 蓝线维持在高位的平台宽度，就是精确的 CP 符号数！一旦跌落即为噪声。\n');
end
fprintf('======================================================\n');

%% 15. 终极实锤：CP 与 PSS 尾部的逐点波形比对 (DNA 级鉴定)
fprintf('\n======================================================\n');
fprintf('     🧬 第15步：CP 与尾部波形逐点“DNA 对比” 🧬\n');
fprintf('======================================================\n');

% 我们提取 PSS 前方 64 个符号 (包含 32 个预期的 CP 和 32 个纯噪声区，方便看对比边界)
test_window = 64;
cp_length_assumed = 32; % 我们要验证的界限

% 1. 提取 PSS 前方的 64 个样本 (前置区)
front_idx = (sym_start_idx - test_window) : (sym_start_idx - 1);
front_data = y_costas_in(front_idx);

% 2. 提取 PSS 尾部对应的 64 个样本 (尾部拷贝源)
tail_idx = (sym_start_idx + N_syms - test_window) : (sym_start_idx + N_syms - 1);
tail_data = y_costas_in(tail_idx);

% 3. 精准补偿时间跨度 (1024个符号) 带来的残余固定相位旋转
% 我们利用预期是 CP 的那 32 个点来计算平均相位差
cp_front_part = front_data(end - cp_length_assumed + 1 : end);
cp_tail_part  = tail_data(end - cp_length_assumed + 1 : end);
phase_offset = angle(mean(cp_front_part .* conj(cp_tail_part)));

% 将尾部数据进行相位对齐，使其与前端波形处于同一坐标系
tail_data_aligned = tail_data * exp(1j * phase_offset);

% 4. 计算逐点相位误差 (用于直观展示对齐度)
point_phase_diff = rad2deg(angle(front_data .* conj(tail_data_aligned)));

% 5. 绘图展示 (终极视觉暴击)
x_axis_cp = -test_window : -1;

figure('Name', 'CP 长度逐点鉴定：波形与相位对比', 'Position', [150, 100, 1100, 800]);

% --- 子图 1: I 支路波形叠加对比 ---
subplot(3, 1, 1);
plot(x_axis_cp, real(front_data), 'b-', 'LineWidth', 1.5); hold on;
plot(x_axis_cp, real(tail_data_aligned), 'r--', 'LineWidth', 1.5);
xline(-cp_length_assumed, 'k-', 'CP 理论边界 (-32)', 'LineWidth', 2);
title('I 支路 (实部) 波形逐点叠加对比');
legend('前置区域波形 (Front)', '相位对齐后的尾部波形 (Tail)', 'Location', 'best');
ylabel('I 支路幅度'); grid on; xlim([-test_window, 0]);

% --- 子图 2: Q 支路波形叠加对比 ---
subplot(3, 1, 2);
plot(x_axis_cp, imag(front_data), 'b-', 'LineWidth', 1.5); hold on;
plot(x_axis_cp, imag(tail_data_aligned), 'r--', 'LineWidth', 1.5);
xline(-cp_length_assumed, 'k-', 'CP 理论边界 (-32)', 'LineWidth', 2);
title('Q 支路 (虚部) 波形逐点叠加对比');
ylabel('Q 支路幅度'); grid on; xlim([-test_window, 0]);

% --- 子图 3: 逐点相位误差 ---
subplot(3, 1, 3);
plot(x_axis_cp, point_phase_diff, 'k.-', 'LineWidth', 1.2, 'MarkerSize', 8); hold on;
xline(-cp_length_assumed, 'k-', 'CP 理论边界 (-32)', 'LineWidth', 2);
yline(45, 'r:'); yline(-45, 'r:');
title('前置区域与尾部波形的逐点相位差 (Phase Difference)');
xlabel('相对于 PSS 起点的符号索引'); ylabel('相位差 (度)');
grid on; xlim([-test_window, 0]); ylim([-180 180]); set(gca, 'YTick', -180:45:180);

fprintf('  >> 波形比对图表已生成！\n');
fprintf('  >> 请观察 -32 点左右的波形贴合度与相位散布情况。\n');
fprintf('======================================================\n');

% --- 兜底读取函数 ---
function [iq_data, meta] = iq_read_int16_le(filename, start_idx, num_samples)
    fid = fopen(filename, 'rb');
    if fid == -1, error('无法打开文件，请检查 burst_1.dat 路径'); end
    fseek(fid, start_idx * 4, 'bof');
    raw_data = fread(fid, num_samples * 2, 'int16');
    fclose(fid);
    I = raw_data(1:2:end); Q = raw_data(2:2:end);
    iq_data = complex(double(I), double(Q));
    meta = [];
end