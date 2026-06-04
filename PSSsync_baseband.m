%% ========================================================================
%% 零中频 (基带) 信号同步独立分析脚本 (完美修复版)
%% ========================================================================
clear; clc; close all;

% ===== 自动获取当前脚本所在路径，并切换过去 =====
if ~isdeployed
    currentFilePath = fileparts(mfilename('fullpath'));
    cd(currentFilePath);
    fprintf('已自动切换到脚本所在目录: %s\n', pwd);
end

%% 1. 生成本地理想 PSS 纯基带序列
fprintf('\nStep 1: 生成本地 PSS 纯基带序列...\n');
fs_symbol = 60e6;          
fs_target = 409.6e6;       
fc_nominal = 0; % 标称中心频率

% 解调出的 m 序列真实数据
hex_str = 'BD3CD0148871751F84CED8C1BE32AC96'; 
hex_chars = char(hex_str);
bits = zeros(1, 128);
for i = 1:length(hex_chars)
    bin_str = dec2bin(hex2dec(hex_chars(i)), 4);
    bits((i-1)*4 + 1 : i*4) = bin_str - '0'; 
end

d_phi = (bits == 1) * (-pi/2) + (bits == 0) * (pi/2);
phase_accum = cumsum(d_phi);
m_syms = exp(1j * phase_accum);

% -+++++++ 的极性结构
pss_base = [-m_syms, m_syms, m_syms, m_syms, m_syms, m_syms, m_syms, m_syms];

% 48 符号与 CP 取反结构
N_cp_syms = 48; 
cp_syms = -pss_base(end - N_cp_syms + 1 : end);
pss_base_with_cp = [cp_syms, pss_base];

% 重采样至目标采样率 (得到纯基带波形，不挂载射频载波)
[P, Q] = rat(fs_target / fs_symbol);
pss_up = resample(pss_base_with_cp, P, Q);

% 归一化本地基带序列
pss_local_bb = pss_up / max(abs(pss_up));
pss_local_bb = pss_local_bb(:); % 确保为列向量 (Nx1)

fprintf('          成功生成基带序列，长度: %d 样点\n', length(pss_local_bb));

%% 2. 读取实测数据
fprintf('\nStep 2: 读取原始实测数据...\n');
inFile = 'baseband_sync_output.iq'; 
fInfo = dir(inFile);
if isempty(fInfo)
    error('\n[致命错误]: 找不到文件 "%s"！请检查文件是否存在。\n', inFile);
end
totalSamples = floor(fInfo.bytes / 4); 
readLen = min(totalSamples, 20e6); 
fid = fopen(inFile, 'r', 'ieee-le');
raw_data = fread(fid, readLen * 2, 'int16');
fclose(fid);

x_raw = double(raw_data(1:2:end)) + 1j * double(raw_data(2:2:end));
x_raw = x_raw - mean(x_raw);
x_raw = x_raw / mean(abs(x_raw));
x_raw = x_raw(:); % 确保为列向量 (Mx1)

fprintf('          成功读取 %d 个采样点 (约 %.2f 毫秒).\n', length(x_raw), length(x_raw)/fs_target*1000);

%% 3. 寻找能量起振点与频偏盲估
fprintf('\nStep 3: 寻找首帧能量起振点并执行频偏盲估...\n');
% --- 能量检测 ---
pwr = abs(x_raw).^2;
smoothed_pwr = movmean(pwr, 50);
noise_floor = mean(smoothed_pwr(1:min(2000, length(x_raw))));
steady_signal = max(smoothed_pwr); 
threshold_pwr = noise_floor + 0.1 * (steady_signal - noise_floor);
energy_start_idx = find(smoothed_pwr > threshold_pwr, 1, 'first');
if isempty(energy_start_idx), energy_start_idx = 1; end

% --- M=4 次幂频偏盲估 ---
cfo_est_len = min(20000, length(x_raw) - energy_start_idx);
if cfo_est_len > 4096
    x_for_cfo = x_raw(energy_start_idx + 1000 : energy_start_idx + 1000 + cfo_est_len - 1);
    M = 4; Nfft_cfo = 2^17; 
    Z_spec = fftshift(fft(x_for_cfo.^M, Nfft_cfo));
    f_axis_cfo = (-Nfft_cfo/2 : Nfft_cfo/2 - 1) / Nfft_cfo * fs_target;
    
    target_f = fc_nominal * M;
    target_f_mapped = mod(target_f + fs_target/2, fs_target) - fs_target/2;
    [~, center_idx] = min(abs(f_axis_cfo - target_f_mapped));
    search_range = 200; 
    local_spec = abs(Z_spec(max(1, center_idx-search_range) : min(Nfft_cfo, center_idx+search_range)));
    [~, local_max_idx] = max(local_spec);
    actual_max_idx = center_idx - search_range + local_max_idx - 1;
    
    f_measured = f_axis_cfo(actual_max_idx);
    delta_f = (f_measured - target_f_mapped) / M;
    fprintf('          侦测到全局残余频偏 (CFO): %.2f kHz\n', delta_f / 1e3);
else
    delta_f = 0;
    fprintf('          信号过短，跳过盲频偏估计。\n');
end

%% 4. 【核心】下变频至纯基带 (包含频偏补偿)
fprintf('\nStep 4: 正在进行带有频偏补偿的基带下变频...\n');
t_global = (0:length(x_raw)-1).' / fs_target;

% 将动态频偏与标称频率合并，一步到位搬移至 0 Hz 完美基带
f_down_actual = fc_nominal + delta_f; 
baseband_rx = x_raw .* exp(-1j * 2 * pi * f_down_actual * t_global);

% ================= 新增：提取并输出精确下变频参数 =================
Actual_Downconvert_Hz = f_down_actual; 
fprintf('          [参数导出] 目标信号标称载波: %d Hz\n', fc_nominal);
fprintf('          [参数导出] 实测残余频偏值: %+.3f Hz\n', delta_f);
fprintf('          [参数导出] 实际执行下变频: %.3f Hz (%.6f MHz)\n', Actual_Downconvert_Hz, Actual_Downconvert_Hz / 1e6);
% =================================================================

%% 5. 纯基带匹配滤波与归一化计算
fprintf('\nStep 5: 启动基带滑动归一化互相关扫描...\n');

% 【关键修复点】：因为 pss_local_bb 是列向量，必须使用 flipud (上下翻转) 
% 进行时间反褶，否则卷积退化为普通的互相关，能量无法积攒！
h_matched_bb = flipud(conj(pss_local_bb));            

% 核心滤波与能量计算
corr_out_bb = filter(h_matched_bb, 1, baseband_rx);   
corr_mag_sq_bb = abs(corr_out_bb).^2;                 

% 能量严格归一化 (防止幅度波动干扰判决)
E_local_bb = sum(abs(pss_local_bb).^2); 
win_ones_bb = ones(length(pss_local_bb), 1); % 这里也改为列向量适配
E_rx_moving_bb = filter(win_ones_bb, 1, abs(baseband_rx).^2);
E_rx_moving_bb(E_rx_moving_bb < 1e-10) = 1e-10; % 防除零保护      

% 计算 0~100% 的归一化相似度
corr_norm_pct_bb = (corr_mag_sq_bb ./ (E_local_bb .* E_rx_moving_bb)) * 100;

%% 6. 自动识别基带同步峰值并绘图
fprintf('\nStep 6: 提取同步位置并绘图...\n');
min_peak_height_bb = 50; 
min_peak_dist_bb = length(pss_local_bb); 
[pks_bb, locs_bb] = findpeaks(corr_norm_pct_bb, 'MinPeakHeight', min_peak_height_bb, 'MinPeakDistance', min_peak_dist_bb);

fprintf('          扫描完成！共探测到 %d 个有效的连续 Burst！\n', length(locs_bb));
for i = 1:length(locs_bb)
    fprintf('            - Burst %02d: 匹配度 %.2f%%, 位置(Sample) %d\n', i, pks_bb(i), locs_bb(i));
end

% 绘制波形图
figure('Position', [150, 150, 1200, 400], 'Name', 'Baseband Synchronization');
t_axis_bb = (0 : length(corr_norm_pct_bb)-1) * (1 / fs_target) * 1e3; 
plot(t_axis_bb, corr_norm_pct_bb, 'Color', [0.1 0.5 0.8], 'LineWidth', 1.5); hold on;
if ~isempty(locs_bb)
    plot(t_axis_bb(locs_bb), pks_bb, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 6);
end
yline(min_peak_height_bb, 'k--', '检测门限', 'LabelHorizontalAlignment', 'left');
title('纯基带 (Zero-IF) 下的归一化互相关同步结果');
xlabel('时间 Time (ms)'); ylabel('归一化相似度 (%)'); 
grid on; xlim([0, t_axis_bb(end)]); ylim([0 110]);
fprintf('\n🏁 基带同步分析执行完毕！\n');