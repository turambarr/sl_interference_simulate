 clear; clc; close all;
% ===== 自动获取当前脚本所在路径，并切换过去 =====
if ~isdeployed
    currentFilePath = fileparts(mfilename('fullpath'));
    cd(currentFilePath);
    fprintf('已自动切换到脚本所在目录: %s\n', pwd);
end
% ==============================================

%% 1. 生成本地理想 PSS 序列 (使用实测序列与 48符号 CP)
fprintf('Step 1: 生成本地 PSS 序列...\n');
fs_symbol = 60e6;          
fs_target = 409.6e6;       
fc = 63e6;                 

% 解调出的m序列真实数据
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

[P, Q] = rat(fs_target / fs_symbol);
pss_up = resample(pss_base_with_cp, P, Q);
t_local = (0 : length(pss_up)-1) / fs_target;
pss_local = pss_up .* exp(1j * 2 * pi * fc * t_local);
pss_local = pss_local / max(abs(pss_local));

%% ========================================================================
%% [连续信号全局分析扫描区]
%% ========================================================================
%% 2. 读取包含多个连续 Burst 的实测数据
fprintf('Step 2: 读取包含多个连续 Burst 的实测数据...\n');
inFile = 'target_signal.dat'; % 确保文件名正确
fInfo = dir(inFile);
if isempty(fInfo)
    error('\n[致命错误]: 找不到文件 "%s"！请检查文件名是否正确。\n', inFile);
end
totalSamples = floor(fInfo.bytes / 4); 
readLen = min(totalSamples, 20e6); 
fid = fopen(inFile, 'r', 'ieee-le');
raw_data = fread(fid, readLen * 2, 'int16');
fclose(fid);
x_raw = double(raw_data(1:2:end)) + 1j * double(raw_data(2:2:end));
x_raw = x_raw - mean(x_raw);
x_raw = x_raw / mean(abs(x_raw));
fprintf('          成功读取 %d 个采样点 (约 %.2f 毫秒).\n', length(x_raw), length(x_raw)/fs_target*1000);

%% 2.8 寻找第一个物理能量起振点 (用于后续频偏盲估)
fprintf('Step 2.8: 寻找首个 Burst 的能量起振点...\n');
pwr = abs(x_raw).^2;
smoothed_pwr = movmean(pwr, 50);
noise_floor = mean(smoothed_pwr(1:min(2000, length(x_raw))));
steady_signal = max(smoothed_pwr); 
threshold_pwr = noise_floor + 0.1 * (steady_signal - noise_floor);
energy_start_idx = find(smoothed_pwr > threshold_pwr, 1, 'first');
if isempty(energy_start_idx), energy_start_idx = 1; end

%% 2.9 基于首个 Burst 进行频偏盲估与全局补偿
fprintf('Step 2.9: 提取首个 Burst 执行 M=4 次幂频偏盲估并全局补偿...\n');
cfo_est_len = min(20000, length(x_raw) - energy_start_idx);
if cfo_est_len > 4096
    x_for_cfo = x_raw(energy_start_idx + 1000 : energy_start_idx + 1000 + cfo_est_len - 1);
    M = 4; Nfft_cfo = 2^17; 
    Z_spec = fftshift(fft(x_for_cfo.^M, Nfft_cfo));
    f_axis_cfo = (-Nfft_cfo/2 : Nfft_cfo/2 - 1) / Nfft_cfo * fs_target;
    
    target_f = 63e6 * M;
    target_f_mapped = mod(target_f + fs_target/2, fs_target) - fs_target/2;
    [~, center_idx] = min(abs(f_axis_cfo - target_f_mapped));
    search_range = 200; 
    local_spec = abs(Z_spec(max(1, center_idx-search_range) : min(Nfft_cfo, center_idx+search_range)));
    [~, local_max_idx] = max(local_spec);
    actual_max_idx = center_idx - search_range + local_max_idx - 1;
    
    f_measured = f_axis_cfo(actual_max_idx);
    delta_f = (f_measured - target_f_mapped) / M;
    
    fprintf('          全局粗略频偏: %.2f kHz\n', delta_f / 1e3);
    
    t_global = (0:length(x_raw)-1).' / fs_target; 
    x_comp = x_raw .* exp(-1j * 2 * pi * delta_f * t_global);
else
    fprintf('          首个信号过短，跳过盲频偏补偿。\n');
    x_comp = x_raw;
    delta_f = 0; 
end

%% 3. 全局滑动互相关与严格归一化计算
fprintf('\nStep 3: 启动全局滑动归一化互相关扫描...\n');
h_matched = fliplr(conj(pss_local)); 
corr_out = filter(h_matched, 1, x_comp); 
corr_mag_sq = abs(corr_out).^2; 
E_local = sum(abs(pss_local).^2); 
win_ones = ones(1, length(pss_local));
E_rx_moving = filter(win_ones, 1, abs(x_comp).^2);
E_rx_moving(E_rx_moving < 1e-10) = 1e-10; 
corr_norm_pct_real = (corr_mag_sq ./ (E_local .* E_rx_moving)) * 100;

%% 4. 自动识别所有的 Burst 同步峰值
fprintf('Step 4: 提取所有 Burst 同步位置...\n');
min_peak_height = 50; 
min_peak_dist = length(pss_local); 
[pks_real, locs_real] = findpeaks(corr_norm_pct_real, 'MinPeakHeight', min_peak_height, 'MinPeakDistance', min_peak_dist);
num_frames = length(locs_real);
if num_frames == 0
    error('信号太短或未发现足够的 Burst！');
end
fprintf('          扫描完成！共探测到 %d 个有效的连续 Burst！\n', num_frames);

%% ========================================================================
%% 单频压制干扰注入 (CP干扰)
%% ========================================================================
%% 5. 对每个物理帧叠加 0dB 的单频压制干扰
fprintf('\nStep 5: [电子战攻击] 正在对每个真实帧执行 0dB 单频压制干扰注入...\n');
% 干扰波参数设定
pss_len = length(pss_local);
N_cp_samples = round(N_cp_syms * fs_target / fs_symbol);
f_jam_bb = fs_symbol / 128; % 468.75 kHz 特征单音
N_jam_syms = 1200;          % 干扰长度，保证覆盖整个 PSS
jam_len_samples = round(N_jam_syms * fs_target / fs_symbol);

% 获取真实包络幅度，设定 0dB 压制幅度
first_real_start = locs_real(1) - pss_len + 1;
real_burst_amp = max(abs(x_comp(max(1, first_real_start) : locs_real(1))));
jam_power_ratio = real_burst_amp * 1.0; % 1.0 即 0dB 等功率压制

fprintf('          目标包络幅度约为: %.2f | 单频压制幅度锁定为: %.2f\n', real_burst_amp, jam_power_ratio);

spoof_signal = zeros(size(x_raw)); % 存放生成的压制干扰

for i = 1:num_frames
    % 精确定位当前帧的 PSS 主体起始位置，用于锚定 0 度相位
    real_start_idx = locs_real(i) - pss_len + 1;
    pss_body_start_idx = real_start_idx + N_cp_samples;
    
    % 让干扰稍微提前一点开始，确保覆盖全长
    jam_start_idx = max(1, real_start_idx - round(50 * fs_target / fs_symbol));
    jam_end_idx = min(length(x_raw), jam_start_idx + jam_len_samples - 1);
    
    if jam_start_idx > 0 && jam_end_idx <= length(x_raw)
        % 【核心】相对时间轴，锚定到 pss_body_start_idx 时相位为 0
        t_rel = ( (jam_start_idx : jam_end_idx).' - pss_body_start_idx ) / fs_target;
        jam_bb = exp(1j * 2 * pi * f_jam_bb * t_rel); % 单频基带
        
        % 调制到射频载波并注入 (附加频偏补偿以挂载到物理信道)
        jam_rf = jam_power_ratio * jam_bb .* exp(1j * 2 * pi * (fc + delta_f) * t_global(jam_start_idx : jam_end_idx));
        
        spoof_signal(jam_start_idx : jam_end_idx) = spoof_signal(jam_start_idx : jam_end_idx) + jam_rf;
    end
end

% 污染信道：真实信号 + 压制信号
x_spoofed_raw = x_raw + spoof_signal;
fprintf('          成功注入所有单频压制干扰！\n');

%% 6. 重新扫描被压制干扰后的信道
fprintf('Step 6: 扫描 0dB 被压制后的信道以验证同步瘫痪...\n');
x_spoofed_comp = x_spoofed_raw .* exp(-1j * 2 * pi * delta_f * t_global);
corr_out_spoofed = filter(h_matched, 1, x_spoofed_comp); 
corr_mag_sq_spoofed = abs(corr_out_spoofed).^2; 
E_rx_moving_spoofed = filter(win_ones, 1, abs(x_spoofed_comp).^2);
E_rx_moving_spoofed(E_rx_moving_spoofed < 1e-10) = 1e-10; 
corr_norm_pct_spoofed = (corr_mag_sq_spoofed ./ (E_local .* E_rx_moving_spoofed)) * 100;

% 获取被干扰后的所有帧残存峰值高度
pks_jammed = corr_norm_pct_spoofed(locs_real);

%% 7. 绘制干扰效果对比图
fprintf('Step 7: 绘制单频压制效果对比图...\n');
figure('Position', [100, 100, 1400, 800], 'Name', 'Electronic Warfare: Single-Tone Suppression Jamming');
t_axis = (0 : length(corr_norm_pct_real)-1) * (1 / fs_target) * 1e3; 

% --- 上图：攻击前的纯净信道 ---
subplot(2, 1, 1);
plot(t_axis, corr_norm_pct_real, 'Color', [0.4 0.6 0.8], 'LineWidth', 1.5); hold on;
plot(t_axis(locs_real), pks_real, 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 8);
yline(min_peak_height, 'r--', '检测门限 (50%)', 'LineWidth', 1.5, 'LabelHorizontalAlignment', 'right');
title('【压制前】无干扰环境下的连续目标信号归一化相关峰');
ylabel('归一化相似度 (%)'); grid on; xlim([0, t_axis(end)]); ylim([0 110]);

% 在原图上标记真实高度
for i = 1:num_frames
    text(t_axis(locs_real(i)), pks_real(i) + 6, sprintf('%.1f%%', pks_real(i)), 'Color', 'b', 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
end

% --- 下图：单频攻击后的污染信道 ---
subplot(2, 1, 2);
plot(t_axis, corr_norm_pct_spoofed, 'Color', [0.8 0.4 0.0], 'LineWidth', 1.5); hold on;
plot(t_axis(locs_real), pks_jammed, 'rx', 'MarkerSize', 10, 'LineWidth', 2);
yline(min_peak_height, 'r--', '检测门限 (50%)', 'LineWidth', 1.5, 'LabelHorizontalAlignment', 'right');
title('【压制后】注入 0dB 等功率单频干扰 (CP干扰) 后的同步瘫痪结果');
xlabel('时间 Time (ms)'); ylabel('归一化相似度 (%)'); grid on; xlim([0, t_axis(end)]); ylim([0 110]);

% 在下图标注残存高度
for i = 1:num_frames
    text(t_axis(locs_real(i)), pks_jammed(i) + 6, sprintf('残存峰\n%.1f%%', pks_jammed(i)), 'Color', 'r', 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
end

%% ========================================================================
%% 8. 【核心新增】批量生成并保存不同干信比的污染信号文件 (单频压制干扰)
%% ========================================================================
fprintf('\nStep 8: 正在分别生成 -3dB, 0dB, +3dB 干信比的靶标文件...\n');
JSR_list = [-3, 0, 3];

for idx = 1:length(JSR_list)
    JSR_dB = JSR_list(idx);
    
    % 将 JSR (dB) 转换为线性幅度倍数
    amp_ratio = 10^(JSR_dB / 20);
    
    % 构造当前 JSR 下的受污染原始射频信号 (物理域叠加)
    % 这里的 spoof_signal 在 Step 5 中已经是 0dB 状态，因此直接乘上 amp_ratio 即可
    current_jammed_raw = x_raw + spoof_signal * amp_ratio;
    
    % 动态范围归一化与 int16 量化防溢出 (留 5% 裕量)
    max_val = max(abs(current_jammed_raw));
    x_out_scaled = (current_jammed_raw / max_val) * 32767 * 0.95; 
    
    % 拆分为交织的 I/Q 格式
    out_data = zeros(2 * length(x_out_scaled), 1, 'int16');
    out_data(1:2:end) = int16(round(real(x_out_scaled)));
    out_data(2:2:end) = int16(round(imag(x_out_scaled)));
    
    % 动态生成文件名，替换原文件的 .dat 后缀 (加入 CP_jammed 标识)
    [filepath, name, ext] = fileparts(inFile);
    outFile = fullfile(filepath, sprintf('%s_multiframe_jammed_%ddB%s', name, JSR_dB, ext));
    
    % 小端序写入
    fid_out = fopen(outFile, 'w', 'ieee-le');
    fwrite(fid_out, out_data, 'int16');
    fclose(fid_out);
    
    fprintf('          ✅ 已成功保存靶标文件: %s\n', outFile);
end

fprintf('================================================================\n');
fprintf('🏁 压制式单频干扰 (CP干扰) 评估及靶标生成全流程执行完毕！\n');
fprintf('================================================================\n');

%% ========================================================================
%% 9. 【新增】批量下变频靶标信号至基带，并保存为 R&S .wv 格式
%% ========================================================================
fprintf('\nStep 9: [基带靶标] 正在批量下变频并生成 .wv 格式基带文件...\n');

% 目标下变频中心频率 (Hz) - 保持精准频偏
f_down = 63132812; 

% 遍历 Step 8 中设定的干信比列表
for idx = 1:length(JSR_list)
    JSR_dB = JSR_list(idx);
    
    % 将 JSR (dB) 转换为线性幅度倍数
    amp_ratio = 10^(JSR_dB / 20);
    
    % 重新构造当前 JSR 下的受污染原始射频信号 (中频)
    % 注：spoof_signal 在 Step 5 中已经是 0dB 全局信号
    current_jammed_raw = x_raw + spoof_signal * amp_ratio;
    
    % 1. 执行数字下变频 (DDC)
    % 乘以负频偏复指数，将中频信号平移至 0 Hz (纯基带)
    baseband_jammed = current_jammed_raw .* exp(-1j * 2 * pi * f_down * t_global);
    
    % 2. 动态范围归一化与 16-bit 量化 (留 5% 裕量防削顶)
    max_val_bb = max(abs(baseband_jammed));
    baseband_scaled = (baseband_jammed / max_val_bb) * 32767 * 0.95;
    
    % 拆分为交织的 I/Q 格式 [I, Q, I, Q...]
    out_data_bb = zeros(2 * length(baseband_scaled), 1, 'int16');
    out_data_bb(1:2:end) = int16(round(real(baseband_scaled)));
    out_data_bb(2:2:end) = int16(round(imag(baseband_scaled)));
    
    % 3. 动态生成文件名并保存为标准 R&S .wv 波形文件
    [filepath, name, ~] = fileparts(inFile);
    wvFile = fullfile(filepath, sprintf('%s_multiframe_jammed_%ddB_baseband.wv', name, JSR_dB));
    
    fid_wv = fopen(wvFile, 'w');
    if fid_wv == -1
        error('无法创建文件 %s', wvFile);
    end
    
    % --- 写入 R&S .wv 必备的 ASCII 文件头 ---
    fprintf(fid_wv, '{TYPE: SMU-WV,0}');
    fprintf(fid_wv, '{COMMENT: Downconverted Multi-Frame Jammed Signal. Shift: -%d Hz}', f_down);
    fprintf(fid_wv, '{CLOCK: %d}', fs_target); % 采样率 409.6 MHz
    
    num_samples = length(baseband_scaled);
    fprintf(fid_wv, '{SAMPLES: %d}', num_samples);
    
    % WAVEFORM 标签：字节数 = 样本数 * 4 (I/Q各2字节) + 1 (#字符)
    waveform_bytes = num_samples * 4 + 1; 
    fprintf(fid_wv, '{WAVEFORM-%d: #}', waveform_bytes);
    
    % --- 写入小端序的 16-bit I/Q 二进制波形数据 ---
    fwrite(fid_wv, out_data_bb, 'int16', 'ieee-le');
    fclose(fid_wv);
    
    fprintf('          ✅ [%+3d dB] 多帧基带 .wv 靶标已保存: %s\n', JSR_dB, wvFile);
end

fprintf('================================================================\n');
fprintf('🏁 压制式单频干扰 (CP干扰) 评估及多帧基带靶标生成全流程执行完毕！\n');

%% ========================================================================
%% 10. 【新增】将多帧干扰前 (Clean) 和各干信比干扰后 (Jammed) 的信号下变频并保存为 .iq
%% ========================================================================
fprintf('\nStep 10: 正在将多帧干扰前后的信号下变频至基带并保存为纯二进制 .iq 格式...\n');

% 目标下变频中心频率 (Hz)
f_down = 63132812; 
[filepath, name, ~] = fileparts(inFile);

% --- 1. 处理并保存干扰前的干净信号 (Clean) ---
fprintf('          正在处理干扰前的原始信号 (Clean)...\n');

% 乘以负频偏复指数，将干净信号中心频率向左平移，降至 0 Hz
baseband_clean = x_raw .* exp(-1j * 2 * pi * f_down * t_global);

% 动态范围归一化与 16-bit 量化 (预留 5% 裕量防削顶)
max_val_clean = max(abs(baseband_clean));
baseband_clean_scaled = (baseband_clean / max_val_clean) * 32767 * 0.95;

% 拆分为交织的 I/Q 格式 (I_0, Q_0, I_1, Q_1...)
out_data_clean = zeros(2 * length(baseband_clean_scaled), 1, 'int16');
out_data_clean(1:2:end) = int16(round(real(baseband_clean_scaled)));
out_data_clean(2:2:end) = int16(round(imag(baseband_clean_scaled)));

iqFile_clean = fullfile(filepath, sprintf('%s_multiframe_baseband_clean.iq', name));
fid_iq_clean = fopen(iqFile_clean, 'w', 'ieee-le'); % 强制指定小端序
if fid_iq_clean == -1
    error('无法创建文件 %s', iqFile_clean);
end
fwrite(fid_iq_clean, out_data_clean, 'int16');
fclose(fid_iq_clean);
fprintf('          ✅ [Clean] 下变频完成！已保存至: %s\n', iqFile_clean);

% --- 2. 批量处理并保存各干信比的干扰信号 (Jammed) ---
fprintf('          正在批量处理各干信比的受干扰信号 (Jammed)...\n');
for idx = 1:length(JSR_list)
    JSR_dB = JSR_list(idx);
    
    % 计算干扰幅度倍数并构造当前射频污染信号 (中频)
    amp_ratio = 10^(JSR_dB / 20);
    % 注意：多帧脚本中，spoof_signal 在 Step 5 已是全局 0dB 状态
    current_jammed_raw = x_raw + spoof_signal * amp_ratio;
    
    % 执行数字下变频 (DDC)
    baseband_jammed = current_jammed_raw .* exp(-1j * 2 * pi * f_down * t_global);
    
    % 动态范围归一化与 16-bit 量化
    max_val_bb = max(abs(baseband_jammed));
    baseband_scaled = (baseband_jammed / max_val_bb) * 32767 * 0.95;
    
    % 拆分为交织的 I/Q 格式
    out_data_bb = zeros(2 * length(baseband_scaled), 1, 'int16');
    out_data_bb(1:2:end) = int16(round(real(baseband_scaled)));
    out_data_bb(2:2:end) = int16(round(imag(baseband_scaled)));
    
    % 动态生成文件名并保存
    iqFile_jammed = fullfile(filepath, sprintf('%s_multiframe_jammed_%ddB_baseband.iq', name, JSR_dB));
    fid_iq = fopen(iqFile_jammed, 'w', 'ieee-le');
    if fid_iq == -1
        error('无法创建文件 %s', iqFile_jammed);
    end
    fwrite(fid_iq, out_data_bb, 'int16');
    fclose(fid_iq);
    
    fprintf('          ✅ [%+3d dB] 下变频完成！已保存至: %s\n', JSR_dB, iqFile_jammed);
end

fprintf('================================================================\n');
fprintf('🏁 压制式单频干扰多帧基带 .iq 文件转换全流程执行完毕！\n');

%% ========================================================================
%% 11. 【新增】批量保存中频 (IF) 的多帧被干扰信号 (-6dB 至 -16dB) 为 .iq 文件
%% ========================================================================
fprintf('\nStep 11: 正在批量保存中频 (IF) 的多帧被干扰信号为纯二进制 .iq 格式...\n');

% 定义需要额外生成的干信比列表
JSR_list_IF_extra = [-6, -8, -10, -12, -14, -16];

for idx = 1:length(JSR_list_IF_extra)
    JSR_target_dB_IF = JSR_list_IF_extra(idx);
    
    % 将 JSR (dB) 转换为线性电压放大系数
    amp_ratio_IF = 10^(JSR_target_dB_IF / 20);

    % 物理域叠加，构造中频被干扰信号 
    % (不执行下变频，保持 409.6MHz 采样率和原始 63MHz 载波)
    % 注：spoof_signal 在 Step 5 已是全局多帧的单频压制干扰
    if_jammed_signal = x_raw + spoof_signal * amp_ratio_IF;

    % 动态范围归一化与 16-bit 量化 (预留 5% 裕量防削顶)
    max_val_if = max(abs(if_jammed_signal));
    if_scaled = (if_jammed_signal / max_val_if) * 32767 * 0.95;

    % 拆分为交织的 I/Q 格式 (I_0, Q_0, I_1, Q_1...)
    out_data_if = zeros(2 * length(if_scaled), 1, 'int16');
    out_data_if(1:2:end) = int16(round(real(if_scaled)));
    out_data_if(2:2:end) = int16(round(imag(if_scaled)));

    % 动态生成文件名，并加上 _IF_jammed 标识
    [filepath, name, ~] = fileparts(inFile);
    iqFile_IF = fullfile(filepath, sprintf('%s_multiframe_IF_jammed_%ddB.iq', name, JSR_target_dB_IF));

    % 创建并以 IEEE 小端序写入纯二进制数据
    fid_iq_IF = fopen(iqFile_IF, 'w', 'ieee-le'); 
    if fid_iq_IF == -1
        error('无法创建文件 %s', iqFile_IF);
    end

    % 直接写入 16-bit 交织二进制波形数据
    fwrite(fid_iq_IF, out_data_if, 'int16');
    fclose(fid_iq_IF);

    fprintf('          ✅ 中频 (IF) 多帧被干扰信号 [%3ddB] 已成功保存至: %s\n', JSR_target_dB_IF, iqFile_IF);
end
fprintf('================================================================\n');
fprintf('🏁 扩展干信比 (-6 ~ -16dB) 多帧中频 .iq 文件批量生成完毕！\n');

%% ========================================================================
%% 12. 【新增】验证环节：回读中频靶标文件并验证压制干扰是否生效
%% ========================================================================
fprintf('\nStep 12: 正在回读保存的中频 .iq 文件，并验证压制效果...\n');

% 选择一个刚保存的文件进行验证 (这里以 -10dB 为例，您可以随意更改)
verify_JSR = -10;
[filepath, name, ~] = fileparts(inFile);
verifyFile = fullfile(filepath, sprintf('%s_multiframe_IF_jammed_%ddB.iq', name, verify_JSR));

fprintf('          正在读取并验证文件: %s\n', verifyFile);

if ~isfile(verifyFile)
    error('找不到要验证的文件: %s。请确认前面的保存步骤已成功执行。', verifyFile);
end

% --- 1. 读取文件并恢复复数格式 ---
fid_verify = fopen(verifyFile, 'r', 'ieee-le');
raw_verify = fread(fid_verify, Inf, 'int16');
fclose(fid_verify);

x_verify_raw = double(raw_verify(1:2:end)) + 1j * double(raw_verify(2:2:end));
fprintf('          ✅ 成功回读 %d 个复数样点。\n', length(x_verify_raw));

% --- 2. 频偏补偿 (模拟接收机侧的处理前提) ---
% 因为保存的是原始含有频偏的中频信号，需要先补偿频偏才能做完美互相关
t_global_v = (0:length(x_verify_raw)-1).' / fs_target;
x_verify_comp = x_verify_raw .* exp(-1j * 2 * pi * delta_f * t_global_v);

% --- 3. 重新执行滑动互相关，验证压制效果 ---
corr_out_v = filter(h_matched, 1, x_verify_comp); 
corr_mag_sq_v = abs(corr_out_v).^2; 
E_rx_moving_v = filter(win_ones, 1, abs(x_verify_comp).^2);
E_rx_moving_v(E_rx_moving_v < 1e-10) = 1e-10; 
corr_norm_pct_v = (corr_mag_sq_v ./ (E_local .* E_rx_moving_v)) * 100;

% 提取真实位置处的相关峰残留高度
pks_verify_jammed = corr_norm_pct_v(locs_real);

% --- 4. 绘图直观验证 (频域特征 + 同步瘫痪特征) ---
figure('Position', [150, 150, 1200, 800], 'Name', sprintf('Verification: %ddB Jamming Effect', verify_JSR));

% (1) 频域验证：功率谱密度 PSD
subplot(2, 1, 1);
nfft = 8192;
[pxx, f_axis_v] = pwelch(x_verify_raw, hanning(1024), 512, nfft, fs_target, 'centered');
plot(f_axis_v / 1e6, 10*log10(pxx), 'LineWidth', 1.5, 'Color', [0.3 0.2 0.6]);
title(sprintf('回读文件 [%ddB] 的功率谱密度 (PSD) - 验证中频载波与干扰特征', verify_JSR));
xlabel('频率 (MHz)'); ylabel('功率密度 (dB/Hz)');
grid on; axis tight;
xline(63, 'r--', '标称载频 63MHz', 'LabelVerticalAlignment', 'bottom', 'LineWidth', 1.5);
if f_axis_v(1)/1e6 < 0
    xlim([min(f_axis_v)/1e6, max(f_axis_v)/1e6]);
else
    xlim([0, fs_target/2/1e6]);
end

% (2) 效果验证：互相关同步峰
subplot(2, 1, 2);
t_axis_v = (0 : length(corr_norm_pct_v)-1) * (1 / fs_target) * 1e3; 
plot(t_axis_v, corr_norm_pct_v, 'Color', [0.8 0.4 0.0], 'LineWidth', 1.5); hold on;

% 标记残存峰
plot(t_axis_v(locs_real), pks_verify_jammed, 'rx', 'MarkerSize', 10, 'LineWidth', 2);
yline(min_peak_height, 'r--', '检测门限 (50%)', 'LineWidth', 1.5, 'LabelHorizontalAlignment', 'right');

title(sprintf('【效果验证】读取 [%ddB] 中频文件后的同步接收结果', verify_JSR));
xlabel('时间 Time (ms)'); ylabel('归一化相似度 (%)'); 
grid on; xlim([0, t_axis_v(end)]); ylim([0 110]);

% 在图标注残存高度
for i = 1:num_frames
    text(t_axis_v(locs_real(i)), pks_verify_jammed(i) + 8, ...
         sprintf('残存\n%.1f%%', pks_verify_jammed(i)), ...
         'Color', 'r', 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
end

fprintf('          📈 已绘制验证结果！\n');
fprintf('          观察图2，若残存峰值低于门限 (50%%) 或相较于纯净环境大幅下降，则证明 %ddB 压制成功！\n', verify_JSR);
fprintf('================================================================\n');