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

% 👑 坚守解调出的第一行真实数据
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

% 👑 坚守 -+++++++ 的极性结构
pss_base = [-m_syms, m_syms, m_syms, m_syms, m_syms, m_syms, m_syms, m_syms];

% 👑 坚守 48 符号与 CP 取反结构
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

% ⚠️ 【请在这里修改为你的连续信号文件名】 ⚠️
inFile = 'sigtest1.iq'; 
startSample = 0;      % 0-based 复采样点起始位置
readLenReq  = 20e6;   % 读取复采样点数（可按需调整）

fInfo = dir(inFile);
if isempty(fInfo)
    error('\n[致命错误]: 找不到文件 "%s"！请检查文件名是否正确。\n', inFile);
end

% 按你的 .iq 读取流程（int16 LE, I/Q 交织）
bytesPerComplexSample = 4;
totalSamples = floor(fInfo.bytes / bytesPerComplexSample);

if startSample < 0 || startSample >= totalSamples
    error('startSample 越界：start=%d, total=%d', startSample, totalSamples);
end

maxReadable = totalSamples - startSample;
readLen = min(round(readLenReq), maxReadable);
if readLen <= 0
    error('readLen 非法：%d', readLen);
end

[x_raw, meta] = iq_read_int16_le(inFile, startSample, readLen);
if meta.numSamplesRead <= 0
    error('未读取到有效 IQ 数据。');
end

x_raw = double(x_raw(:));
x_raw = x_raw - mean(x_raw);
x_raw = x_raw / mean(abs(x_raw));

fprintf('          成功读取 %d 个采样点 (Start=%d, Req=%d, 实读=%d, 约 %.2f 毫秒).\n', ...
    length(x_raw), startSample, round(readLenReq), meta.numSamplesRead, length(x_raw)/fs_target*1000);

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
    
    t_global = (0:length(x_raw)-1).' / fs_target; % 注意这里改为列向量适配
    x_comp = x_raw .* exp(-1j * 2 * pi * delta_f * t_global);
else
    fprintf('          首个信号过短，跳过盲频偏补偿。\n');
    x_comp = x_raw;
end

%% 3. 全局滑动互相关与严格归一化计算 (雷达级扫描)
fprintf('\nStep 3: 启动全局滑动归一化互相关扫描 (雷达全开，请稍候)...\n');

% [第一步]：计算未经归一化的原始互相关
h_matched = fliplr(conj(pss_local)); 
corr_out = filter(h_matched, 1, x_comp); 
corr_mag_sq = abs(corr_out).^2; 

% [第二步]：计算滑动归一化因子
E_local = sum(abs(pss_local).^2); 
win_ones = ones(1, length(pss_local));
E_rx_moving = filter(win_ones, 1, abs(x_comp).^2);
E_rx_moving(E_rx_moving < 1e-10) = 1e-10; 

% [第三步]：得到全局 0 到 1 之间的绝对相似度系数
corr_norm = corr_mag_sq ./ (E_local .* E_rx_moving);
corr_norm_pct = corr_norm * 100;

%% 4. 自动识别所有的 Burst 同步峰值
fprintf('Step 4: 提取所有 Burst 同步位置...\n');

% 设置寻峰参数：
% 1. 相似度必须大于 50%
% 2. 两个 Burst 之间至少间隔 1 个 PSS 的长度 (防止同一个峰的旁瓣被误判)
min_peak_height = 50; 
min_peak_dist = length(pss_local); 

[pks, locs] = findpeaks(corr_norm_pct, 'MinPeakHeight', min_peak_height, 'MinPeakDistance', min_peak_dist);

% locs 对应匹配模板末端位置；同时计算模板起始位置（全局复采样点索引，0-based）
L_pss = length(pss_local);
end_pos_global = startSample + (locs - 1);
start_pos_global = end_pos_global - (L_pss - 1);

fprintf('          🎯 扫描完成！共探测到 %d 个有效的连续 Burst！\n', length(pks));
for i = 1:length(pks)
    fprintf('            - Burst %02d: 匹配度 %.2f%%, 末端=%d, 起始=%d\n', ...
        i, pks(i), end_pos_global(i), start_pos_global(i));
end

%% 5. 绘制全局归一化相关谱与时域包络对比图
fprintf('Step 5: 绘制连续 Burst 全局同步结果图...\n');

figure('Position', [100, 100, 1400, 600], 'Name', 'Global Continuous Bursts Synchronization');
sample_axis = startSample + (0 : length(corr_norm_pct)-1); % X轴：复采样点索引

% ==========================================
% 上子图：接收信号的原始包络能量 (用于直观看到信号存在)
% ==========================================
subplot(2, 1, 1);
plot(sample_axis, sqrt(smoothed_pwr), 'Color', [0.6 0.6 0.6], 'LineWidth', 1);
title(['连续信号包络图（复采样点）']);
ylabel('幅度 (Amplitude)');
grid on; axis tight;
xlim([sample_axis(1), sample_axis(end)]);

% ==========================================
% 下子图：全局滑动归一化相似度 (雷达寻峰图)
% ==========================================
subplot(2, 1, 2);
plot(sample_axis, corr_norm_pct, 'b', 'LineWidth', 1.2); 
hold on;

% 在所有识别到的峰值上打上红色的星星标记
if ~isempty(locs)
    plot(sample_axis(locs), pks, 'r*', 'MarkerSize', 10, 'LineWidth', 1.5);
    
    % 在每个星星上方标注匹配度
    for i = 1:length(pks)
        text(sample_axis(locs(i)), pks(i) + 5, sprintf('%.1f%%', pks(i)), ...
            'Color', 'r', 'FontSize', 8, 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
    end
end

% 绘制门限线
yline(min_peak_height, 'g--', 'LineWidth', 1.5, 'Label', sprintf('检测门限 (%d%%)', min_peak_height));

title(sprintf('全局归一化相关系数 (共发现 %d 个 Burst)', length(pks)));
xlabel('复采样点索引 (Complex Sample Index)');
ylabel('归一化相关系数 Match (%)');
ylim([0, 100]);
grid on; 
xlim([sample_axis(1), sample_axis(end)]);

fprintf('\n================================================================\n');
fprintf('处理完毕！请查看弹出的【全局连续突发同步图】。\n');
fprintf('================================================================\n');