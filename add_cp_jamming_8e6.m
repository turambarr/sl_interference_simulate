% 对 8e6.iq 加 SSS CP干扰的独立脚本
clear; clc; close all;

inFile = '8e6.iq';
outFile = 'cp_interference_0dB.iq';

fs_source = 409.6e6;       
fs_target = 60e6;          

% 加入解调时等效的中心频率与补偿
center_nominal_hz = 63.5e6; 
cfo_hz = -367188; 
freq_shift_hz = center_nominal_hz + cfo_hz; % 实频：63.132812 MHz

JSR_dB = 0;               

N = 1024;    
L = 48;      
W = 4;       
delay_L = 0;
disturb_L = L + delay_L * 2;

%% 构建 60MHz 干扰局部
disturb_list = [];
for i = 1 : 2 : disturb_L
    disturb_list = cat(2, disturb_list, exp(1j * 2 * pi * i / disturb_L) * sinc((i - disturb_L / 2) / W));
    disturb_list = cat(2, disturb_list, exp(1j * 2 * pi * (i + 1) / disturb_L) * sinc((i - disturb_L / 2) / W));
end

jam_sym_60m = zeros(N + L, 1);
for j = 0 : L - 1
    if mod(j, 2) == 0
        jam_sym_60m(1 + j) = disturb_list(j + 1);
    else
        jam_sym_60m(1 + j) = -disturb_list(j + 1);
    end
    jam_sym_60m(1 + N + j) = disturb_list(j + 1); 
end

%% 上采样到 409.6MHz 并上变频
[P, Q] = rat(fs_source / fs_target);
jam_sym_409m = resample(jam_sym_60m, P, Q);

% 增加与解调端对齐的低通滤波器（35MHz带宽截断），使干扰信号的频谱包络与真实波形完全对齐，避免带外泄露使其看起来“不像”
Wn = 35e6 / (fs_source/2);
lpf_order = 30;
b_lpf = fir1(lpf_order, Wn);
jam_sym_409m = filtfilt(b_lpf, 1, jam_sym_409m);

t_jam = (0 : length(jam_sym_409m) - 1).' / fs_source;
jam_sym_rf_base = jam_sym_409m .* exp(1j * 2 * pi * freq_shift_hz * t_jam);
L_jam_rf = length(jam_sym_rf_base);

%% 生成 PSS 模板 (保持和解调同步阶段一致的长序列，且载频强制在 63MHz 以对齐它的原装匹配滤波器)
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

pss_base = [-m_syms, m_syms, m_syms, m_syms, m_syms, m_syms, m_syms, m_syms];
N_cp_syms = 48;
cp_syms = -pss_base(end - N_cp_syms + 1 : end);
pss_base_with_cp = [cp_syms, pss_base];

pss_up = resample(pss_base_with_cp, P, Q);
t_local = (0:length(pss_up)-1) / fs_source;

% 同步寻找原代码是写死 63e6 作为本地 fc 的
fc_sync = 63e6;
pss_local = pss_up .* exp(1j * 2 * pi * fc_sync * t_local);
pss_local = pss_local / max(abs(pss_local)); 

%% 互相关寻找所有的 PSS 同步峰 (严格复刻 `sync_then_sss_demod.m` 中的归一化匹配)
fprintf('Reading %s...\n', inFile);
fid = fopen(inFile, 'rb');
raw = fread(fid, inf, 'int16');
fclose(fid);
x_iq = raw(1:2:end) + 1j * raw(2:2:end);
x_raw = double(x_iq);
x_sync = x_raw - mean(x_raw);
x_sync = x_sync / (mean(abs(x_sync)) + eps);

fprintf('Running exact normalized matched filter (from sync script)... \n');
h_matched = fliplr(conj(pss_local));
corr_out = filter(h_matched, 1, x_sync);
corr_mag_sq = abs(corr_out).^2;

E_local = sum(abs(pss_local).^2);
win_ones = ones(1, length(pss_local));
E_rx_moving = filter(win_ones, 1, abs(x_sync).^2);
E_rx_moving(E_rx_moving < 1e-10) = 1e-10;

corr_norm_pct = (corr_mag_sq ./ (E_local .* E_rx_moving)) * 100;

min_peak_height = 50; % 解调代码中往往是大于这个值
min_peak_dist = length(pss_local);
[~, locs_end] = findpeaks(corr_norm_pct, 'MinPeakHeight', min_peak_height, 'MinPeakDistance', min_peak_dist);

L_pss_409m = length(pss_local);
start_pos_list = locs_end - (L_pss_409m - 1);
fprintf('Find %d PSS peaks using normalized matched filter.\n', length(start_pos_list));

%% SSS CP 干扰注入
% 修改：由于 PSS 包含 1024 长度的底和 48 长度的 CP，总长为 1072
% 因此 SSS 的起始点相对于 PSS 起始点的偏移应当为 1072 而不是 1024！
offset_409m = round(1072 * (fs_source / fs_target));
x_out = x_raw;

for i = 1:length(start_pos_list)
    p_start = start_pos_list(i);
    if p_start < 1 || p_start + offset_409m + L_jam_rf > length(x_raw)
        continue;
    end
    
    jam_inject_range = p_start + offset_409m : p_start + offset_409m + L_jam_rf - 1;
    
    % 为了使得干扰信号具有全局连续的射频相位，我们在这一步直接根据局部的绝对时间 t_abs 生成射频载波
    t_abs = (jam_inject_range - 1).' / fs_source;
    jam_sym_rf_local = jam_sym_409m .* exp(1j * 2 * pi * freq_shift_hz * t_abs);
    
    local_sig = x_raw(jam_inject_range);
    sig_pwr = sum(abs(local_sig).^2) / L_jam_rf;
    jam_pwr = sum(abs(jam_sym_rf_local).^2) / L_jam_rf;
    
    % Jammer Power = Signal Power * 10^(JSR / 10)
    target_jam_pwr = sig_pwr * (10^(JSR_dB / 10));
    scale = sqrt(target_jam_pwr / jam_pwr);
    
    x_out(jam_inject_range) = x_out(jam_inject_range) + jam_sym_rf_local * scale;
    fprintf('  Peak %d injected (JSR=%.1fdB).\n', i, JSR_dB);
end

%% CP 相关性压制效果对比展示 (取第一个可用峰位进行前后对比)
if ~isempty(start_pos_list)
    fprintf('生成并显示第一个同步峰处的时域结构与 CP 压制效果对比图...\n');
    p_start = start_pos_list(1);
    
    % 扩大截取范围：包含 PSS符号、SSS符号以及前后的长裕量
    margin_409m = round(500 * (fs_source / fs_target));
    test_start = p_start - margin_409m;
    test_end = p_start + offset_409m + L_jam_rf + margin_409m;
    
    if test_start > 0 && test_end <= length(x_raw)
        test_range = test_start : test_end;
        t_test = (test_range - 1).' / fs_source;
        
        sig_clean_rf = x_raw(test_range);
        sig_jammed_rf = x_out(test_range);
        
        % ============= 图 1：时域幅度结构与加扰位置标注 =============
        figure('Name', 'Time Domain & Positions', 'Position', [100 500 900 350]);
        idx_rel = test_range - p_start; % 以 PSS 起点为 0 坐标
        
        plot(idx_rel, abs(sig_clean_rf), 'b', 'LineWidth', 1.5); hold on;
        plot(idx_rel, abs(sig_jammed_rf), 'r--', 'LineWidth', 1.5);
        
        y_lims = ylim; y_max = y_lims(2);
        
        % 标注 PSS 区间 (绿色半透明)
        patch([0, offset_409m, offset_409m, 0], [0, 0, y_max, y_max], 'g', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
        text(offset_409m/2, y_max*0.9, 'PSS 同步符号', 'HorizontalAlignment', 'center', 'Color', [0 0.5 0], 'FontWeight', 'bold');
        
        % 标注 加密/干扰 SSS 区间 (红色半透明)
        jam_end_idx = offset_409m + L_jam_rf - 1;
        patch([offset_409m, jam_end_idx, jam_end_idx, offset_409m], [0, 0, y_max, y_max], 'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
        text(offset_409m + L_jam_rf/2, y_max*0.9, 'SSS 被加干扰部分', 'HorizontalAlignment', 'center', 'Color', 'r', 'FontWeight', 'bold');
        
        legend('原始信号幅度', '注入干扰后幅度', 'Location', 'northeast');
        title('时域射频包络展示：同步 PSS 标识与干扰 SSS 注入段');
        xlabel('相对于 PSS 符头的物理采样偏移点 (409.6MHz)');
        ylabel('绝对包络强度 Abs(x)');
        grid on;
        
        % ============= 图 2：全景基带 CP 相关峰受损对比 =============
        % 1. 下变频
        sig_clean_ddc = sig_clean_rf .* exp(-1j * 2 * pi * freq_shift_hz * t_test);
        sig_jammed_ddc = sig_jammed_rf .* exp(-1j * 2 * pi * freq_shift_hz * t_test);
        
        % 2. 低通滤波
        sig_clean_lpf = filtfilt(b_lpf, 1, sig_clean_ddc);
        sig_jammed_lpf = filtfilt(b_lpf, 1, sig_jammed_ddc);
        
        % 3. 降采样恢复到 60MHz基带
        sig_clean_bb = resample(sig_clean_lpf, Q, P);
        sig_jammed_bb = resample(sig_jammed_lpf, Q, P);
        
        % 4. 计算 CP 滑动自相关 (N=1024, L=48)
        corr_len = length(sig_clean_bb) - N - L;
        cp_corr_clean = zeros(corr_len, 1);
        cp_corr_jammed = zeros(corr_len, 1);
        
        for k = 1:corr_len
            window_idx = k : k + L - 1;
            % Clean
            R_c = sum(sig_clean_bb(window_idx) .* conj(sig_clean_bb(window_idx + N)));
            P_c = sum(abs(sig_clean_bb(window_idx)).^2 + abs(sig_clean_bb(window_idx + N)).^2);
            cp_corr_clean(k) = 2 * abs(R_c) / (P_c + eps);
            % Jammed
            R_j = sum(sig_jammed_bb(window_idx) .* conj(sig_jammed_bb(window_idx + N)));
            P_j = sum(abs(sig_jammed_bb(window_idx)).^2 + abs(sig_jammed_bb(window_idx + N)).^2);
            cp_corr_jammed(k) = 2 * abs(R_j) / (P_j + eps);
        end
        
        figure('Name', 'CP Suppression Target Area', 'Position', [100 50 900 350]);
        % 将 X 坐标对齐到 PSS 起点(基带单位点)
        bb_offset = round(margin_409m * fs_target / fs_source); 
        bb_idx = (1:corr_len) - bb_offset; 
        
        plot(bb_idx, cp_corr_clean, 'b', 'LineWidth', 1.5); hold on;
        plot(bb_idx, cp_corr_jammed, 'r', 'LineWidth', 1.5);
        
        % 在相关图上也盖上位置提示块
        y_lims2 = ylim; y_max2 = y_lims2(2);
        bb_pss_len = round(offset_409m * fs_target / fs_source);
        bb_jam_len = round(L_jam_rf * fs_target / fs_source);
        
        patch([0, bb_pss_len, bb_pss_len, 0], [0, 0, y_max2, y_max2], 'g', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
        patch([bb_pss_len, bb_pss_len+bb_jam_len, bb_pss_len+bb_jam_len, bb_pss_len], [0, 0, y_max2, y_max2], 'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
        
        legend('原始 CP 滑动相关结果', '加压制干扰后 CP 结果', 'Location', 'northeast');
        title('CP 压制效果全局透视对比：被定位打击的 SSS 相关峰');
        xlabel('相对于 PSS 起始点的基带相对位置 (60MHz Baseband Samples)');
        ylabel('归一化 CP 相关强度 (2|R|/P)');
        grid on;
        drawnow;
    end
end

%% 输出 IQ
out_data = zeros(2 * length(x_out), 1);
out_data(1:2:end) = real(x_out);
out_data(2:2:end) = imag(x_out);
out_data = min(max(out_data, -32768), 32767); % 限制在 int16 范围防溢出翻转
out_data_int16 = int16(round(out_data));

f_out = fopen(outFile, 'wb');
fwrite(f_out, out_data_int16, 'int16');
fclose(f_out);
fprintf('Done! Saved %s\n', outFile);
