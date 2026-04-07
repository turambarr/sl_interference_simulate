% sss_demodulation.m
% SSS 单点解调脚本 (OFDM + 4QAM)
% 基于 GARDNER 模块估计出的 SRO 和 CFO 进行预修正

clear; clc; close all;

%% 1. 参数设置
inFile = 'sigtest1.iq';

% ===== 手动填写的 3 个关键参数 =====
read_start_sample = 14471; % [参数1] 文件中开始读取的点（原始采样点）
read_length = 6992*3+1000;                  % [参数2] 读取长度（原始采样点个数）
sss_decode_start_idx = 1024+48;           % [参数3] 从重采样后序列 x_sro 的第几个点开始作为 SSS 解调基准，也就是说SSS从此处开始

% 原始采样率
fs_source = 409.6e6;
% 目标符号/采样率 (OFDM基带)
fs_target = 60e6;

% --- 之前估计的同步参数 ---
% 请根据 gardner_farrow_timing_recovery.m 的输出回填这里
% 示例值：
sro_ppm  = 0;       % 采样率偏差 (ppm)
cfo_hz   = -367188;    % 载波频偏 (Hz)

% SSS 长度: 1个 OFDM 符号 = 1024 点 (假设无 CP 或与 pss_test.m 一致)
N_fft = 1024;
target_offset = 10; % 指定你想要测试的偏移量
freq_shift_hz = 63.5e6; % 频谱归基带向左搬移 63MHz

% 解调模式：true=解调全部1024子载波；false=仅解调有效子载波
demod_all_subcarriers = true;

% 判决前整体旋转角（用于把落在坐标轴上的点转回4QAM象限中心）
decision_rotate_deg = 45;

% ===== BER 计算参数 =====
enable_ber_eval = true;
% 预设正确值（4组）
ber_ref_hex_list = {
    '80C46B5267759FD6505644FC1F0052124876B17F9FA9EA13D3C2A8E274700BE9F511FBCDF1B237B1E6E46FB1A80A762952FCA4FCDB45612E5C6641BC200043442D23B538172067919D89A3D640557B004A5507460C3DC60DEDE45BB18D5FC9D5E944683265F591D6585671A960FF43002E13B478182063925276497C15006A10178064119F89AE831FAABF442C33B9F812205B9272769C29EAFC2E57F1E270347E7388504C853A85EF97AA0E79A3FDD5CC122832DE4EA2E7191A815164CC8CE96EEF4BC01B139D2E5D75B280250380B83D208AC7BD8ACD3A263D2B583218CAA1EB37998E6EF65A2982FFCE5624BB8FCD6F5F4483E3EC5384ED75493BE9CD690B', ... % Ref-1
    '159DC2F4CEEF3ABCF5FCDDA97A55F474D1EC27EA3A038076B6940184EDE55283AF77A29BA7246E278C8DCA270150EC43F4A90DA9B2DFC748F9CCD7294555D6DD4B462F617E45CE373B1306BCD5FFE255D0FF5EDC596B9C5B8B8DF2271BFA93BF83DDC164CFAF37BCF1FCE703C5AAD65548762DE17145C634F4ECD3E97F55C0757E15CD773A1308167A002ADD496623A17445F234E4EC394380A948FEA784E56DE8E611F5D91F601F8A3E0058E306ABBF99744164B8D8048E737017F7CD991983C88AD29572763B48FBEF24154F5615216B45109E2B109B604C6B42F164719007826E3318C8ACF04314AA98FC4D221A9BCAFADD168689F61D8BEFD362839AC352', ... % Ref-2
    '7F3B94AD988A6029AFA9BB03E0FFADEDB7894E80605615EC2C3D571D8B8FF4160AEE04320E4DC84E191B904E57F589D6AD035B0324BA9ED1A399BE43DFFFBCBBD2DC4AC7E8DF986E62765C29BFAA84FFB5AAF8B9F3C239F2121BA44E72A0362A16BB97CD9A0A6E29A7A98E569F00BCFFD1EC4B87E7DF9C6DAD89B683EAFF95EFE87F9BEE6076517CE05540BBD3CC4607EDDFA46D8D8963D61503D1A80E1D8FCB818C77AFB37AC57A106855F1865C022A33EDD7CD21B15D18E6E57EAE9B3373169110B43FE4EC62D1A28A4D7FDAFC7F47C2DF7538427532C5D9C2D4A7CDE7355E14C866719109A5D67D0031A9DB44703290A0BB7C1C13AC7B128AB6C4163096F4', ... % Ref-3
    'EA623D0B3110C5430A03225685AA0B8B2E13D815C5FC7F89496BFE7B121AAD7C50885D6458DB91D8737235D8FEAF13BC0B56F2564D2038B7063328D6BAAA2922B4B9D09E81BA31C8C4ECF9432A001DAA2F00A123A69463A474720DD8E4056C407C223E9B3050C8430E0318FC3A5529AAB789D21E8EBA39CB0B132C1680AA3F8A81EA3288C5ECF7E985FFD522B699DC5E8BBA0DCB1B13C6BC7F56B701587B1A921719EE0A26E09FE075C1FFA71CF95440668BBE9B4727FB718C8FE8083266E67C37752D6A8D89C4B70410DBEAB0A9EADE94BAEF61D4EF649FB394BD0E9B8E6FF87D91CCE737530FBCEB556703B2DDE564350522E9797609E274102C9D7C643CAD'  ... % Ref-4
};
unstable_bit_positions = [];      % 1-based bit索引
unstable_nibble_positions = [1, 2, 3, 4, 257, 510, 511, 512];   % 重点标记的零子载波

%% 2. 读取原始数据
fprintf('Loading file: %s from %d, len %d...\n', inFile, read_start_sample, read_length);
[x_raw, ~] = iq_read_int16_le(inFile, read_start_sample, read_length);
x_raw = double(x_raw);

% 归一化
x_raw = x_raw - mean(x_raw);
x_raw = x_raw / mean(abs(x_raw));

%% 3. 频谱下变频与 CFO 修正 (分步独立处理)
t_vec = (0:length(x_raw)-1) / fs_source;

% 第一步：在 409.6MHz 原始采样率下，向左搬移 63MHz (DDC 归基带)
fprintf('Shifting spectrum left by %.1f MHz at 409.6MHz...\n', freq_shift_hz/1e6);
x_shifted = x_raw .* exp(-1j * 2 * pi * freq_shift_hz * t_vec).'; 

% 第二步：单独进行 CFO 频偏补偿
fprintf('Applying separately CFO Correction: %.2f Hz...\n', cfo_hz);
x_cfo = x_shifted .* exp(-1j * 2 * pi * cfo_hz * t_vec).'; 

%% 4. SRO 修正 (采用 Farrow 分数阶延迟插值消除边界效应)
fprintf('Applying SRO Correction (Farrow Interpolation): %.2f ppm...\n', sro_ppm);
fs_eff = fs_source * (1 + sro_ppm/1e6);

% --- 新增：使用零相移(Zero-Phase)低通滤波器抗混叠 ---
% 目标下降到 60MHz 采样率，放宽截止频率为 35MHz
fprintf('Applying Zero-Phase Anti-aliasing LPF (35MHz)...\n');
Wn = 35e6 / (fs_source / 2);
lpf_order = 30;
b_lpf = fir1(lpf_order, Wn);
% 消除群时延(相位偏移)，过滤 CFO 修正后的信号
x_cfo_filtered = filtfilt(b_lpf, 1, x_cfo);

% 使用三次 Lagrange (Cubic Lagrange) Farrow 结构进行重采样
T_in = 1 / fs_eff;
T_out = 1 / fs_target;
% 避免两端由于缺少足量插值基准点而越界，舍弃最后微小的尾部
t_out = 0 : T_out : (length(x_cfo_filtered)-3)*T_in; 

% 映射到输入序列的虚拟索引 (1-based)
idx_frac = t_out / T_in + 1; 
idx_base = floor(idx_frac);
mu = idx_frac - idx_base;

% 确保 base 索引不会越界 (Farrow Cubic 需要 base-1 到 base+2 共4个点)
valid_mask = (idx_base >= 2) & (idx_base <= length(x_cfo_filtered)-2);
idx_base = idx_base(valid_mask);
mu = mu(valid_mask);

% Farrow 立方插值滤波器系数 (Cubic Lagrange)
h0 = -(mu - 1) .* (mu - 2) .* mu / 6;
h1 =  (mu - 1) .* (mu + 1) .* (mu - 2) / 2;
h2 = -(mu + 1) .* mu .* (mu - 2) / 2;
h3 =  (mu + 1) .* (mu - 1) .* mu / 6;

% 卷积求和（无边界失真效应）
x_sro = h0 .* x_cfo_filtered(idx_base - 1) + ...
        h1 .* x_cfo_filtered(idx_base) + ...
        h2 .* x_cfo_filtered(idx_base + 1) + ...
        h3 .* x_cfo_filtered(idx_base + 2);
x_sro = x_sro(:); % 转为列向量

%% 5. 定位 SSS 符号与时域分析
fprintf('Analyzing Time Domain Signal...\n');

% 这里使用你手动填写的解调起始点（x_sro 坐标）
base_start_idx = sss_decode_start_idx;

fig_const = figure('Position', [300, 300, 600, 600], 'Name', 'SSS Demodulation Result (Single Offset)');

% 单点偏移测试（不做遍历）
offset = target_offset;
sss_start_idx_60 = base_start_idx + offset;
        
    % 防止越界
if sss_start_idx_60 < 1 || sss_start_idx_60 + N_fft > length(x_sro)
    error('Offset %d 导致下标越界，无法计算。', offset);
end

x_sss_time = x_sro(sss_start_idx_60 : sss_start_idx_60 + N_fft - 1);
    
%% 6. OFDM 解调 (FFT) 与 小数定时偏差/相位纠正
x_sss_freq = fft(x_sss_time, N_fft) / sqrt(N_fft);  
    
% 在这里我们先提取出有效子载波用来做相位拟合 (去除空载波对拟合的影响)
% --- 修正：采用幅度判断法自动鉴别有效载波边缘，而不是写死 824 ---
% 计算频点功率，找出平均功率，滤除深陷低于某个阈值的零点或衰落点
pwr_bins = abs(x_sss_freq).^2;
mean_pwr = mean(pwr_bins);
threshold_pwr = mean_pwr * 0.1; % 设定阈值（例如低于均值 1/10 的当成是噪声空载波）
    
% 找到功率大于阈值的点作为有效子载波 (不把处于基带 0Hz 附近或边缘的无效点算进来)
valid_mask_power = pwr_bins > threshold_pwr;
    
% 进一步滤除最最中心的直流点(DC, 也就是 MATLAB index 1 附近的高风险带)
% 我们可以用一段保护带来防住零载波
dc_guard = 3; % 左右丢弃几个子载波当做保护段免受本振泄露干扰
valid_mask_power(1:dc_guard) = false;
valid_mask_power(N_fft-dc_guard+1:N_fft) = false;
    
% 获取有效载波的索引 (逻辑 1-1024)
valid_idxs = find(valid_mask_power);
    
% 构建基于物理频率的相对频率轴(-512 to 511)
% N_fft/2+1 到 N_fft 是负频点
rel_freq = zeros(N_fft, 1);
rel_freq(1:N_fft/2) = (0:N_fft/2 - 1).'; 
rel_freq(N_fft/2+1:N_fft) = (-N_fft/2:-1).';
    
freq_indices = rel_freq(valid_idxs);
syms_valid = x_sss_freq(valid_idxs);

% --- 严重 BUG 修复：必须对频率进行从小到大排序 ---
% 原始的 FFT 抽取顺序是先从 0 到 +511，再从 -512 到 -1
% 如果不排序（也就是从最高正频率突然跳变到最低负频率），unwrap 解卷绕时会产生巨大的 2*pi 错位！
% 这将彻底毁掉后面的 polyfit 线性拟合，导致所有的相位补偿全错
[freq_indices, sort_idx] = sort(freq_indices);
syms_valid = syms_valid(sort_idx);
    
% 注意：由于我们对 valid_idxs 对应的 syms_valid 进行了重排，后续根据 valid_idxs 原顺序截取 hard bits 会不对应
% 所以最后用来打印 Bits 的 valid_idxs 也必须一起重排，才能正确找到对应的子载波
valid_idxs = valid_idxs(sort_idx);

% ===== 使用 sweep 中验证过的相位估计逻辑：分段 unwrap + 分段 polyfit =====
syms_pow4 = syms_valid.^4;
df = diff(freq_indices);

% 以 df>1 的位置作为频率断点（例如跨越 DC 缺口）
seg_start = [1; find(df > 1) + 1];
seg_end   = [find(df > 1); length(freq_indices)];

slope_list = [];
w_list = [];
for k = 1:length(seg_start)
    idx_seg = seg_start(k):seg_end(k);
    if numel(idx_seg) < 8
        continue;
    end

    f_seg = freq_indices(idx_seg);
    ph_seg = unwrap(angle(syms_pow4(idx_seg)));
    p_seg = polyfit(f_seg, ph_seg, 1);

    slope_list(end+1,1) = p_seg(1); %#ok<SAGROW>
    w_list(end+1,1) = numel(idx_seg); %#ok<SAGROW>
end

if isempty(slope_list)
    error('Offset %d 连续有效载波不足，无法完成分段相位拟合。', offset);
end

% 用段长加权平均估计 4 次方域斜率
slope4 = sum(slope_list .* w_list) / sum(w_list);

% 已知 slope4 后，用复均值估计 4 次方域常量相位项（对噪声更稳）
pow4_detrended = syms_pow4 .* exp(-1j * slope4 * freq_indices);
phase0_4 = angle(mean(pow4_detrended));

% 计算补偿用相位：除以 4 恢复出原始偏移
full_freq_indices = [(0:N_fft/2-1), (-N_fft/2:-1)].';
phase_correction = (slope4 * full_freq_indices + phase0_4) / 4;
    
% 5. 应用补偿
x_sss_freq_corr = x_sss_freq .* exp(-1j * phase_correction);

% 提取所有子载波用于最终显示
syms_all = x_sss_freq_corr(:);
    
%% 7. 绘图与解调判决
clf(fig_const); % 清除上一帧的图
    
% --- 绘制星座图 ---
subplot(1, 2, 1);
syms_all_plot = syms_all .* exp(1j * decision_rotate_deg * pi/180);
plot(real(syms_all_plot), imag(syms_all_plot), 'b.', 'MarkerSize', 8);
hold on; grid on; axis square;
    
th = 0:0.01:2*pi;
plot(cos(th), sin(th), 'k:', 'LineWidth', 0.5); % 单位圆
    
title(sprintf('SSS Constellation (Rotated %d°)\nStart Offset: %d 点\nResidual STO slope: %.4f', decision_rotate_deg, offset, slope4/4));
xlabel('I'); ylabel('Q');
xlim([-2 2]); ylim([-2 2]);
    
% --- 提取有效载波并进行 4QAM (QPSK) 硬判决 ---
subplot(1, 2, 2);
% 从已做相位补偿的全频段数据中提取待解调子载波
if demod_all_subcarriers
    syms_payload = x_sss_freq_corr(:);  % 全部 1024 子载波
    payload_label = sprintf('All carriers: %d', N_fft);
else
    syms_payload = x_sss_freq_corr(valid_idxs);
    payload_label = sprintf('Valid carriers: %d', length(valid_idxs));
end
    
% 执行 4QAM/QPSK 极性硬判决
% 先整体旋转 decision_rotate_deg，再考虑 QPSK 的 90/180/270 度歧义
rot_angles_deg = decision_rotate_deg + [0, 90, 180, 270];
hex_candidates = cell(1, 4);
demod_bits_candidates = cell(1, 4);

for r = 1:4
    % 旋转补偿候选：decision_rotate_deg + {0,90,180,270} 度
    syms_payload_rot = syms_payload .* exp(1j * decision_rotate_deg * pi/180) .* exp(-1j * (r-1) * pi/2);

    bits_I = real(syms_payload_rot) < 0;
    bits_Q = imag(syms_payload_rot) < 0;

    % 将 I 和 Q 合并成一溜长比特数据流 (I,Q交织存放)
    demod_bits = zeros(length(syms_payload_rot)*2, 1);
    demod_bits(1:2:end) = bits_I;
    demod_bits(2:2:end) = bits_Q;
    demod_bits_candidates{r} = demod_bits;

    % --- 将解调出的所有比特转化为完整的 Hex 字符串 ---
    full_hex_str = '';
    for i_hex = 1:4:length(demod_bits)
        if i_hex+3 > length(demod_bits)
            % 如果最后有不足 4 比特的零碎数据，补齐处理
            chunk = demod_bits(i_hex:end);
            val = 0;
            for b = 1:length(chunk)
                val = val + chunk(b) * 2^(length(chunk)-b);
            end
            full_hex_str = [full_hex_str, dec2hex(val)]; %#ok<AGROW>
            break;
        end
        chunk = demod_bits(i_hex:i_hex+3);
        val = chunk(1)*8 + chunk(2)*4 + chunk(3)*2 + chunk(4)*1;
        full_hex_str = [full_hex_str, dec2hex(val)]; %#ok<AGROW>
    end

    hex_candidates{r} = full_hex_str;
end
    
% 在副图仅打印前 32 个 hex 字符作预览即可
axis off;
title('Demodulated Hard Bits');
hex_snip_0   = hex_candidates{1}(1:min(20, length(hex_candidates{1})));
hex_snip_90  = hex_candidates{2}(1:min(20, length(hex_candidates{2})));
hex_snip_180 = hex_candidates{3}(1:min(20, length(hex_candidates{3})));
hex_snip_270 = hex_candidates{4}(1:min(20, length(hex_candidates{4})));
text(0.05, 0.78, sprintf('Offset = %d\n0°  : %s...\n90° : %s...\n180°: %s...\n270°: %s...', ...
    offset, hex_snip_0, hex_snip_90, hex_snip_180, hex_snip_270), 'FontSize', 11, 'Interpreter', 'none');
text(0.05, 0.34, sprintf('Total bits: %d\n%s * 2 bit/sym', length(demod_bits_candidates{1}), payload_label), 'FontSize', 12);
    
% 无论 offset 是多少，都直接打印最终结果
fprintf('\n>>> [Offset = %d] SSS Demodulation Completed! <<<\n', offset);
fprintf('Extracted %d bits per candidate (%s).\n', length(demod_bits_candidates{1}), payload_label);
fprintf('====== HEX CANDIDATES (QPSK 4-way rotation ambiguity) ======\n');
for r = 1:4
    fprintf('[%3d deg] %s\n', rot_angles_deg(r), hex_candidates{r});
end
fprintf('============================================================\n');

%% 8. BER 评估（与4组预设正确值对比，剔除0子载波不稳位）
if enable_ber_eval
    ref_filled = cellfun(@(s) ~isempty(strtrim(s)), ber_ref_hex_list);
    if any(ref_filled)
        nCand = numel(hex_candidates);
        nRef = numel(ber_ref_hex_list);
        ber_mat = nan(nCand, nRef);

        fprintf('\n================ BER EVALUATION ================\n');
        fprintf('Unstable drop: bit_pos=%d, nibble_pos=%d\n', ...
            numel(unstable_bit_positions), numel(unstable_nibble_positions));

        for i = 1:nCand
            bits_est = demod_bits_candidates{i};
            for j = 1:nRef
                ref_hex = strtrim(ber_ref_hex_list{j});
                if isempty(ref_hex)
                    continue;
                end
                
                % 本地函数将hex转换为bits
                bits_ref = hex_to_bits_local(ref_hex);
                nUse = min(length(bits_est), length(bits_ref));
                if nUse <= 0
                    continue;
                end

                keep_mask = true(nUse,1);

                % 剔除不稳 bit 位
                ub = unstable_bit_positions(:);
                ub = ub(ub >= 1 & ub <= nUse);
                keep_mask(ub) = false;

                % 剔除不稳 nibble 位（每个hex字符对应4bit）
                un = unstable_nibble_positions(:);
                un = un(un >= 1);
                for k = 1:numel(un)
                    b1 = 4*(un(k)-1) + 1;
                    b2 = 4*un(k);
                    if b1 > nUse
                        continue;
                    end
                    b2 = min(b2, nUse);
                    keep_mask(b1:b2) = false;
                end

                bits_e = bits_est(1:nUse);
                bits_r = bits_ref(1:nUse);
                nKeep = nnz(keep_mask);
                if nKeep <= 0
                    continue;
                end

                err = nnz(bits_e(keep_mask) ~= bits_r(keep_mask));
                ber_mat(i,j) = err / nKeep;

                fprintf('Cand[%3d deg] vs Ref[%d]: BER=%.6g (err=%d/%d)\n', ...
                    rot_angles_deg(i), j, ber_mat(i,j), err, nKeep);
            end
        end

        % 打印每个候选的最优匹配参考
        for i = 1:nCand
            row = ber_mat(i,:);
            if all(isnan(row))
                continue;
            end
            [v, idx] = min(row);
            fprintf('Best for Cand[%3d deg] -> Ref[%d], BER=%.6g\n', rot_angles_deg(i), idx, v);
        end
        fprintf('================================================\n');
    end
end

fprintf('\n解调流程结束。\n');

%% ===== 辅助函数 =====
function bits_out = hex_to_bits_local(hex_str)
    % 辅助函数：将一行HEX字符串解析回0/1列向量
    len_hex = length(hex_str);
    bits_out = zeros(len_hex * 4, 1);
    
    for idx_char = 1:len_hex
        val = hex2dec(hex_str(idx_char));
        bin_str = dec2bin(val, 4);
        bidx = (idx_char-1)*4 + 1;
        bits_out(bidx)   = bin_str(1) - '0';
        bits_out(bidx+1) = bin_str(2) - '0';
        bits_out(bidx+2) = bin_str(3) - '0';
        bits_out(bidx+3) = bin_str(4) - '0';
    end
end

