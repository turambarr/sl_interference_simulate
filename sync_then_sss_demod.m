% sync_then_sss_demod.m
% 同步-解调一体脚本：
% 1) 先复刻 frame_sync_routine 的归一化匹配找峰，自动得到 burst 起始点
% 2) 将该起始点作为 SSS 解调的 read_start_sample
% 3) 执行 SSS 单点解调（支持可选 CFO 补偿）

clear; clc; close all;

%% ===================== 0. 全局参数 =====================
% inFile = 'target_signal_multiframe_jammed_3dB.dat';
% inFile = 'sigtest2.iq';
% inFile = 'target_signal_jammed_-3dB.dat';
inFile = "target_signal_multiframe_IF_jammed_-14dB.iq";

% 输入文件读取设置：
% input_format = 'auto' | 'iq' | 'dat'
% - auto: 按扩展名自动判断（.iq/.dat）
% - iq/dat: 强制按指定格式读取
input_format = 'auto';
dat_header_bytes = 0; % 若 .dat 有文件头可设置为非0

% --- 同步扫描参数（复刻 frame_sync_routine） ---
sync_scan_start = 0;      % 0-based
sync_scan_len   = 20e6;   % 最多读取多少点用于找峰
min_peak_height = 50;     % 峰值门限（%）
peak_select_mode = 'first'; % 'first' or 'strongest'

% 时域全图显示参数
plot_full_time_domain = true;   % true: 绘制同步扫描区完整时域幅度图
time_plot_max_points = 2e6;     % 显示抽点上限（仅影响绘图，不影响处理）

% --- SSS 解调参数 ---
read_length = 6992*3 + 1000;   % 从 read_start_sample 开始读取长度（409.6MHz采样率）
sss_decode_start_idx = 1024 + 48;
target_offset = 0;
enable_isolated_sss_demod_stream = true; % 同步流与SSS解调流分离，避免PSS/周边干扰进入降采样链路
isolated_sss_pre_guard_60 = 0;           % 独立SSS流在SSS前保留的60MHz样点；0表示不读PSS尾部
isolated_sss_post_guard_60 = 2048;       % 独立SSS流在SSS后保留的60MHz样点，用于滤波/重采样余量

fs_source = 409.6e6;
fs_target = 60e6;
center_nominal_hz = 63.5e6; % 名义中心频率（用于显示/对比）

% 频率补偿模式：
% 'manual'      -> 使用手填 center + cfo
% 'blind_pilot' -> 自动盲估“中心单载波频率”，一步补偿（cfo 置0）
freq_comp_mode = 'blind_pilot';

% 手动模式参数
freq_shift_hz_manual = 63.5e6;   % DDC 中心频率
cfo_hz_manual = -367188;         % CFO

% CFO 补偿参数（可选）
enable_cfo_comp = true;   % true: 启用固定 CFO 补偿；false: 跳过

% 盲估模式参数（复用你实验脚本的方法）
blind_target_bw_hz = 60e6;
blind_bw_tol_ratio = 0.35;
blind_occ_bg_win_hz = 2.0e6;
blind_occ_smooth_hz = 1.2e6;
blind_occ_thresh_db = 12.0; % [修改] 从5.5提高到12，确保能切出信号的高原台阶，避开整个底噪带
blind_min_component_bw_hz = 20e6;
blind_max_expected_cfo_hz = 6e6;
blind_pilot_bg_win_hz = 0.9e6;
blind_pilot_min_prom_db = 4.5;
blind_pilot_width_ref_hz = 0.25e6;
blind_est_pre_start_samples = round(5000 * 409.6 / 600);  % 从 PSS 起点前多少点开始盲估（409.6MHz采样率）
blind_est_len = round(40000 * 409.6 / 600);               % 盲估时域长度（409.6MHz采样率）

% 盲估后时域精修（你当前场景峰值很清晰，默认关闭，直接用频域峰值）
blind_refine_enable = false;
blind_refine_lpf_bw_hz = 0.6e6;   % 初始窄带低通带宽（围绕单音）
blind_refine_iters = 2;           % 迭代次数
blind_refine_fir_order = 256;     % 低通 FIR 阶数
blind_refine_max_delta_hz = 120e3; % 精修相对频域粗估的最大允许修正（超限则回退）

sro_ppm = 0;              % 若已知可回填

N_fft = 1024;
demod_all_subcarriers = true;
decision_rotate_deg = 45;

% SSS 解调前的干净上下文构造：
% raw: 在409.6MHz域先构造 [SSS; SSS; SSS]，再做DDC/低通/SRO。
% post_sro: 在SRO后、FFT前构造 [SSS; SSS; SSS]。
enable_raw_sss_repeat_context = true;
enable_sss_repeat_context = false;

% x_raw 归一化方式：
% 'full_mean'        -> 原始整段 read_length 去均值/按 mean(abs) 缩放
% 'sss_local_robust' -> 仅用 SSS 附近局部窗口估计 DC/幅度，降低区外干扰影响
raw_norm_mode = 'sss_local_robust';
raw_norm_guard_60 = 256; % SSS FFT 窗口前后额外纳入的 60MHz 样点数
enable_raw_margin_clip = true;
raw_margin_clip_factor = 1;    % 归一化后，余量区幅度超过该值则限幅
raw_margin_protect_guard_60 = 128; % SSS FFT 窗口前后不做限幅保护区（60MHz样点）

% 主解调低通方式：
% 'filtfilt'    -> 双向零相位滤波，可能把目标区后的强干扰反向带入目标区
% 'single_pass' -> 单向 FIR 滤波，并按 lpf_order/2 补偿固定群延迟，用于对照测试
main_lpf_mode = 'filtfilt';

% 有效载波选择：
% 'dynamic' -> 按当前 SSS 符号功率门限检测
% 'fixed'   -> 使用 fixed_valid_rel_freq_indices 或 fixed_valid_fft_bins
valid_carrier_mode = 'fixed';
valid_carrier_power_threshold_ratio = 0.1;
valid_carrier_dc_guard = 3;
fixed_valid_rel_freq_indices = [-512:-10, 7:511]; % 相对频率编号
fixed_valid_fft_bins = [];         % 例如 1:1024，MATLAB FFT bin 编号（1-based）

% BER 计算参数（与4个候选结果逐一比对）
enable_ber_eval = true;
enable_ref_aided_select = true;    % 使用预设 SSS 序列参与 offset 候选选择
enable_ref_aided_comp = true;      % 使用预设 SSS 序列直接参与相位斜率/公共相位补偿
enable_ref_guided_phase_fit = true; % 用预设序列去调制后，引导 unwrap + 平滑相位拟合
ref_guided_phase_fit_order = 2;     % 1=线性相位，2=允许轻微相位弯曲
ref_guided_fit_max_iters = 3;       % 残差剔除迭代次数
ref_guided_fit_outlier_sigma = 3.0; % 残差剔除强度
enable_ref_segment_anchor_fit = true; % 分段参考锚点拟合，抗局部受扰载波段
ref_anchor_seg_len = 48;              % 每段有效载波数
ref_anchor_min_coh = 0.1;            % 分段相位锚点最低可信度
ref_anchor_fit_order = 2;             % 锚点拟合阶数
ref_anchor_min_segments = 5;          % 最少可信锚点段数
enable_ref_aided_eq = false;       % 强制感较重的逐子载波相位均衡，默认关闭
ref_aided_eq_smooth_bins = 3;      % 越小越“贴参考”，1 会退化成几乎逐子载波相位强制对齐
ref_offset_search = -8:8;        % 在 target_offset 附近搜索的 60MHz 样点偏移
ref_accept_ber_thresh = 0.03;      % 超过该 BER 则判为无法可靠锁定，不硬改结果
% 预设正确值（4组）。可填你在 sss_group_summary_173.csv 的 voted_hex
ber_ref_hex_list = {
    '80C46B5267759FD6505644FC1F0052124876B17F9FA9EA13D3C2A8E274700BE9F511FBCDF1B237B1E6E46FB1A80A762952FCA4FCDB45612E5C6641BC200043442D23B538172067919D89A3D640557B004A5507460C3DC60DEDE45BB18D5FC9D5E944683265F591D6585671A960FF43002E13B478182063925276497C15006A10178064119F89AE831FAABF442C33B9F812205B9272769C29EAFC2E57F1E270347E7388504C853A85EF97AA0E79A3FDD5CC122832DE4EA2E7191A815164CC8CE96EEF4BC01B139D2E5D75B280250380B83D208AC7BD8ACD3A263D2B583218CAA1EB37998E6EF65A2982FFCE5624BB8FCD6F5F4483E3EC5384ED75493BE9CD690B', ... % Ref-1
    '159DC2F4CEEF3ABCF5FCDDA97A55F474D1EC27EA3A038076B6940184EDE55283AF77A29BA7246E278C8DCA270150EC43F4A90DA9B2DFC748F9CCD7294555D6DD4B462F617E45CE373B1306BCD5FFE255D0FF5EDC596B9C5B8B8DF2271BFA93BF83DDC164CFAF37BCF1FCE703C5AAD65548762DE17145C634F4ECD3E97F55C0757E15CD773A1308167A002ADD496623A17445F234E4EC394380A948FEA784E56DE8E611F5D91F601F8A3E0058E306ABBF99744164B8D8048E737017F7CD991983C88AD29572763B48FBEF24154F5615216B45109E2B109B604C6B42F164719007826E3318C8ACF04314AA98FC4D221A9BCAFADD168689F61D8BEFD362839AC352', ... % Ref-2
    '7F3B94AD988A6029AFA9BB03E0FFADEDB7894E80605615EC2C3D571D8B8FF4160AEE04320E4DC84E191B904E57F589D6AD035B0324BA9ED1A399BE43DFFFBCBBD2DC4AC7E8DF986E62765C29BFAA84FFB5AAF8B9F3C239F2121BA44E72A0362A16BB97CD9A0A6E29A7A98E569F00BCFFD1EC4B87E7DF9C6DAD89B683EAFF95EFE87F9BEE6076517CE05540BBD3CC4607EDDFA46D8D8963D61503D1A80E1D8FCB818C77AFB37AC57A106855F1865C022A33EDD7CD21B15D18E6E57EAE9B3373169110B43FE4EC62D1A28A4D7FDAFC7F47C2DF7538427532C5D9C2D4A7CDE7355E14C866719109A5D67D0031A9DB44703290A0BB7C1C13AC7B128AB6C4163096F4', ... % Ref-3
    'EA623D0B3110C5430A03225685AA0B8B2E13D815C5FC7F89496BFE7B121AAD7C50885D6458DB91D8737235D8FEAF13BC0B56F2564D2038B7063328D6BAAA2922B4B9D09E81BA31C8C4ECF9432A001DAA2F00A123A69463A474720DD8E4056C407C223E9B3050C8430E0318FC3A5529AAB789D21E8EBA39CB0B132C1680AA3F8A81EA3288C5ECF7E985FFD522B699DC5E8BBA0DCB1B13C6BC7F56B701587B1A921719EE0A26E09FE075C1FFA71CF95440668BBE9B4727FB718C8FE8083266E67C37752D6A8D89C4B70410DBEAB0A9EADE94BAEF61D4EF649FB394BD0E9B8E6FF87D91CCE737530FBCEB556703B2DDE564350522E9797609E274102C9D7C643CAD'  ... % Ref-4
};
% 不稳位剔除：支持两种方式（二选一或都填）
unstable_bit_positions = [];      % 1-based bit索引
unstable_nibble_positions = [1, 2, 3, 4, 257, 510, 511, 512];   % 来自 sss_group_summary_173.csv 的 marked 位

%% ===================== 1. 构造本地 PSS 模板 =====================
fprintf('Step 1: 生成本地 PSS 模板...\n');
fs_symbol = 60e6;
fc = 63e6;

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

[P, Q] = rat(fs_source / fs_symbol);
pss_up = resample(pss_base_with_cp, P, Q);
t_local = (0:length(pss_up)-1) / fs_source;
pss_local = pss_up .* exp(1j * 2 * pi * fc * t_local);
pss_local = pss_local / max(abs(pss_local));

%% ===================== 2. 全局同步找峰 =====================
fprintf('Step 2: 全局归一化匹配找峰...\n');

fInfo = dir(inFile);
if isempty(fInfo)
    error('找不到文件: %s', inFile);
end

dataBytes = fInfo.bytes;
if strcmpi(input_format, 'dat') || (strcmpi(input_format, 'auto') && endsWith(lower(inFile), '.dat'))
    dataBytes = max(0, fInfo.bytes - dat_header_bytes);
end
totalSamples = floor(dataBytes / 4);
if sync_scan_start < 0 || sync_scan_start >= totalSamples
    error('sync_scan_start 越界: %d / total=%d', sync_scan_start, totalSamples);
end

sync_scan_len = min(round(sync_scan_len), totalSamples - sync_scan_start);
[x_sync, meta_sync] = read_iq_auto_local(inFile, sync_scan_start, sync_scan_len, input_format, dat_header_bytes);
if meta_sync.numSamplesRead <= 0
    error('同步扫描阶段未读取到有效数据。');
end

x_sync = double(x_sync(:));
x_sync = x_sync - mean(x_sync);
x_sync = x_sync / (mean(abs(x_sync)) + eps);

% 频率估计/补偿参数在 PSS 相关确定 read_start_sample 后再设置。
% manual 模式直接使用固定参数；blind_pilot 模式从该起点之后读取一段信号估计。
switch lower(freq_comp_mode)
    case 'manual'
        fprintf('Freq mode: MANUAL, center=%.6f MHz, cfo=%.3f kHz\n', ...
            freq_shift_hz_manual/1e6, cfo_hz_manual/1e3);
    case 'blind_pilot'
        fprintf('Freq mode: BLIND_PILOT, will estimate after PSS sync start is selected.\n');
    otherwise
        error('未知 freq_comp_mode: %s', freq_comp_mode);
end

% 归一化匹配滤波（使用 fftfilt 大幅加速千万级点位与上千Tap滤波器的卷积）
h_matched = fliplr(conj(pss_local));
corr_out = fftfilt(h_matched(:), x_sync(:));
corr_mag_sq = abs(corr_out).^2;

E_local = sum(abs(pss_local).^2);
win_ones = ones(1, length(pss_local));
E_rx_moving = fftfilt(win_ones(:), abs(x_sync(:)).^2);
E_rx_moving(E_rx_moving < 1e-10) = 1e-10;

corr_norm_pct = (corr_mag_sq ./ (E_local .* E_rx_moving)) * 100;

min_peak_dist = length(pss_local);
[pks, locs] = findpeaks(corr_norm_pct, 'MinPeakHeight', min_peak_height, 'MinPeakDistance', min_peak_dist);

read_start_sample_list = [];
mode_names = {};
L_pss = length(pss_local);

if isempty(locs)
    figure('Name', 'Debug: Sync Peak Output', 'Position', [200,200,900,400]);
    plot(corr_norm_pct, 'b');
    yline(min_peak_height, 'r--', 'Threshold');
    title('Debug: 归一化匹配度度全景图 (未检测到峰值)');
    xlabel('Samples'); ylabel('Correlation (%)');
    grid on; axis tight;
    drawnow;
    
    fprintf('\n[警告] 未检测到大于 %d%% 的有效同步峰。已绘制扫描段度量图。\n', min_peak_height);
    fprintf('将执行两种降级解调方案：\n');
    read_start_sample_list(1) = sync_scan_start;
    mode_names{1} = '降级方案 1: 扫描起始点当作起点';
    
    [~, max_idx] = max(corr_norm_pct);
    fallback_pos2 = sync_scan_start + (max_idx - 1) - (L_pss - 1);
    read_start_sample_list(2) = max(0, fallback_pos2);
    mode_names{2} = '降级方案 2: 最大峰位置当作起点';
    
    fprintf('  1) 起始点=%d (%s)\n', read_start_sample_list(1), mode_names{1});
    fprintf('  2) 起始点=%d (%s)\n', read_start_sample_list(2), mode_names{2});
else
    end_pos_global = sync_scan_start + (locs - 1);
    start_pos_global = end_pos_global - (L_pss - 1);

    fprintf('  检测到 %d 个峰：\n', length(locs));
    for i = 1:length(locs)
        fprintf('    Peak %02d: %.2f%%, 末端=%d, 起始=%d\n', i, pks(i), end_pos_global(i), start_pos_global(i));
    end

    switch lower(peak_select_mode)
        case 'first'
            pick_idx = 1;
        case 'strongest'
            [~, pick_idx] = max(pks);
        otherwise
            error('未知 peak_select_mode: %s', peak_select_mode);
    end

    read_start_sample = start_pos_global(pick_idx);
    if read_start_sample < 0
        error('选中的起始点为负数: %d（请检查峰选择）', read_start_sample);
    end

    fprintf('\n>>> 选中峰 #%d, 匹配度 %.2f%%\n', pick_idx, pks(pick_idx));
    fprintf('>>> 自动设置 read_start_sample = %d\n', read_start_sample);
    
    read_start_sample_list(1) = read_start_sample;
    mode_names{1} = '正常同步峰起点';
end

% 完整时域图（同步扫描区）
if plot_full_time_domain
    sample_axis_td = sync_scan_start + (0:length(x_sync)-1);
    stride = max(1, ceil(length(sample_axis_td) / time_plot_max_points));
    idx_disp = 1:stride:length(sample_axis_td);

    figure('Position', [80, 80, 1400, 520], 'Name', 'Full Time Domain (Sync Scan Range)');
    plot(sample_axis_td(idx_disp), abs(x_sync(idx_disp)), 'Color', [0.55 0.55 0.55], 'LineWidth', 0.8);
    hold on; grid on;

    % 标注所有检测峰（起始/末端）
    if exist('start_pos_global', 'var')
        for ii = 1:length(start_pos_global)
            xline(start_pos_global(ii), 'b--', 'LineWidth', 0.8);
            xline(end_pos_global(ii), 'r--', 'LineWidth', 0.8);
        end
    end

    % 标注最终选中/降级测试起始点
    for r_idx = 1:length(read_start_sample_list)
        xline(read_start_sample_list(r_idx), 'g-', 'LineWidth', 1.8);
    end

    xlabel('复采样点索引 (Complex Sample Index)');
    ylabel('|x| (Normalized Amplitude)');
    title(sprintf('Full Time Domain (N=%d, stride=%d) | Blue=Start, Red=End, Green=Selected Start', length(x_sync), stride));
    xlim([sample_axis_td(1), sample_axis_td(end)]);
end

% 同步结果图
figure('Position', [80, 80, 1400, 520], 'Name', 'Sync Peaks (Normalized Correlation)');
sample_axis_sync = sync_scan_start + (0:length(corr_norm_pct)-1);
plot(sample_axis_sync, corr_norm_pct, 'b', 'LineWidth', 1.1); hold on; grid on;
plot(sync_scan_start + (locs - 1), pks, 'r*', 'MarkerSize', 8, 'LineWidth', 1.2);
yline(min_peak_height, 'g--', 'LineWidth', 1.2);
xlabel('复采样点索引 (Complex Sample Index)');
ylabel('归一化相关系数 Match (%)');
title(sprintf('Sync Peaks (detected=%d, starts to test=%d)', length(locs), length(read_start_sample_list)));

%% ===================== 3. SSS 解调（支持多模式/降级方案） =====================
for run_idx = 1:length(read_start_sample_list)
    read_start_sample = read_start_sample_list(run_idx);
    current_mode = mode_names{run_idx};
    
    fprintf('\n========================================================================\n');
    fprintf('Step 3: 使用 [%s] 起始点 (read_start_sample=%d) 进行 SSS 解调...\n', current_mode, read_start_sample);
    fprintf('========================================================================\n');
    
if enable_isolated_sss_demod_stream
    demod_stream_start_60 = max(1, sss_decode_start_idx - isolated_sss_pre_guard_60);
    demod_stream_end_60 = sss_decode_start_idx + N_fft - 1 + isolated_sss_post_guard_60;
    demod_stream_start_offset_raw = floor((demod_stream_start_60 - 1) * fs_source / fs_target);
    demod_stream_len_raw = ceil((demod_stream_end_60 - demod_stream_start_60 + 1) * fs_source / fs_target) + 8;
    demod_read_start_sample = read_start_sample + demod_stream_start_offset_raw;
    demod_read_length = min(round(demod_stream_len_raw), totalSamples - demod_read_start_sample);
    sss_decode_start_idx_eff = max(1, sss_decode_start_idx - demod_stream_start_60 + 1);
else
    demod_read_start_sample = read_start_sample;
    demod_read_length = read_length;
    sss_decode_start_idx_eff = sss_decode_start_idx;
end

if demod_read_start_sample < 0 || demod_read_start_sample >= totalSamples
    error('SSS 解调流读取起点越界: %d / total=%d', demod_read_start_sample, totalSamples);
end
if demod_read_length <= 0
    error('SSS 解调流读取长度无效: %d', demod_read_length);
end

fprintf(['Loading demod stream: %s from %d, len %d ' ...
         '(isolated=%d, sync_start=%d, sss_idx_eff=%d)...\n'], ...
    inFile, demod_read_start_sample, demod_read_length, ...
    enable_isolated_sss_demod_stream, read_start_sample, sss_decode_start_idx_eff);
[x_raw, ~] = read_iq_auto_local(inFile, demod_read_start_sample, demod_read_length, input_format, dat_header_bytes);
x_raw = double(x_raw);

if isempty(x_raw)
    error('SSS 解调阶段读取为空。');
end

if enable_raw_sss_repeat_context
    raw_center_start_60 = sss_decode_start_idx_eff + target_offset;
    raw_center_end_60 = raw_center_start_60 + N_fft - 1;
    raw_center_start = max(1, floor((raw_center_start_60 - 1) * fs_source / fs_target) + 1);
    raw_center_end = min(length(x_raw), ceil((raw_center_end_60 - 1) * fs_source / fs_target) + 1);
    if raw_center_start < 1 || raw_center_end > length(x_raw) || raw_center_end <= raw_center_start
        error('raw SSS repeat context 中心窗口越界: raw=%d:%d, len=%d', ...
            raw_center_start, raw_center_end, length(x_raw));
    end

    x_raw_center = x_raw(raw_center_start:raw_center_end);
    raw_center_len = length(x_raw_center);
    x_raw = [x_raw_center; x_raw_center; x_raw_center];
    sss_decode_start_idx_eff = floor((raw_center_len + 1 - 1) * fs_target / fs_source) + 1;
    fprintf(['Raw SSS repeat context: enable=1, center_raw=%d:%d (%d samples), ' ...
             'new_len=%d, sss_idx_eff=%d\n'], ...
        raw_center_start, raw_center_end, raw_center_len, length(x_raw), sss_decode_start_idx_eff);
else
    fprintf('Raw SSS repeat context: enable=0\n');
end

sss_start_idx_60_for_norm = sss_decode_start_idx_eff + target_offset;
norm_start_60 = max(1, sss_start_idx_60_for_norm - raw_norm_guard_60);
norm_end_60 = sss_start_idx_60_for_norm + N_fft - 1 + raw_norm_guard_60;
norm_start_raw = max(1, floor((norm_start_60 - 1) * fs_source / fs_target) + 1);
norm_end_raw = min(length(x_raw), ceil((norm_end_60 - 1) * fs_source / fs_target) + 1);
[x_raw, raw_norm_info] = normalize_iq_window_local( ...
    x_raw, raw_norm_mode, norm_start_raw, norm_end_raw);
fprintf(['Raw normalization: mode=%s, raw_win=%d:%d (%d samples), ' ...
         'sss60_win=%d:%d, dc=%.4g%+.4gj, scale=%.4g, ' ...
         'fullMeanAbs=%.4g, localMeanAbs=%.4g\n'], ...
    raw_norm_info.mode, raw_norm_info.start_idx, raw_norm_info.end_idx, raw_norm_info.n_win, ...
    norm_start_60, norm_end_60, real(raw_norm_info.dc), imag(raw_norm_info.dc), ...
    raw_norm_info.scale, raw_norm_info.full_mean_abs, raw_norm_info.local_mean_abs);

clip_protect_start_60 = max(1, sss_start_idx_60_for_norm - raw_margin_protect_guard_60);
clip_protect_end_60 = sss_start_idx_60_for_norm + N_fft - 1 + raw_margin_protect_guard_60;
clip_protect_start_raw = max(1, floor((clip_protect_start_60 - 1) * fs_source / fs_target) + 1);
clip_protect_end_raw = min(length(x_raw), ceil((clip_protect_end_60 - 1) * fs_source / fs_target) + 1);
[x_raw, raw_clip_info] = clip_iq_margin_local( ...
    x_raw, enable_raw_margin_clip, raw_margin_clip_factor, ...
    clip_protect_start_raw, clip_protect_end_raw);
fprintf(['Raw margin clip: enable=%d, clip=%.3g, protect_raw=%d:%d, ' ...
         'margin_samples=%d, clipped=%d, maxBefore=%.4g, maxAfter=%.4g\n'], ...
    raw_clip_info.enable, raw_clip_info.clip_level, raw_clip_info.protect_start, ...
    raw_clip_info.protect_end, raw_clip_info.n_margin, raw_clip_info.n_clipped, ...
    raw_clip_info.max_before, raw_clip_info.max_after);

% === 频率补偿参数选择：手动 or 基于当前 PSS 起点盲估 ===
switch lower(freq_comp_mode)
    case 'manual'
        freq_shift_hz_used = freq_shift_hz_manual;
        cfo_hz_used = cfo_hz_manual;

    case 'blind_pilot'
        blind_read_start = max(0, read_start_sample - round(blind_est_pre_start_samples));
        blind_read_len = min(round(blind_est_len), totalSamples - blind_read_start);
        if blind_read_len < 2048
            error('PSS 起点附近剩余样本过短，无法盲估中心频率: pss_start=%d, blind_start=%d, remain=%d', ...
                read_start_sample, blind_read_start, totalSamples - blind_read_start);
        end

        fprintf('Estimating center pilot around PSS start: pss_start=%d, blind_start=%d, len=%d...\n', ...
            read_start_sample, blind_read_start, blind_read_len);
        [x_blind, meta_blind] = read_iq_auto_local(inFile, blind_read_start, blind_read_len, input_format, dat_header_bytes);
        if meta_blind.numSamplesRead < 2048
            error('盲估阶段读取样本过短: got=%d', meta_blind.numSamplesRead);
        end

        x_blind = double(x_blind(:));
        x_blind = x_blind - mean(x_blind);
        x_blind = x_blind / (mean(abs(x_blind)) + eps);

        [f_center_pilot_hz, blind_info] = estimate_center_pilot_blind_from_signal( ...
            x_blind, fs_source, blind_target_bw_hz, blind_bw_tol_ratio, ...
            blind_occ_bg_win_hz, blind_occ_smooth_hz, blind_occ_thresh_db, blind_min_component_bw_hz, ...
            blind_max_expected_cfo_hz, blind_pilot_bg_win_hz, blind_pilot_min_prom_db, blind_pilot_width_ref_hz, ...
            blind_refine_enable, blind_refine_lpf_bw_hz, blind_refine_iters, blind_refine_fir_order, ...
            blind_refine_max_delta_hz);

        freq_shift_hz_used = f_center_pilot_hz;
        cfo_hz_used = 0; % 一步补偿

        fprintf(['Freq mode: BLIND_PILOT, start=%d, center=%.6f MHz ' ...
                 '(equiv CFO vs nominal=%.3f kHz), pilotProm=%.2f dB, ' ...
                 'refineDelta=%.3f kHz, refineUsed=%d\n'], ...
            blind_read_start, freq_shift_hz_used/1e6, ...
            (freq_shift_hz_used-center_nominal_hz)/1e3, ...
            blind_info.prom_db, blind_info.refine_delta_hz/1e3, blind_info.refine_used);
end

%% 3.1 DDC + 可选 CFO 补偿
t_vec = (0:length(x_raw)-1) / fs_source;
fprintf('Shifting spectrum left by %.6f MHz at %.1fMHz...\n', freq_shift_hz_used/1e6, fs_source/1e6);
x_shifted = x_raw .* exp(-1j * 2 * pi * freq_shift_hz_used * t_vec).';

if enable_cfo_comp
    fprintf('Applying CFO correction: %.2f Hz...\n', cfo_hz_used);
    x_base = x_shifted .* exp(-1j * 2 * pi * cfo_hz_used * t_vec).';
else
    fprintf('CFO correction disabled.\n');
    x_base = x_shifted;
end

%% 3.2 SRO 修正（Farrow）
fprintf('Applying SRO Correction (Farrow): %.2f ppm...\n', sro_ppm);
fs_eff = fs_source * (1 + sro_ppm/1e6);

Wn = 32e6 / (fs_source/2);
lpf_order = 30;
b_lpf = fir1(lpf_order, Wn);

switch lower(main_lpf_mode)
    case 'filtfilt'
        fprintf('Main LPF mode: filtfilt (zero-phase bidirectional), order=%d\n', lpf_order);
        x_filt = filtfilt(b_lpf, 1, x_base);

    case 'single_pass'
        lpf_group_delay = lpf_order / 2;
        if abs(lpf_group_delay - round(lpf_group_delay)) > eps
            error('single_pass 模式要求 lpf_order 为偶数，以便补偿整数群延迟。当前 lpf_order=%d', lpf_order);
        end
        lpf_group_delay = round(lpf_group_delay);
        fprintf('Main LPF mode: single_pass FIR, order=%d, group_delay=%d input samples\n', ...
            lpf_order, lpf_group_delay);
        x_filt_delayed = filter(b_lpf, 1, x_base);
        x_filt = [x_filt_delayed(lpf_group_delay+1:end); zeros(lpf_group_delay, 1)];

    otherwise
        error('未知 main_lpf_mode: %s', main_lpf_mode);
end
x_filt = x_filt(:);

T_in = 1 / fs_eff;
T_out = 1 / fs_target;
t_out = 0 : T_out : (length(x_filt)-3)*T_in;

idx_frac = t_out / T_in + 1;
idx_base = floor(idx_frac);
mu = idx_frac - idx_base;

valid_mask = (idx_base >= 2) & (idx_base <= length(x_filt)-2);
idx_base = idx_base(valid_mask);
mu = mu(valid_mask);

h0 = -(mu - 1) .* (mu - 2) .* mu / 6;
h1 =  (mu - 1) .* (mu + 1) .* (mu - 2) / 2;
h2 = -(mu + 1) .* mu .* (mu - 2) / 2;
h3 =  (mu + 1) .* (mu - 1) .* mu / 6;

x_sro = h0 .* x_filt(idx_base - 1) + ...
        h1 .* x_filt(idx_base) + ...
        h2 .* x_filt(idx_base + 1) + ...
        h3 .* x_filt(idx_base + 2);
x_sro = x_sro(:);

%% 3.3 SSS 频域解调（可选参考辅助 offset 选择）
base_start_idx = sss_decode_start_idx_eff;
if enable_ref_aided_select
    offset_candidates = unique(target_offset + ref_offset_search(:).');
else
    offset_candidates = target_offset;
end

best_trial = [];
target_trial = [];
best_score = inf;
fprintf('Reference-aided offset search: enable=%d, candidates=%s, acceptBER=%.4g\n', ...
    enable_ref_aided_select, mat2str(offset_candidates), ref_accept_ber_thresh);

for off_try = offset_candidates
    try
        trial = demod_sss_once_local( ...
            x_sro, base_start_idx, off_try, target_offset, N_fft, ...
            (enable_sss_repeat_context && ~enable_raw_sss_repeat_context), valid_carrier_mode, ...
            valid_carrier_power_threshold_ratio, valid_carrier_dc_guard, ...
            fixed_valid_rel_freq_indices, fixed_valid_fft_bins, ...
            demod_all_subcarriers, decision_rotate_deg, ...
            enable_ref_aided_comp, ber_ref_hex_list, ...
            enable_ref_guided_phase_fit, ref_guided_phase_fit_order, ...
            ref_guided_fit_max_iters, ref_guided_fit_outlier_sigma, ...
            enable_ref_segment_anchor_fit, ref_anchor_seg_len, ...
            ref_anchor_min_coh, ref_anchor_fit_order, ref_anchor_min_segments, ...
            enable_ref_aided_eq, ref_aided_eq_smooth_bins, ...
            unstable_bit_positions, unstable_nibble_positions, false);

        [trial.ref_best_ber, trial.ref_best_rot_idx, trial.ref_best_ref_idx, trial.ref_best_err, trial.ref_best_n] = ...
            score_bits_against_refs_local( ...
                trial.demod_bits_candidates, ber_ref_hex_list, ...
                unstable_bit_positions, unstable_nibble_positions);
        if trial.ref_best_rot_idx >= 1 && trial.ref_best_rot_idx <= numel(trial.rot_angles_deg)
            trial.ref_best_rot_deg = trial.rot_angles_deg(trial.ref_best_rot_idx);
        else
            trial.ref_best_rot_deg = nan;
        end

        if enable_ref_aided_select && ~isnan(trial.ref_best_ber)
            score = trial.ref_best_ber;
        else
            score = off_try ~= target_offset;
        end

        fprintf(['  offset=%+d: scoreBER=%.6g, rot=%g deg, ref=%d, err=%d/%d, ' ...
                 'comp=%s, compBER=%.6g, compCoh=%.3f, slope_coh=%.3f\n'], ...
            off_try, trial.ref_best_ber, trial.ref_best_rot_deg, trial.ref_best_ref_idx, ...
            trial.ref_best_err, trial.ref_best_n, ...
            trial.comp_mode, trial.comp_score_ber, trial.comp_coh, trial.slope_coherence);

        if score < best_score
            best_score = score;
            best_trial = trial;
        end
        if off_try == target_offset
            target_trial = trial;
        end
    catch ME
        fprintf('  offset=%+d: failed (%s)\n', off_try, ME.message);
    end
end

if isempty(best_trial)
    error('所有 SSS offset 候选均解调失败。');
end

if enable_ref_aided_select
    ref_aided_accept = ~isnan(best_trial.ref_best_ber) && best_trial.ref_best_ber <= ref_accept_ber_thresh;
    if ref_aided_accept
        fprintf(['Reference-aided select: ACCEPT offset=%+d, rot=%g deg, ref=%d, ' ...
                 'BER=%.6g (err=%d/%d)\n'], ...
            best_trial.offset, best_trial.ref_best_rot_deg, best_trial.ref_best_ref_idx, ...
            best_trial.ref_best_ber, best_trial.ref_best_err, best_trial.ref_best_n);
    else
        fprintf(['Reference-aided select: REJECT best_offset=%+d, rot=%g deg, ref=%d, ' ...
                 'BER=%.6g (err=%d/%d), threshold=%.6g\n'], ...
            best_trial.offset, best_trial.ref_best_rot_deg, best_trial.ref_best_ref_idx, ...
            best_trial.ref_best_ber, best_trial.ref_best_err, best_trial.ref_best_n, ...
            ref_accept_ber_thresh);
        if ~isempty(target_trial)
            fprintf('Reference-aided select: fallback to original target_offset=%+d.\n', target_offset);
            best_trial = target_trial;
        end
    end
end

offset = best_trial.offset;
sss_start_idx_60 = best_trial.sss_start_idx_60;
x_sss_freq_corr = best_trial.x_sss_freq_corr;
syms_payload = best_trial.syms_payload;
payload_label = best_trial.payload_label;
rot_angles_deg = best_trial.rot_angles_deg;
hex_candidates = best_trial.hex_candidates;
demod_bits_candidates = best_trial.demod_bits_candidates;

fprintf(['Valid carriers: mode=%s, count=%d/%d, segments=%d, longest=%d, ' ...
         'meanP=%.4g, medianP=%.4g, maxP=%.4g, threshold=%.4g\n'], ...
    valid_carrier_mode, best_trial.valid_count, N_fft, ...
    best_trial.valid_segments, best_trial.valid_longest, ...
    best_trial.mean_pwr, best_trial.median_pwr, best_trial.max_pwr, best_trial.threshold_pwr);
fprintf('Valid carrier rel-freq ranges: %s\n', best_trial.valid_ranges);
fprintf('SSS phase fit: slope_coh=%.3f, phase0_coh=%.3f, adjacent_pairs=%d\n', ...
    best_trial.slope_coherence, best_trial.phase0_coherence, best_trial.adjacent_pairs);
fprintf('SSS repeat context: raw=%d, post_sro=%d, center_offset=%+d\n', ...
    enable_raw_sss_repeat_context, (enable_sss_repeat_context && ~enable_raw_sss_repeat_context), target_offset);
fprintf('SSS compensation selected: %s, compBER=%.6g, compRef=%d, compRotIdx=%d, compCoh=%.3f\n', ...
    best_trial.comp_mode, best_trial.comp_score_ber, best_trial.comp_score_ref_idx, ...
    best_trial.comp_score_rot_idx, best_trial.comp_coh);
fprintf('Reference-aided final best BER: %.6g (err=%d/%d), offset=%+d\n', ...
    best_trial.ref_best_ber, best_trial.ref_best_err, best_trial.ref_best_n, best_trial.offset);

%% 3.5 可视化与输出
figure('Position', [300, 300, 640, 560], 'Name', 'SSS Demod (Sync-driven read_start_sample)');
subplot(1,2,1);
syms_plot = x_sss_freq_corr(:) .* exp(1j * decision_rotate_deg * pi/180);
plot(real(syms_plot), imag(syms_plot), 'b.', 'MarkerSize', 8); hold on; grid on; axis square;
th = 0:0.01:2*pi;
plot(cos(th), sin(th), 'k:', 'LineWidth', 0.5);
title(sprintf('Constellation (Rot %d°)', decision_rotate_deg));
xlabel('I'); ylabel('Q'); xlim([-2 2]); ylim([-2 2]);

subplot(1,2,2); axis off;
text(0.02, 0.82, sprintf('Auto read_start = %d\nOffset=%d\n%s', read_start_sample, offset, payload_label), 'FontSize', 11);
text(0.02, 0.54, sprintf('0°  : %s...', hex_candidates{1}(1:min(22,end))), 'FontSize', 10, 'Interpreter','none');
text(0.02, 0.42, sprintf('90° : %s...', hex_candidates{2}(1:min(22,end))), 'FontSize', 10, 'Interpreter','none');
text(0.02, 0.30, sprintf('180°: %s...', hex_candidates{3}(1:min(22,end))), 'FontSize', 10, 'Interpreter','none');
text(0.02, 0.18, sprintf('270°: %s...', hex_candidates{4}(1:min(22,end))), 'FontSize', 10, 'Interpreter','none');

fprintf('\n>>> [Auto Sync + SSS Demod] 完成 <<<\n');
fprintf('Auto read_start_sample = %d\n', read_start_sample);
fprintf('Extracted bits per candidate = %d (%s)\n', length(syms_payload)*2, payload_label);
fprintf('====== HEX CANDIDATES (QPSK 4-way rotation ambiguity) ======\n');
for r = 1:4
    fprintf('[%3d deg] %s\n', rot_angles_deg(r), hex_candidates{r});
end
fprintf('============================================================\n');

%% 3.6 BER 评估（与4组预设正确值对比，支持剔除不稳位）
if enable_ber_eval
    % 检查参考是否已填写
    ref_filled = cellfun(@(s) ~isempty(strtrim(s)), ber_ref_hex_list);
    if ~any(ref_filled)
        fprintf('\n[BER] 未填写 ber_ref_hex_list，跳过 BER 计算。\n');
    else
        nCand = numel(hex_candidates);
        nRef = numel(ber_ref_hex_list);
        ber_mat = nan(nCand, nRef);
        used_bits_mat = zeros(nCand, nRef);

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
                used_bits_mat(i,j) = nKeep;

                fprintf('Cand[%3d deg] vs Ref[%d]: BER=%.6g (err=%d/%d)\n', ...
                    rot_angles_deg(i), j, ber_mat(i,j), err, nKeep);
            end
        end

        % 每个候选最优匹配参考
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
end % end for run_idx

function trial = demod_sss_once_local( ...
    x_sro, base_start_idx, offset, context_center_offset, N_fft, ...
    enable_sss_repeat_context, valid_carrier_mode, ...
    valid_carrier_power_threshold_ratio, valid_carrier_dc_guard, ...
    fixed_valid_rel_freq_indices, fixed_valid_fft_bins, ...
    demod_all_subcarriers, decision_rotate_deg, ...
    enable_ref_aided_comp, ref_hex_list, ...
    enable_ref_guided_phase_fit, ref_guided_phase_fit_order, ...
    ref_guided_fit_max_iters, ref_guided_fit_outlier_sigma, ...
    enable_ref_segment_anchor_fit, ref_anchor_seg_len, ...
    ref_anchor_min_coh, ref_anchor_fit_order, ref_anchor_min_segments, ...
    enable_ref_aided_eq, ref_aided_eq_smooth_bins, ...
    unstable_bit_positions, unstable_nibble_positions, verbose)

if nargin < 29
    verbose = false;
end

sss_start_idx_60 = base_start_idx + offset;
if sss_start_idx_60 < 1 || sss_start_idx_60 + N_fft - 1 > length(x_sro)
    error('SSS 起始点越界: idx=%d, len=%d', sss_start_idx_60, length(x_sro));
end

if enable_sss_repeat_context
    context_center_start = base_start_idx + context_center_offset;
    if context_center_start < 1 || context_center_start + N_fft - 1 > length(x_sro)
        error('SSS repeat context 中心窗口越界: idx=%d, len=%d', context_center_start, length(x_sro));
    end

    x_sss_center = x_sro(context_center_start : context_center_start + N_fft - 1);
    x_sss_repeat = [x_sss_center; x_sss_center; x_sss_center];
    repeat_start = N_fft + 1 + (offset - context_center_offset);
    if repeat_start < 1 || repeat_start + N_fft - 1 > length(x_sss_repeat)
        error('SSS repeat context offset 越界: offset=%d, center=%d', offset, context_center_offset);
    end
    x_sss_time = x_sss_repeat(repeat_start : repeat_start + N_fft - 1);
else
    x_sss_time = x_sro(sss_start_idx_60 : sss_start_idx_60 + N_fft - 1);
end
x_sss_freq = fft(x_sss_time, N_fft) / sqrt(N_fft);

rel_freq = zeros(N_fft,1);
rel_freq(1:N_fft/2) = (0:N_fft/2 - 1).';
rel_freq(N_fft/2+1:N_fft) = (-N_fft/2:-1).';

pwr_bins = abs(x_sss_freq).^2;
mean_pwr = mean(pwr_bins);
median_pwr = median(pwr_bins);
max_pwr = max(pwr_bins);
threshold_pwr = mean_pwr * valid_carrier_power_threshold_ratio;
valid_mask_power = false(N_fft, 1);

switch lower(valid_carrier_mode)
    case 'dynamic'
        valid_mask_power = pwr_bins > threshold_pwr;
        dc_guard = valid_carrier_dc_guard;
        if dc_guard > 0
            valid_mask_power(1:dc_guard) = false;
            valid_mask_power(N_fft-dc_guard+1:N_fft) = false;
        end

    case 'fixed'
        if ~isempty(fixed_valid_fft_bins)
            fixed_bins = unique(round(fixed_valid_fft_bins(:)));
            fixed_bins = fixed_bins(fixed_bins >= 1 & fixed_bins <= N_fft);
            valid_mask_power(fixed_bins) = true;
        elseif ~isempty(fixed_valid_rel_freq_indices)
            fixed_rel = unique(round(fixed_valid_rel_freq_indices(:)));
            valid_mask_power = ismember(rel_freq, fixed_rel);
        else
            error('valid_carrier_mode=fixed 时，需要设置 fixed_valid_rel_freq_indices 或 fixed_valid_fft_bins。');
        end

    otherwise
        error('未知 valid_carrier_mode: %s', valid_carrier_mode);
end

valid_idxs = find(valid_mask_power);
if isempty(valid_idxs)
    error('有效载波检测结果为空。');
end

freq_indices = rel_freq(valid_idxs);
syms_valid = x_sss_freq(valid_idxs);
[freq_indices, sort_idx] = sort(freq_indices);
syms_valid = syms_valid(sort_idx);
valid_idxs = valid_idxs(sort_idx);

df_valid = diff(freq_indices);
seg_start_valid = [1; find(df_valid > 1) + 1];
seg_end_valid = [find(df_valid > 1); length(freq_indices)];
seg_len_valid = seg_end_valid - seg_start_valid + 1;

syms_pow4 = syms_valid.^4;
syms_pow4_unit = syms_pow4 ./ (abs(syms_pow4) + eps);
df = diff(freq_indices);
adjacent_mask = (df == 1);
phase_steps = conj(syms_pow4_unit(1:end-1)) .* syms_pow4_unit(2:end);
phase_steps = phase_steps(adjacent_mask);
if numel(phase_steps) < 8
    error('连续有效载波不足，无法完成相位斜率估计。');
end

phase_step_sum = sum(phase_steps);
slope4 = angle(phase_step_sum);
slope_coherence = abs(phase_step_sum) / numel(phase_steps);
pow4_detrended = syms_pow4_unit .* exp(-1j * slope4 * freq_indices);
phase0_4 = angle(sum(pow4_detrended));
phase0_coherence = abs(mean(pow4_detrended));

full_freq_indices = [(0:N_fft/2-1), (-N_fft/2:-1)].';
phase_correction = (slope4 * full_freq_indices + phase0_4) / 4;
rot_angles_deg = decision_rotate_deg + [0, 90, 180, 270];

comp_trials = {};
x_sss_freq_corr_blind = x_sss_freq .* exp(-1j * phase_correction);
[syms_payload_blind, payload_label_blind] = select_payload_syms_local( ...
    x_sss_freq_corr_blind, valid_idxs, demod_all_subcarriers, N_fft);
[hex_blind, bits_blind] = demod_qpsk_candidates_local( ...
    syms_payload_blind, decision_rotate_deg);
[ber_blind, rot_blind, ref_blind, err_blind, n_blind] = score_bits_against_refs_local( ...
    bits_blind, ref_hex_list, unstable_bit_positions, unstable_nibble_positions);
comp_trials{end+1} = struct( ...
    'mode', 'blind_pow4', ...
    'x_corr', x_sss_freq_corr_blind, ...
    'syms_payload', syms_payload_blind, ...
    'payload_label', payload_label_blind, ...
    'hex_candidates', {hex_blind}, ...
    'demod_bits_candidates', {bits_blind}, ...
    'score_ber', ber_blind, ...
    'score_rot_idx', rot_blind, ...
    'score_ref_idx', ref_blind, ...
    'score_err', err_blind, ...
    'score_n', n_blind, ...
    'slope', slope4 / 4, ...
    'phase0', phase0_4 / 4, ...
    'coh', slope_coherence);

if enable_ref_aided_comp
    for ref_idx = 1:numel(ref_hex_list)
        ref_hex = strtrim(ref_hex_list{ref_idx});
        if isempty(ref_hex)
            continue;
        end
        ref_bits = hex_to_bits_local(ref_hex);
        if length(ref_bits) < 2*N_fft
            continue;
        end

        for rot_idx = 1:numel(rot_angles_deg)
            ref_syms_full = bits_to_qpsk_pre_rotation_local( ...
                ref_bits(1:2*N_fft), decision_rotate_deg, rot_idx);
            ref_syms_valid = ref_syms_full(valid_idxs);

            z_ref = x_sss_freq(valid_idxs) .* conj(ref_syms_valid(:));
            z_ref = z_ref ./ (abs(z_ref) + eps);
            phase_steps_ref = conj(z_ref(1:end-1)) .* z_ref(2:end);
            phase_steps_ref = phase_steps_ref(adjacent_mask);
            if numel(phase_steps_ref) < 8
                continue;
            end

            phase_step_sum_ref = sum(phase_steps_ref);
            slope_ref = angle(phase_step_sum_ref);
            coh_ref = abs(phase_step_sum_ref) / numel(phase_steps_ref);
            phase0_ref = angle(sum(z_ref .* exp(-1j * slope_ref * freq_indices)));

            phase_correction_ref = slope_ref * full_freq_indices + phase0_ref;
            x_corr_ref = x_sss_freq .* exp(-1j * phase_correction_ref);
            [syms_payload_ref, payload_label_ref] = select_payload_syms_local( ...
                x_corr_ref, valid_idxs, demod_all_subcarriers, N_fft);
            [hex_ref, bits_ref] = demod_qpsk_candidates_local( ...
                syms_payload_ref, decision_rotate_deg);
            [ber_ref, rot_ref, best_ref_idx, err_ref, n_ref] = score_bits_against_refs_local( ...
                bits_ref, ref_hex_list, unstable_bit_positions, unstable_nibble_positions);

            comp_trials{end+1} = struct( ...
                'mode', sprintf('ref_aided_ref%d_rot%d', ref_idx, rot_angles_deg(rot_idx)), ...
                'x_corr', x_corr_ref, ...
                'syms_payload', syms_payload_ref, ...
                'payload_label', payload_label_ref, ...
                'hex_candidates', {hex_ref}, ...
                'demod_bits_candidates', {bits_ref}, ...
                'score_ber', ber_ref, ...
                'score_rot_idx', rot_ref, ...
                'score_ref_idx', best_ref_idx, ...
                'score_err', err_ref, ...
                'score_n', n_ref, ...
                'slope', slope_ref, ...
                'phase0', phase0_ref, ...
                'coh', coh_ref);
        end
    end
end

if enable_ref_guided_phase_fit
    for ref_idx = 1:numel(ref_hex_list)
        ref_hex = strtrim(ref_hex_list{ref_idx});
        if isempty(ref_hex)
            continue;
        end
        ref_bits = hex_to_bits_local(ref_hex);
        if length(ref_bits) < 2*N_fft
            continue;
        end

        for rot_idx = 1:numel(rot_angles_deg)
            ref_syms_full = bits_to_qpsk_pre_rotation_local( ...
                ref_bits(1:2*N_fft), decision_rotate_deg, rot_idx);

            [phase_correction_fit, fit_info] = fit_ref_guided_phase_local( ...
                x_sss_freq, ref_syms_full, rel_freq, valid_idxs, ...
                ref_guided_phase_fit_order, ref_guided_fit_max_iters, ref_guided_fit_outlier_sigma);

            x_corr_fit = x_sss_freq .* exp(-1j * phase_correction_fit);
            [syms_payload_fit, payload_label_fit] = select_payload_syms_local( ...
                x_corr_fit, valid_idxs, demod_all_subcarriers, N_fft);
            [hex_fit, bits_fit] = demod_qpsk_candidates_local( ...
                syms_payload_fit, decision_rotate_deg);
            [ber_fit, rot_fit, best_ref_idx_fit, err_fit, n_fit] = score_bits_against_refs_local( ...
                bits_fit, ref_hex_list, unstable_bit_positions, unstable_nibble_positions);

            comp_trials{end+1} = struct( ...
                'mode', sprintf('ref_guided_fit_ref%d_rot%d_ord%d', ref_idx, rot_angles_deg(rot_idx), ref_guided_phase_fit_order), ...
                'x_corr', x_corr_fit, ...
                'syms_payload', syms_payload_fit, ...
                'payload_label', payload_label_fit, ...
                'hex_candidates', {hex_fit}, ...
                'demod_bits_candidates', {bits_fit}, ...
                'score_ber', ber_fit, ...
                'score_rot_idx', rot_fit, ...
                'score_ref_idx', best_ref_idx_fit, ...
                'score_err', err_fit, ...
                'score_n', n_fit, ...
                'slope', fit_info.slope, ...
                'phase0', fit_info.phase0, ...
                'coh', fit_info.coherence);
        end
    end
end

if enable_ref_segment_anchor_fit
    for ref_idx = 1:numel(ref_hex_list)
        ref_hex = strtrim(ref_hex_list{ref_idx});
        if isempty(ref_hex)
            continue;
        end
        ref_bits = hex_to_bits_local(ref_hex);
        if length(ref_bits) < 2*N_fft
            continue;
        end

        for rot_idx = 1:numel(rot_angles_deg)
            ref_syms_full = bits_to_qpsk_pre_rotation_local( ...
                ref_bits(1:2*N_fft), decision_rotate_deg, rot_idx);

            try
                [phase_correction_anchor, anchor_info] = fit_ref_segment_anchor_phase_local( ...
                    x_sss_freq, ref_syms_full, rel_freq, valid_idxs, ...
                    ref_anchor_seg_len, ref_anchor_min_coh, ref_anchor_fit_order, ref_anchor_min_segments);
            catch
                continue;
            end

            x_corr_anchor = x_sss_freq .* exp(-1j * phase_correction_anchor);
            [syms_payload_anchor, payload_label_anchor] = select_payload_syms_local( ...
                x_corr_anchor, valid_idxs, demod_all_subcarriers, N_fft);
            [hex_anchor, bits_anchor] = demod_qpsk_candidates_local( ...
                syms_payload_anchor, decision_rotate_deg);
            [ber_anchor, rot_anchor, best_ref_idx_anchor, err_anchor, n_anchor] = score_bits_against_refs_local( ...
                bits_anchor, ref_hex_list, unstable_bit_positions, unstable_nibble_positions);

            comp_trials{end+1} = struct( ...
                'mode', sprintf('ref_anchor_ref%d_rot%d_seg%d_ord%d', ref_idx, rot_angles_deg(rot_idx), ref_anchor_seg_len, ref_anchor_fit_order), ...
                'x_corr', x_corr_anchor, ...
                'syms_payload', syms_payload_anchor, ...
                'payload_label', payload_label_anchor, ...
                'hex_candidates', {hex_anchor}, ...
                'demod_bits_candidates', {bits_anchor}, ...
                'score_ber', ber_anchor, ...
                'score_rot_idx', rot_anchor, ...
                'score_ref_idx', best_ref_idx_anchor, ...
                'score_err', err_anchor, ...
                'score_n', n_anchor, ...
                'slope', anchor_info.slope, ...
                'phase0', anchor_info.phase0, ...
                'coh', anchor_info.mean_coh);
        end
    end
end

if enable_ref_aided_eq
    for ref_idx = 1:numel(ref_hex_list)
        ref_hex = strtrim(ref_hex_list{ref_idx});
        if isempty(ref_hex)
            continue;
        end
        ref_bits = hex_to_bits_local(ref_hex);
        if length(ref_bits) < 2*N_fft
            continue;
        end

        for rot_idx = 1:numel(rot_angles_deg)
            ref_syms_full = bits_to_qpsk_pre_rotation_local( ...
                ref_bits(1:2*N_fft), decision_rotate_deg, rot_idx);

            h_est = x_sss_freq(:) .* conj(ref_syms_full(:));
            h_unit = h_est ./ (abs(h_est) + eps);
            h_smooth = smooth_freq_response_local(h_unit, rel_freq, ref_aided_eq_smooth_bins);
            x_corr_eq = x_sss_freq .* conj(h_smooth);

            [syms_payload_eq, payload_label_eq] = select_payload_syms_local( ...
                x_corr_eq, valid_idxs, demod_all_subcarriers, N_fft);
            [hex_eq, bits_eq] = demod_qpsk_candidates_local( ...
                syms_payload_eq, decision_rotate_deg);
            [ber_eq, rot_eq, best_ref_idx_eq, err_eq, n_eq] = score_bits_against_refs_local( ...
                bits_eq, ref_hex_list, unstable_bit_positions, unstable_nibble_positions);

            comp_trials{end+1} = struct( ...
                'mode', sprintf('ref_eq_ref%d_rot%d_sm%d', ref_idx, rot_angles_deg(rot_idx), ref_aided_eq_smooth_bins), ...
                'x_corr', x_corr_eq, ...
                'syms_payload', syms_payload_eq, ...
                'payload_label', payload_label_eq, ...
                'hex_candidates', {hex_eq}, ...
                'demod_bits_candidates', {bits_eq}, ...
                'score_ber', ber_eq, ...
                'score_rot_idx', rot_eq, ...
                'score_ref_idx', best_ref_idx_eq, ...
                'score_err', err_eq, ...
                'score_n', n_eq, ...
                'slope', nan, ...
                'phase0', nan, ...
                'coh', abs(mean(h_smooth)));
        end
    end
end

best_comp = comp_trials{1};
for ci = 2:numel(comp_trials)
    cand = comp_trials{ci};
    if isnan(best_comp.score_ber) || (~isnan(cand.score_ber) && cand.score_ber < best_comp.score_ber)
        best_comp = cand;
    end
end

trial = struct();
trial.offset = offset;
trial.sss_start_idx_60 = sss_start_idx_60;
trial.x_sss_freq_corr = best_comp.x_corr;
trial.syms_payload = best_comp.syms_payload;
trial.payload_label = best_comp.payload_label;
trial.rot_angles_deg = rot_angles_deg;
trial.hex_candidates = best_comp.hex_candidates;
trial.demod_bits_candidates = best_comp.demod_bits_candidates;
trial.comp_mode = best_comp.mode;
trial.comp_score_ber = best_comp.score_ber;
trial.comp_score_rot_idx = best_comp.score_rot_idx;
trial.comp_score_ref_idx = best_comp.score_ref_idx;
trial.comp_score_err = best_comp.score_err;
trial.comp_score_n = best_comp.score_n;
trial.comp_slope = best_comp.slope;
trial.comp_phase0 = best_comp.phase0;
trial.comp_coh = best_comp.coh;
trial.valid_count = numel(valid_idxs);
trial.valid_segments = numel(seg_len_valid);
trial.valid_longest = max(seg_len_valid);
trial.valid_ranges = format_index_ranges_local(freq_indices, 12);
trial.mean_pwr = mean_pwr;
trial.median_pwr = median_pwr;
trial.max_pwr = max_pwr;
trial.threshold_pwr = threshold_pwr;
trial.slope_coherence = slope_coherence;
trial.phase0_coherence = phase0_coherence;
trial.adjacent_pairs = numel(phase_steps);

if verbose
    fprintf('SSS offset=%+d: carriers=%d, slope_coh=%.3f, phase0_coh=%.3f\n', ...
        offset, trial.valid_count, slope_coherence, phase0_coherence);
end
end

function [best_ber, best_rot_idx, best_ref_idx, best_err, best_n] = score_bits_against_refs_local( ...
    demod_bits_candidates, ref_hex_list, unstable_bit_positions, unstable_nibble_positions)

best_ber = nan;
best_rot_idx = 0;
best_ref_idx = 0;
best_err = 0;
best_n = 0;

for i = 1:numel(demod_bits_candidates)
    bits_est = demod_bits_candidates{i};
    for j = 1:numel(ref_hex_list)
        ref_hex = strtrim(ref_hex_list{j});
        if isempty(ref_hex)
            continue;
        end

        bits_ref = hex_to_bits_local(ref_hex);
        nUse = min(length(bits_est), length(bits_ref));
        if nUse <= 0
            continue;
        end

        keep_mask = true(nUse,1);
        ub = unstable_bit_positions(:);
        ub = ub(ub >= 1 & ub <= nUse);
        keep_mask(ub) = false;

        un = unstable_nibble_positions(:);
        un = un(un >= 1);
        for k = 1:numel(un)
            b1 = 4*(un(k)-1) + 1;
            b2 = 4*un(k);
            if b1 > nUse
                continue;
            end
            keep_mask(b1:min(b2,nUse)) = false;
        end

        nKeep = nnz(keep_mask);
        if nKeep <= 0
            continue;
        end

        bits_e = bits_est(1:nUse);
        bits_r = bits_ref(1:nUse);
        err = nnz(bits_e(keep_mask) ~= bits_r(keep_mask));
        ber = err / nKeep;

        if isnan(best_ber) || ber < best_ber
            best_ber = ber;
            best_rot_idx = i;
            best_ref_idx = j;
            best_err = err;
            best_n = nKeep;
        end
    end
end
end

function [syms_payload, payload_label] = select_payload_syms_local(x_sss_freq_corr, valid_idxs, demod_all_subcarriers, N_fft)
if demod_all_subcarriers
    syms_payload = x_sss_freq_corr(:);
    payload_label = sprintf('All carriers: %d', N_fft);
else
    syms_payload = x_sss_freq_corr(valid_idxs);
    payload_label = sprintf('Valid carriers: %d', length(valid_idxs));
end
end

function [hex_candidates, demod_bits_candidates] = demod_qpsk_candidates_local(syms_payload, decision_rotate_deg)
hex_candidates = cell(1,4);
demod_bits_candidates = cell(1,4);
for r = 1:4
    syms_payload_rot = syms_payload .* exp(1j * decision_rotate_deg * pi/180) .* exp(-1j * (r-1) * pi/2);
    bits_I = real(syms_payload_rot) < 0;
    bits_Q = imag(syms_payload_rot) < 0;

    demod_bits = zeros(length(syms_payload_rot)*2, 1);
    demod_bits(1:2:end) = bits_I;
    demod_bits(2:2:end) = bits_Q;
    demod_bits_candidates{r} = demod_bits;
    hex_candidates{r} = bits_to_hex_local(demod_bits);
end
end

function syms = bits_to_qpsk_pre_rotation_local(bits, decision_rotate_deg, rot_idx)
bits = bits(:);
n_sym = floor(length(bits) / 2);
bits = bits(1:2*n_sym);
bits_I = bits(1:2:end);
bits_Q = bits(2:2:end);

syms_rot = (1 - 2*double(bits_I)) + 1j * (1 - 2*double(bits_Q));
syms_rot = syms_rot / sqrt(2);

syms = syms_rot .* exp(-1j * decision_rotate_deg * pi/180) .* exp(1j * (rot_idx-1) * pi/2);
syms = syms(:);
end

function [phase_full, info] = fit_ref_guided_phase_local( ...
    x_freq, ref_syms_full, rel_freq, valid_idxs, fit_order, max_iters, outlier_sigma)
% 用已知参考符号去掉 QPSK 数据相位，再对残余相位做 unwrap + 低阶拟合。
x_freq = x_freq(:);
ref_syms_full = ref_syms_full(:);
rel_freq = rel_freq(:);
valid_idxs = valid_idxs(:);

z = x_freq(valid_idxs) .* conj(ref_syms_full(valid_idxs));
z = z ./ (abs(z) + eps);
f = rel_freq(valid_idxs);
f_scale = max(abs(f));
if f_scale < eps
    f_scale = 1;
end
fn = f / f_scale;

ph = unwrap(angle(z));
mask = true(size(ph));
fit_order = max(1, min(round(fit_order), 3));
max_iters = max(1, round(max_iters));

for it = 1:max_iters
    if nnz(mask) <= fit_order + 2
        mask(:) = true;
    end
    p = polyfit(fn(mask), ph(mask), fit_order);
    ph_fit = polyval(p, fn);
    res_circ = angle(exp(1j * (ph - ph_fit)));
    med_res = median(res_circ(mask));
    mad_res = median(abs(res_circ(mask) - med_res)) + eps;
    sigma = 1.4826 * mad_res;
    new_mask = abs(res_circ - med_res) <= max(pi/10, outlier_sigma * sigma);
    if nnz(new_mask) <= fit_order + 2 || isequal(new_mask, mask)
        break;
    end
    mask = new_mask;
end

rel_scale = rel_freq / f_scale;
phase_full = polyval(p, rel_scale);
phase_full = phase_full(:);

res_final = angle(exp(1j * (ph - polyval(p, fn))));
info = struct();
info.coherence = abs(mean(exp(1j * res_final(mask))));
info.n_used = nnz(mask);
info.phase0 = polyval(p, 0);
if fit_order >= 1
    info.slope = p(end-1) / f_scale;
else
    info.slope = 0;
end
end

function [phase_full, info] = fit_ref_segment_anchor_phase_local( ...
    x_freq, ref_syms_full, rel_freq, valid_idxs, seg_len, min_coh, fit_order, min_segments)
% 分段用参考序列去调制，形成可信相位锚点，再用锚点拟合平滑补偿曲线。
x_freq = x_freq(:);
ref_syms_full = ref_syms_full(:);
rel_freq = rel_freq(:);
valid_idxs = valid_idxs(:);
seg_len = max(8, round(seg_len));
fit_order = max(1, min(round(fit_order), 3));
min_segments = max(fit_order + 1, round(min_segments));

f_valid = rel_freq(valid_idxs);
[f_sorted, ord] = sort(f_valid);
idx_sorted = valid_idxs(ord);

anchors_f = [];
anchors_phase = [];
anchors_coh = [];

seg_start = 1;
while seg_start <= numel(idx_sorted)
    seg_end = min(numel(idx_sorted), seg_start + seg_len - 1);
    seg_idx = idx_sorted(seg_start:seg_end);
    seg_f = f_sorted(seg_start:seg_end);

    z = x_freq(seg_idx) .* conj(ref_syms_full(seg_idx));
    z = z ./ (abs(z) + eps);
    z_mean = mean(z);
    coh = abs(z_mean);

    if coh >= min_coh && numel(seg_idx) >= max(8, floor(seg_len/3))
        anchors_f(end+1,1) = mean(seg_f); %#ok<AGROW>
        anchors_phase(end+1,1) = angle(z_mean); %#ok<AGROW>
        anchors_coh(end+1,1) = coh; %#ok<AGROW>
    end

    seg_start = seg_end + 1;
end

if numel(anchors_f) < min_segments
    error('可信相位锚点不足: got=%d, need=%d', numel(anchors_f), min_segments);
end

[anchors_f, ord_anchor] = sort(anchors_f);
anchors_phase = unwrap(anchors_phase(ord_anchor));
anchors_coh = anchors_coh(ord_anchor);

f_scale = max(abs(anchors_f));
if f_scale < eps
    f_scale = 1;
end
fn_anchor = anchors_f / f_scale;
fit_order = min(fit_order, numel(anchors_f)-1);

p = weighted_polyfit_local(fn_anchor, anchors_phase, fit_order, anchors_coh);
phase_full = polyval(p, rel_freq / f_scale);
phase_full = phase_full(:);

res = angle(exp(1j * (anchors_phase - polyval(p, fn_anchor))));
info = struct();
info.n_anchor = numel(anchors_f);
info.mean_coh = mean(anchors_coh);
info.fit_res_rms = sqrt(mean(res.^2));
info.phase0 = polyval(p, 0);
if fit_order >= 1
    info.slope = p(end-1) / f_scale;
else
    info.slope = 0;
end
end

function p = weighted_polyfit_local(x, y, order, w)
x = x(:);
y = y(:);
w = max(w(:), eps);
V = zeros(numel(x), order + 1);
for k = 0:order
    V(:, order + 1 - k) = x.^k;
end
sw = sqrt(w);
p = (V .* sw) \ (y .* sw);
p = p(:).';
end

function h_smooth = smooth_freq_response_local(h, rel_freq, smooth_bins)
h = h(:);
rel_freq = rel_freq(:);
smooth_bins = max(1, round(smooth_bins));
if mod(smooth_bins, 2) == 0
    smooth_bins = smooth_bins + 1;
end

[~, ord] = sort(rel_freq);
h_sorted = h(ord);

if smooth_bins > 1
    h_pad = [repmat(h_sorted(1), floor(smooth_bins/2), 1); ...
             h_sorted; ...
             repmat(h_sorted(end), floor(smooth_bins/2), 1)];
    h_sm_sorted = movmean(h_pad, smooth_bins);
    h_sm_sorted = h_sm_sorted(floor(smooth_bins/2)+1:end-floor(smooth_bins/2));
else
    h_sm_sorted = h_sorted;
end

h_sm_sorted = h_sm_sorted ./ (abs(h_sm_sorted) + eps);
h_smooth = zeros(size(h));
h_smooth(ord) = h_sm_sorted;
end

function hex_str = bits_to_hex_local(bits)
bits = bits(:);
hex_str = '';
for i_hex = 1:4:length(bits)
    chunk = bits(i_hex:min(i_hex+3, length(bits)));
    val = 0;
    for b = 1:length(chunk)
        val = val + chunk(b) * 2^(length(chunk)-b);
    end
    hex_str = [hex_str, dec2hex(val)]; %#ok<AGROW>
end
end

function [f_center_pilot_hz, info] = estimate_center_pilot_blind_from_signal( ...
    x, Fs, target_bw_hz, bw_tol_ratio, occ_bg_win_hz, occ_smooth_hz, occ_thresh_db, min_component_bw_hz, ...
    max_expected_cfo_hz, pilot_bg_win_hz, pilot_min_prom_db, pilot_width_ref_hz, ...
    refine_enable, refine_lpf_bw_hz, refine_iters, refine_fir_order, refine_max_delta_hz)

% 参数固定成稳健默认，避免主流程过重
segLen = 32768;
overlapRatio = 0.75;
nfft = 65536;

N = length(x);
if N < 2048
    error('盲估中心频率失败：样本过短。');
end
if N < segLen * 2
    segLen = 2^floor(log2(max(1024, floor(N/2))));
    segLen = min(segLen, N);
    nfft = max(nfft, 2^nextpow2(segLen));
end

[P_med, ~, f] = robust_welch_psd_local(x, Fs, segLen, overlapRatio, nfft);
% 注意：此时 P_med 已经是对齐 plot_psd.m 的单量纲振幅比值了
Pdb = 20*log10(P_med + eps);

bin_hz = Fs / length(f);
occ_bg_bins = make_odd_local(max(9, round(occ_bg_win_hz / bin_hz)));
occ_sm_bins = make_odd_local(max(5, round(occ_smooth_hz / bin_hz)));

P_occ_base = movmedian(Pdb, occ_bg_bins);
P_occ_env = movmean(P_occ_base, occ_sm_bins);
noise_floor_db = prctile(P_occ_env, 20);
occ_thr_db = noise_floor_db + occ_thresh_db;
occ_mask = P_occ_env > occ_thr_db;

components = mask_to_components_local(occ_mask);
cand = [];
for i = 1:size(components,1)
    iL = components(i,1); iR = components(i,2);
    bw = f(iR) - f(iL);
    if bw < min_component_bw_hz
        continue;
    end

    bw_err = abs(bw - target_bw_hz) / target_bw_hz;
    in_mean = mean(P_occ_env(iL:iR));

    guard = round(1.5e6 / bin_hz);
    oL = max(1, iL-guard):max(1, iL-1);
    oR = min(length(P_occ_env), iR+1):min(length(P_occ_env), iR+guard);
    if isempty(oL), outL = noise_floor_db; else, outL = mean(P_occ_env(oL)); end
    if isempty(oR), outR = noise_floor_db; else, outR = mean(P_occ_env(oR)); end
    edge_drop = in_mean - 0.5*(outL + outR);

    score = 3.0*(1 - min(bw_err,1)) + 0.12*(in_mean - noise_floor_db) + 0.18*max(edge_drop,0);
    if bw_err <= bw_tol_ratio
        cand = [cand; iL, iR, score]; %#ok<AGROW>
    end
end

if isempty(cand)
    % 滑窗兜底
    Nf = length(f);
    bw_bins = max(9, round(target_bw_hz / bin_hz));
    half_bw = floor(bw_bins/2);
    out_half = min(floor(1.5*bw_bins), floor((Nf-1)/2));
    score_sw = -inf(Nf,1);

    for k = (out_half+2):(Nf-out_half-1)
        i1 = k-half_bw; i2 = k+half_bw;
        if i1 < 1 || i2 > Nf, continue; end
        in_band = P_occ_env(i1:i2);

        l1 = max(1, k-out_half); l2 = max(1, i1-1);
        r1 = min(Nf, i2+1);      r2 = min(Nf, k+out_half);
        if l2 < l1 || r2 < r1, continue; end

        out_vals = [P_occ_env(l1:l2); P_occ_env(r1:r2)];
        contrast = mean(in_band) - mean(out_vals);
        score_sw(k) = contrast - 0.15*std(in_band);
    end
    [~, k_best] = max(score_sw);
    f_band_center = f(k_best);
else
    [~, id] = max(cand(:,3));
    iL = cand(id,1); iR = cand(id,2);
    f_band_center = (f(iL) + f(iR)) / 2;
end

search_half_hz = max(1.5e6, max_expected_cfo_hz);
search_half_hz = min(search_half_hz, target_bw_hz/2 - 0.5e6);

[~, f_subbin, pinfo] = detect_center_pilot_once_local( ...
    Pdb, f, f_band_center, search_half_hz, pilot_bg_win_hz, pilot_min_prom_db, pilot_width_ref_hz);

f_center_pilot_hz = f_subbin;

% 时域精修：对粗估频率下变频后窄带滤波，再用相位斜率估计残余频偏
if refine_enable
    f_refined = refine_tone_freq_time_local(x, Fs, f_subbin, refine_lpf_bw_hz, refine_iters, refine_fir_order);
else
    f_refined = f_subbin;
end

refine_delta = f_refined - f_subbin;
if abs(refine_delta) > refine_max_delta_hz
    % 防止时域精修在非理想单音条件下“跑飞”
    f_center_pilot_hz = f_subbin;
    refine_used = 0;
else
    f_center_pilot_hz = f_refined;
    refine_used = 1;
end
info = pinfo;
info.f_band_center = f_band_center;
info.occ_threshold_db = occ_thr_db;
info.f_fft_subbin_hz = f_subbin;
info.f_refined_hz = f_refined;
info.refine_delta_hz = refine_delta;
info.refine_used = refine_used;
end

function [P_med, P_mean, f] = robust_welch_psd_local(x, Fs, segLen, overlapRatio, nfft)
% 完全抛弃 Welch 分段 PSD，直接1:1复刻 plot_psd.m 暴力大点数FFT，赚取单音极大的相干增益
x = x(:);
N = length(x);

% 直接取 2^18 或全长（看齐 plot_psd.m 的 nfftEff = 2^nextpow2(N)）
nfftEff = 2^nextpow2(N);
w = ones(N, 1);

X = fftshift(fft(x .* w, nfftEff));
% 直接用振幅幅度计算dB（不除以 Fs，刻意让单音靠着 N 的成倍暴涨拔高），与 plot_psd.m 完全一致
P_out = abs(X); 
P_out = (P_out / max(P_out)); % 归一化，使得那根最强的针顶格

% 强制把这全尺寸频谱输送出去
P_med = P_out;
P_mean = P_out;
f = ((-nfftEff/2):(nfftEff/2-1)).' * (Fs / nfftEff);
end

function [f_coarse, f_subbin, info] = detect_center_pilot_once_local(Pdb, f, f_center, search_half_hz, bg_win_hz, prom_min_db, width_ref_hz)
df = abs(f(2)-f(1));
idx = abs(f - f_center) <= search_half_hz;
f_roi = f(idx); 
y = Pdb(idx);

if isempty(y)
    error('搜索区间为空。');
end

% 【核心修改】：去掉了所有鲁棒性的多维度记分（抛弃对比背景底噪和峰宽），
% 直接暴力寻找指定带宽窗口内，PSD 幅值最高的那一根单音主峰！
[max_val, ksel] = max(y);
f_coarse = f_roi(ksel);

% 保留三点抛物线插值，用于在离散频率格点之间再提取一步小数精确度
if ksel <= 1 || ksel >= length(y)
    f_subbin = f_coarse;
else
    y1 = y(ksel-1); y2 = y(ksel); y3 = y(ksel+1);
    den = (y1 - 2*y2 + y3);
    delta = (abs(den)<1e-12)*0 + (abs(den)>=1e-12)*(0.5*(y1-y3)/den);
    f_subbin = f_roi(ksel) + max(min(delta, 1), -1) * df;
end
info = struct('prom_db', max_val, 'width_hz', 0, 'score', max_val);
end

function comps = mask_to_components_local(mask)
mask = logical(mask(:));
d = diff([false; mask; false]);
st = find(d == 1);
ed = find(d == -1) - 1;
comps = [st, ed];
end

function n = make_odd_local(n)
n = round(n);
if mod(n,2)==0
    n = n + 1;
end
n = max(3, n);
end

function bits = hex_to_bits_local(hex_str)
% HEX字符串 -> bit列向量（MSB first，每个hex字符4bit）
if isempty(hex_str)
    bits = zeros(0,1);
    return;
end

hex_str = upper(regexprep(hex_str, '\s+', ''));
hex_str = regexprep(hex_str, '[^0-9A-F]', '');
if isempty(hex_str)
    bits = zeros(0,1);
    return;
end

vals = zeros(length(hex_str),1);
for i = 1:length(hex_str)
    vals(i) = hex2dec(hex_str(i));
end

bits = zeros(4*length(vals),1);
for i = 1:length(vals)
    v = vals(i);
    bits(4*(i-1)+1) = bitget(v,4);
    bits(4*(i-1)+2) = bitget(v,3);
    bits(4*(i-1)+3) = bitget(v,2);
    bits(4*(i-1)+4) = bitget(v,1);
end
end

function f_refined = refine_tone_freq_time_local(x, Fs, f_init_hz, lpf_bw_hz, iters, fir_order)
% 基于相位增量的时域残余频偏精修
x = x(:);
N = length(x);
n = (0:N-1).';
f_refined = f_init_hz;

for it = 1:max(1, iters)
    z = x .* exp(-1j * 2*pi * f_refined * n / Fs);

    % 窄带低通，隔离中心单音
    bw_now = lpf_bw_hz / (2^(it-1));
    Wn = min(0.95, max(1e-5, bw_now / (Fs/2)));
    b = fir1(fir_order, Wn);
    zf = filtfilt(b, 1, z);

    % 去掉滤波过渡段，降低边界影响
    trim = min(length(zf)-2, max(8, fir_order));
    if trim*2 >= length(zf)-1
        z_use = zf;
    else
        z_use = zf(trim+1:end-trim);
    end

    if length(z_use) < 4
        break;
    end

    % 残余频率：arg(sum(conj(z[n])z[n+1])) * Fs/(2pi)
    r = conj(z_use(1:end-1)) .* z_use(2:end);
    dphi = angle(sum(r));
    df_hz = dphi * Fs / (2*pi);

    % 迭代修正
    f_refined = f_refined + df_hz;

    % 已很小则提前结束
    if abs(df_hz) < 5
        break;
    end
end

end

function [x_norm, info] = normalize_iq_window_local(x, mode, win_start, win_end)
% 对 IQ 做可切换归一化；局部鲁棒模式只让目标窗口附近样本参与统计。
x = x(:);
n = length(x);
win_start = max(1, min(n, round(win_start)));
win_end = max(win_start, min(n, round(win_end)));
mode_l = lower(mode);

info = struct();
info.mode = mode_l;
info.start_idx = win_start;
info.end_idx = win_end;
info.n_win = win_end - win_start + 1;
info.full_mean_abs = mean(abs(x));
info.local_mean_abs = mean(abs(x(win_start:win_end)));

switch mode_l
    case 'full_mean'
        dc = mean(x);
        scale = mean(abs(x - dc));

    case 'sss_local_robust'
        x_win = x(win_start:win_end);
        dc = median(real(x_win)) + 1j * median(imag(x_win));
        scale = median(abs(x_win - dc));
        if scale < eps
            scale = mean(abs(x_win - dc));
        end

    otherwise
        error('未知 raw_norm_mode: %s', mode);
end

if scale < eps
    scale = mean(abs(x - mean(x)));
end
if scale < eps
    scale = 1;
end

x_norm = (x - dc) / scale;
info.dc = dc;
info.scale = scale;
end

function [x_out, info] = clip_iq_margin_local(x, enable_clip, clip_level, protect_start, protect_end)
% 只对保护区外的余量样本做幅度限幅，保留样本相位。
x_out = x(:);
n = length(x_out);
protect_start = max(1, min(n, round(protect_start)));
protect_end = max(protect_start, min(n, round(protect_end)));

margin_mask = true(n, 1);
margin_mask(protect_start:protect_end) = false;
mag_before = abs(x_out);

info = struct();
info.enable = logical(enable_clip);
info.clip_level = clip_level;
info.protect_start = protect_start;
info.protect_end = protect_end;
info.n_margin = nnz(margin_mask);
info.max_before = max(mag_before);

if enable_clip && clip_level > 0 && any(margin_mask)
    margin_idx = find(margin_mask);
    mag_margin = mag_before(margin_idx);
    clip_mask = mag_margin > clip_level;
    if any(clip_mask)
        idx_clip = margin_idx(clip_mask);
        x_out(idx_clip) = x_out(idx_clip) .* (clip_level ./ (mag_before(idx_clip) + eps));
    end
    info.n_clipped = nnz(clip_mask);
else
    info.n_clipped = 0;
end

info.max_after = max(abs(x_out));
end

function s = format_index_ranges_local(idx, max_ranges)
% 将连续整数序列压缩成 a:b 形式，避免诊断打印过长。
if nargin < 2 || isempty(max_ranges)
    max_ranges = 12;
end

idx = unique(idx(:).');
if isempty(idx)
    s = '(empty)';
    return;
end

d = diff(idx);
seg_start = [1, find(d > 1) + 1];
seg_end = [find(d > 1), numel(idx)];
n_seg = numel(seg_start);
n_show = min(n_seg, max_ranges);
parts = cell(1, n_show);

for k = 1:n_show
    a = idx(seg_start(k));
    b = idx(seg_end(k));
    if a == b
        parts{k} = sprintf('%d', a);
    else
        parts{k} = sprintf('%d:%d', a, b);
    end
end

s = strjoin(parts, ', ');
if n_seg > max_ranges
    s = sprintf('%s, ... (%d ranges total)', s, n_seg);
end
end

function [x, meta] = read_iq_auto_local(filename, startSample, numSamples, input_format, dat_header_bytes)
% 支持 .iq / .dat 读取（int16 LE, I/Q 交织）

if nargin < 5 || isempty(dat_header_bytes)
    dat_header_bytes = 0;
end

if nargin < 4 || isempty(input_format)
    input_format = 'auto';
end

[~, ~, ext] = fileparts(filename);
ext = lower(ext);
fmt = lower(input_format);

if strcmp(fmt, 'auto')
    if strcmp(ext, '.dat')
        fmt = 'dat';
    else
        fmt = 'iq';
    end
end

switch fmt
    case 'iq'
        [x, meta] = iq_read_int16_le(filename, startSample, numSamples);

    case 'dat'
        % .dat 按 int16 LE I/Q 交织读取；支持可选文件头
        [x, meta] = iq_read_int16_le(filename, startSample, numSamples, dat_header_bytes);

    otherwise
        error('不支持的 input_format: %s (可选: auto/iq/dat)', input_format);
end
end
