% estimate_center_pilot_blind.m
% 盲估计“目标60MHz占用带中心单载波”频率（工程实验版）
% 输出：中心单载波频率估计值（Hz / MHz）
%
% 两阶段：
% A) 全带稳健频谱 + 连续占用带检测（约60MHz）
% B) 仅在占用带中心附近检测窄峰，并做亚bin精估与多块稳健融合

clear; clc; close all;

%% ======================= 参数区 =======================
inFile = 'sigtest58.iq';
Fs = 409.6e6;

startSample = 0;          % 0-based
Nread = 2^21;             % 建议 >= 1e6
removeMean = true;        % 去直流
normalizeAmp = true;

% --- 先验：目标占用带约60MHz ---
target_bw_hz = 60e6;
bw_tol_ratio = 0.35;      % 候选带宽容差 ±35%

% --- 稳健频谱估计参数（Welch + 中值融合） ---
segLen = 32768;
overlapRatio = 0.75;
nfft = 65536;

% --- 占用带检测参数 ---
occ_bg_win_hz = 2.0e6;    % 抑制窄峰的中值窗口
occ_smooth_hz = 1.2e6;    % 连续包络平滑窗口
occ_thresh_db = 5.5;      % 相对底噪门限(dB)
min_component_bw_hz = 20e6;

% --- 中心附近窄峰检测参数 ---
max_expected_cfo_hz = 6e6;      % 若不确定可放大，但越大越易引入干扰
pilot_bg_win_hz = 0.9e6;        % 局部背景窗口
pilot_min_prom_db = 4.5;        % 峰显著度门限
pilot_width_ref_hz = 0.25e6;    % 惩罚参考峰宽（越窄越好）

% --- 多时间块一致性 ---
num_consistency_blocks = 8;      % 每块独立估计再融合
min_valid_blocks = 3;

plot_debug = true;

%% ======================= 1) 读取IQ =======================
[x, meta] = iq_read_int16_le(inFile, startSample, Nread);
Ngot = meta.numSamplesRead;
if Ngot < 2048
    error('样本过少：读取到 %d 点，至少需要 2048 点。', Ngot);
end

% 若样本较短，自动缩小 Welch 参数（避免直接报错）
if Ngot < segLen * 2
    segLen_new = 2^floor(log2(max(1024, floor(Ngot/2))));
    segLen_new = min(segLen_new, Ngot);
    nfft_new = max(nfft, 2^nextpow2(segLen_new));

    warning(['样本不足以使用原参数，自动调整: segLen %d -> %d, nfft %d -> %d ' ...
             '(N=%d)'], segLen, segLen_new, nfft, nfft_new, Ngot);
    segLen = segLen_new;
    nfft = nfft_new;
end

x = double(x(:));
if removeMean
    x = x - mean(x);
end
if normalizeAmp
    x = x / (mean(abs(x)) + eps);
end

fprintf('Read IQ: %s, start=%d, req=%d, got=%d (segLen=%d, nfft=%d)\n', ...
    inFile, startSample, Nread, Ngot, segLen, nfft);

%% ======================= 2) 全带稳健频谱 =======================
[P_med, P_mean, f] = robust_welch_psd(x, Fs, segLen, overlapRatio, nfft);
Pdb = 10*log10(P_med + eps);

%% ======================= 3) 阶段A：占用带检测 =======================
bin_hz = Fs / nfft;
occ_bg_bins = make_odd(max(9, round(occ_bg_win_hz / bin_hz)));
occ_sm_bins = make_odd(max(5, round(occ_smooth_hz / bin_hz)));

% 先用 movmedian 压窄峰，再做平滑提连续包络
P_occ_base = movmedian(Pdb, occ_bg_bins);
P_occ_env = movmean(P_occ_base, occ_sm_bins);

noise_floor_db = prctile(P_occ_env, 20);
occ_thr_db = noise_floor_db + occ_thresh_db;
occ_mask = P_occ_env > occ_thr_db;

components = mask_to_components(occ_mask);
if isempty(components)
    error('未检测到占用带候选（可尝试降低 occ_thresh_db）。');
end

cand = [];          % 严格满足带宽容差的候选
for i = 1:size(components,1)
    iL = components(i,1); iR = components(i,2);
    bw = f(iR) - f(iL);
    if bw < min_component_bw_hz
        continue;
    end

    bw_err = abs(bw - target_bw_hz) / target_bw_hz;

    in_mean = mean(P_occ_env(iL:iR));
    in_peak = max(P_occ_env(iL:iR));

    % 边缘回落（工程上帮助排除“全带平坦抬高”）
    guard = round(1.5e6 / bin_hz);
    oL = max(1, iL-guard):max(1, iL-1);
    oR = min(length(P_occ_env), iR+1):min(length(P_occ_env), iR+guard);
    if isempty(oL), outL = noise_floor_db; else, outL = mean(P_occ_env(oL)); end
    if isempty(oR), outR = noise_floor_db; else, outR = mean(P_occ_env(oR)); end
    edge_drop = in_mean - 0.5*(outL + outR);

    % 打分：带宽匹配 + 带内能量 + 边缘回落
    score = 3.0*(1 - min(bw_err,1)) + 0.12*(in_mean - noise_floor_db) + 0.18*max(edge_drop,0);

    if bw_err > bw_tol_ratio
        continue;
    end

    cand = [cand; iL, iR, bw, in_mean, in_peak, edge_drop, score]; %#ok<AGROW>
end

if isempty(cand)
    % === 降级策略：固定60MHz滑窗对比打分（更抗“宽平台+窄峰干扰”）===
    Nf = length(f);
    bw_bins = max(9, round(target_bw_hz / bin_hz));
    half_bw = floor(bw_bins/2);
    out_half = min(floor(1.5*bw_bins), floor((Nf-1)/2));

    score_sw = -inf(Nf,1);
    mean_in_all = nan(Nf,1);

    for k = (out_half+2):(Nf-out_half-1)
        i1 = k-half_bw; i2 = k+half_bw;
        if i1 < 1 || i2 > Nf
            continue;
        end

        in_band = P_occ_env(i1:i2);
        mean_in = mean(in_band);

        l1 = max(1, k-out_half); l2 = max(1, i1-1);
        r1 = min(Nf, i2+1);      r2 = min(Nf, k+out_half);

        if l2 < l1 || r2 < r1
            continue;
        end

        out_vals = [P_occ_env(l1:l2); P_occ_env(r1:r2)];
        mean_out = mean(out_vals);

        % 对比打分：带内抬升 + 带内稳定性（抑制纯尖峰）
        contrast = mean_in - mean_out;
        flat_penalty = std(in_band);
        score_sw(k) = contrast - 0.15*flat_penalty;
        mean_in_all(k) = mean_in;
    end

    [best_sw, k_best] = max(score_sw);
    if ~isfinite(best_sw)
        error('严格候选为空，且滑窗法失败（可降低 occ_thresh_db 或增大 Nread）。');
    end

    f_band_center = f(k_best);
    fL = f_band_center - target_bw_hz/2;
    fR = f_band_center + target_bw_hz/2;

    warning(['无候选满足严格60MHz容差，已启用滑窗带选择: center=%.6f MHz, ' ...
             'score=%.2f dB, inMean=%.2f dB'], ...
            f_band_center/1e6, best_sw, mean_in_all(k_best));
else
    [~, id_best] = max(cand(:,7));
    iL = cand(id_best,1); iR = cand(id_best,2);
    fL = f(iL); fR = f(iR);
    f_band_center = (fL + fR)/2;
end

fprintf('\n[Stage-A] Selected band:\n');
fprintf('  fL = %.6f MHz, fR = %.6f MHz, BW = %.3f MHz\n', fL/1e6, fR/1e6, (fR-fL)/1e6);
fprintf('  Band center = %.6f MHz\n', f_band_center/1e6);

%% ======================= 4) 阶段B：中心附近单载波 =======================
search_half_hz = max(1.5e6, max_expected_cfo_hz);
search_half_hz = min(search_half_hz, target_bw_hz/2 - 0.5e6); % 不能扩到带外过多

[pilot_coarse_hz, pilot_subbin_hz, peak_info] = detect_center_pilot_once( ...
    Pdb, f, f_band_center, search_half_hz, pilot_bg_win_hz, pilot_min_prom_db, pilot_width_ref_hz);

fprintf('\n[Stage-B] Pilot (global PSD):\n');
fprintf('  Coarse peak = %.6f MHz\n', pilot_coarse_hz/1e6);
fprintf('  Sub-bin peak = %.6f MHz\n', pilot_subbin_hz/1e6);
fprintf('  Prominence = %.2f dB, Width = %.1f kHz, Score = %.3f\n', ...
    peak_info.prom_db, peak_info.width_hz/1e3, peak_info.score);

%% ======================= 5) 多块一致性融合 =======================
blk_freq = nan(num_consistency_blocks,1);
blk_score = nan(num_consistency_blocks,1);

N = length(x);
blk_len = floor(N / num_consistency_blocks);
if blk_len < segLen
    warning('块太短，自动降低一致性块数。');
    num_consistency_blocks = max(1, floor(N / segLen));
    blk_freq = nan(num_consistency_blocks,1);
    blk_score = nan(num_consistency_blocks,1);
    blk_len = floor(N / num_consistency_blocks);
end

for b = 1:num_consistency_blocks
    s = (b-1)*blk_len + 1;
    e = min(N, b*blk_len);
    xb = x(s:e);
    if length(xb) < segLen
        continue;
    end

    [Pb_med, ~, fb] = robust_welch_psd(xb, Fs, segLen, overlapRatio, nfft);
    Pb_db = 10*log10(Pb_med + eps);

    try
        [~, f_sub, info_b] = detect_center_pilot_once( ...
            Pb_db, fb, f_band_center, search_half_hz, pilot_bg_win_hz, pilot_min_prom_db, pilot_width_ref_hz);
        blk_freq(b) = f_sub;
        blk_score(b) = info_b.score;
    catch
        % 无有效峰，留NaN
    end
end

valid = ~isnan(blk_freq);
if nnz(valid) < min_valid_blocks
    warning('有效时间块不足（%d/%d），退化使用全局估计。', nnz(valid), num_consistency_blocks);
    f_center_pilot_hz = pilot_subbin_hz;
else
    f_valid = blk_freq(valid);
    s_valid = blk_score(valid);

    f_med = median(f_valid);
    mad_val = median(abs(f_valid - f_med));
    gate = max(2*bin_hz, 3*max(mad_val, bin_hz));
    inlier = abs(f_valid - f_med) <= gate;

    if nnz(inlier) < max(2, floor(0.5*nnz(valid)))
        f_center_pilot_hz = f_med;
    else
        w = max(s_valid(inlier), 1e-3);
        f_center_pilot_hz = sum(f_valid(inlier).*w) / sum(w);
    end
end

fprintf('\n================ FINAL ESTIMATION ================\n');
fprintf('Estimated center pilot frequency = %.6f MHz\n', f_center_pilot_hz/1e6);
fprintf('Estimated center pilot frequency = %.3f Hz\n', f_center_pilot_hz);
fprintf('==================================================\n');

%% ======================= 6) 调试图 =======================
if plot_debug
    figure('Name','Blind Center-Pilot Estimation','Position',[120 80 1400 760]);

    subplot(2,2,1);
    plot(f/1e6, Pdb, 'Color',[0.3 0.3 0.8]); hold on;
    plot(f/1e6, P_occ_env, 'k', 'LineWidth', 1.1);
    yline(occ_thr_db, 'g--');
    xline(fL/1e6, 'r--'); xline(fR/1e6, 'r--');
    xline(f_band_center/1e6, 'm-');
    grid on; xlabel('Frequency (MHz)'); ylabel('PSD (dB)');
    title('Stage-A: Occupied Band Detection');
    legend('Robust PSD','Occupancy Envelope','Threshold','fL/fR','Band Center','Location','best');

    subplot(2,2,2);
    idx = abs(f - f_band_center) <= search_half_hz;
    f_roi = f(idx);
    P_roi = Pdb(idx);
    bg_roi = movmedian(P_roi, make_odd(max(9, round(pilot_bg_win_hz / bin_hz))));
    res_roi = P_roi - bg_roi;
    plot(f_roi/1e6, res_roi, 'b'); hold on;
    xline(f_center_pilot_hz/1e6, 'r-', 'LineWidth', 1.2);
    yline(pilot_min_prom_db, 'k--');
    grid on; xlabel('Frequency (MHz)'); ylabel('Residual over local baseline (dB)');
    title('Stage-B: Center-Window Residual Peak');

    subplot(2,2,3);
    stem(1:num_consistency_blocks, blk_freq/1e6, 'filled'); hold on;
    yline(f_center_pilot_hz/1e6, 'r-', 'LineWidth', 1.2);
    grid on; xlabel('Block Index'); ylabel('Estimated Frequency (MHz)');
    title('Multi-block Consistency');

    subplot(2,2,4);
    stem(1:num_consistency_blocks, blk_score, 'filled');
    grid on; xlabel('Block Index'); ylabel('Peak Score');
    title('Block Peak Scores');
end

%% ======================= 本地函数 =======================
function [P_med, P_mean, f] = robust_welch_psd(x, Fs, segLen, overlapRatio, nfft)
    x = x(:);
    N = length(x);
    hop = max(1, round(segLen*(1-overlapRatio)));
    if N < segLen
        error('N(%d) < segLen(%d)', N, segLen);
    end

    nSeg = floor((N-segLen)/hop) + 1;
    P = zeros(nfft, nSeg);

    n = (0:segLen-1).';
    w = 0.5 - 0.5*cos(2*pi*n/(segLen-1)); % Hann
    U = sum(w.^2);

    c = 0;
    for s = 1:hop:(N-segLen+1)
        c = c + 1;
        seg = x(s:s+segLen-1).*w;
        X = fft(seg, nfft);
        Px = (abs(X).^2)/(Fs*U);
        P(:,c) = fftshift(Px);
    end

    P_med = median(P, 2);
    P_mean = mean(P, 2);
    f = ((-nfft/2):(nfft/2-1)).' * (Fs/nfft);
end

function [f_coarse, f_subbin, info] = detect_center_pilot_once(Pdb, f, f_center, search_half_hz, bg_win_hz, prom_min_db, width_ref_hz)
    df = abs(f(2)-f(1));
    idx = abs(f - f_center) <= search_half_hz;
    if nnz(idx) < 9
        error('中心搜索窗口过窄或频轴异常。');
    end

    f_roi = f(idx);
    y = Pdb(idx);

    bg_bins = make_odd(max(9, round(bg_win_hz/df)));
    y_bg = movmedian(y, bg_bins);
    r = y - y_bg; % 残差谱

    % 局部峰检测（不依赖 findpeaks）
    pk = false(size(r));
    pk(2:end-1) = (r(2:end-1) > r(1:end-2)) & (r(2:end-1) >= r(3:end));
    idx_pk = find(pk);
    if isempty(idx_pk)
        error('中心窗口内未发现局部峰。');
    end

    prom = r(idx_pk);

    % 峰宽（以残差半高近似）
    widths = zeros(size(idx_pk));
    dist_center = abs(f_roi(idx_pk) - f_center);
    for i = 1:length(idx_pk)
        k = idx_pk(i);
        halfv = max(0, 0.5*prom(i));
        L = k;
        while L > 1 && r(L) > halfv, L = L - 1; end
        R = k;
        while R < length(r) && r(R) > halfv, R = R + 1; end
        widths(i) = (R - L + 1) * df;
    end

    % 打分：显著度高 + 峰宽窄 + 更靠近带中心
    prom_n = prom / 10;                    % 约每10dB归一化1
    width_n = widths / max(width_ref_hz, df);
    dist_n = dist_center / max(search_half_hz, df);
    score = 2.0*prom_n - 0.6*width_n - 1.2*dist_n;

    % 先去掉显著度太弱
    valid = prom >= prom_min_db;
    if ~any(valid)
        error('未找到满足显著度门限的中心峰。');
    end

    idx_cand = find(valid);
    [~, im] = max(score(idx_cand));
    ksel = idx_pk(idx_cand(im));

    f_coarse = f_roi(ksel);

    % 亚bin抛物线插值（log谱）
    if ksel <= 1 || ksel >= length(y)
        f_subbin = f_coarse;
    else
        y1 = y(ksel-1); y2 = y(ksel); y3 = y(ksel+1);
        den = (y1 - 2*y2 + y3);
        if abs(den) < 1e-12
            delta = 0;
        else
            delta = 0.5 * (y1 - y3) / den;
            delta = max(min(delta, 1), -1);
        end
        f_subbin = f_roi(ksel) + delta * df;
    end

    info = struct();
    info.prom_db = r(ksel);
    info.width_hz = widths(idx_cand(im));
    info.score = score(idx_cand(im));
end

function comps = mask_to_components(mask)
    mask = logical(mask(:));
    d = diff([false; mask; false]);
    st = find(d == 1);
    ed = find(d == -1) - 1;
    comps = [st, ed];
end

function n = make_odd(n)
    n = round(n);
    if mod(n,2)==0
        n = n + 1;
    end
    n = max(3, n);
end
