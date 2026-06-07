function result = sss_ref_locked_demod_core(x_sro, base_start_idx, opts)
%SSS_REF_LOCKED_DEMOD_CORE  Reference-locked SSS demodulation core.
%
% 这个函数提炼 sync_then_sss_demod.m 中当前最有效的 SSS 解调方法：
%   1) 在若干候选 FFT 起点 offset 上提取 SSS 窗口；
%   2) 先做非参考的 QPSK 四次方线性相位补偿；
%   3) 再用已知 SSS 候选序列参与相位斜率/相位曲线/分段锚点拟合；
%   4) 对所有补偿候选做 QPSK 判决，并用 BER 选择；
%   5) 最后做频域子载波整体 circular shift 锁定，解决“星座正确但序列错位”的问题。
%
% 注意：这里不会逐比特改写结果。参考序列只用于估计补偿参数和选择整体锁定状态。
%
% 输入:
%   x_sro          DDC/低通/重采样后的 60 MHz 复数基带序列
%   base_start_idx SSS 标称起点，1-based，单位为 60 MHz 样点
%   opts           参数结构体，常用字段见 default_opts_local()
%
% 输出:
%   result         最优 offset / compensation / carrierShift / HEX / BER 等诊断结果

opts = merge_opts_local(default_opts_local(), opts);
x_sro = x_sro(:);

offset_candidates = unique(opts.target_offset + opts.offset_search(:).');
best = [];
best_score = inf;

for off = offset_candidates
    try
        trial = demod_one_offset_local(x_sro, base_start_idx, off, opts);
    catch ME
        trial = struct();
        trial.failed = true;
        trial.error = ME.message;
        trial.offset = off;
    end

    if isfield(trial, 'failed') && trial.failed
        if opts.verbose
            fprintf('offset=%+d failed: %s\n', off, trial.error);
        end
        continue;
    end

    score = trial.best_ber;
    if isnan(score)
        score = inf;
    end

    if opts.verbose
        fprintf(['offset=%+d: BER=%.6g, err=%d/%d, ref=%d, rot=%d deg, ' ...
                 'comp=%s, carrierShift=%+d, rawBER=%.6g, shiftedBER=%.6g\n'], ...
            trial.offset, trial.best_ber, trial.best_err, trial.best_n, ...
            trial.best_ref_idx, trial.best_rot_deg, trial.comp_mode, ...
            trial.carrier_shift, trial.carrier_shift_raw_ber, trial.carrier_shift_ber);
    end

    if score < best_score
        best_score = score;
        best = trial;
    end
end

if isempty(best)
    error('所有 SSS offset 候选均解调失败。');
end

result = best;
result.accepted = ~isnan(best.best_ber) && best.best_ber <= opts.accept_ber_thresh;
end

function opts = default_opts_local()
opts = struct();
opts.N_fft = 1024;
opts.target_offset = 0;
opts.offset_search = -8:8;
opts.decision_rotate_deg = 45;

% 有效载波。当前实验中固定为两段，避开 DC 附近空洞。
opts.fixed_valid_rel_freq_indices = [-512:-10, 7:511];

% 参考辅助补偿开关。
opts.enable_blind_pow4 = true;
opts.enable_ref_aided_linear = true;
opts.enable_ref_guided_phase_fit = true;
opts.ref_guided_phase_fit_order = 2;
opts.ref_guided_fit_max_iters = 3;
opts.ref_guided_fit_outlier_sigma = 3.0;
opts.enable_ref_segment_anchor_fit = true;
opts.ref_anchor_seg_len = 48;
opts.ref_anchor_min_coh = 0.10;
opts.ref_anchor_fit_order = 2;
opts.ref_anchor_min_segments = 5;

% 最关键的新增锁定：补偿后允许频域载波序列整体 circular shift。
opts.enable_ref_carrier_shift_lock = true;
opts.ref_carrier_shift_search = -32:32;
opts.ref_carrier_shift_min_gain = 0.02;

% BER 和诊断。
opts.ref_hex_list = {};
opts.unstable_bit_positions = [];
opts.unstable_nibble_positions = [];
opts.accept_ber_thresh = 0.03;
opts.verbose = true;
end

function trial = demod_one_offset_local(x_sro, base_start_idx, offset, opts)
N_fft = opts.N_fft;
sss_start = base_start_idx + offset;
if sss_start < 1 || sss_start + N_fft - 1 > length(x_sro)
    error('SSS 起始点越界: idx=%d, len=%d', sss_start, length(x_sro));
end

x_time = x_sro(sss_start:sss_start+N_fft-1);
x_freq = fft(x_time, N_fft) / sqrt(N_fft);

rel_freq = rel_freq_local(N_fft);
valid_mask = ismember(rel_freq, opts.fixed_valid_rel_freq_indices(:));
valid_idxs = find(valid_mask);
if isempty(valid_idxs)
    error('有效载波为空。');
end

freq_indices = rel_freq(valid_idxs);
[freq_indices, ord] = sort(freq_indices);
valid_idxs = valid_idxs(ord);
adjacent_mask = diff(freq_indices) == 1;

rot_angles_deg = opts.decision_rotate_deg + [0, 90, 180, 270];
comp_trials = {};

if opts.enable_blind_pow4
    comp_trials{end+1} = make_blind_pow4_trial_local( ...
        x_freq, rel_freq, valid_idxs, freq_indices, adjacent_mask, opts);
end

if opts.enable_ref_aided_linear || opts.enable_ref_guided_phase_fit || opts.enable_ref_segment_anchor_fit
    for ref_idx = 1:numel(opts.ref_hex_list)
        ref_hex = strtrim(opts.ref_hex_list{ref_idx});
        if isempty(ref_hex)
            continue;
        end

        ref_bits = hex_to_bits_local(ref_hex);
        if length(ref_bits) < 2*N_fft
            continue;
        end

        for rot_idx = 1:numel(rot_angles_deg)
            ref_syms = bits_to_qpsk_pre_rotation_local( ...
                ref_bits(1:2*N_fft), opts.decision_rotate_deg, rot_idx);

            if opts.enable_ref_aided_linear
                comp_trials{end+1} = make_ref_linear_trial_local( ...
                    x_freq, ref_syms, rel_freq, valid_idxs, freq_indices, ...
                    adjacent_mask, opts, ref_idx, rot_angles_deg(rot_idx));
            end

            if opts.enable_ref_guided_phase_fit
                comp_trials{end+1} = make_ref_guided_fit_trial_local( ...
                    x_freq, ref_syms, rel_freq, valid_idxs, opts, ref_idx, rot_angles_deg(rot_idx));
            end

            if opts.enable_ref_segment_anchor_fit
                try
                    comp_trials{end+1} = make_ref_anchor_trial_local( ...
                        x_freq, ref_syms, rel_freq, valid_idxs, opts, ref_idx, rot_angles_deg(rot_idx));
                catch
                    % 锚点不足时跳过该参考/旋转组合。
                end
            end
        end
    end
end

if isempty(comp_trials)
    error('没有可用补偿候选。');
end

best_comp = comp_trials{1};
for k = 2:numel(comp_trials)
    cand = comp_trials{k};
    if isnan(best_comp.score_ber) || (~isnan(cand.score_ber) && cand.score_ber < best_comp.score_ber)
        best_comp = cand;
    end
end

[best_comp, carrier_info] = apply_carrier_shift_lock_local(best_comp, valid_idxs, opts);

trial = struct();
trial.offset = offset;
trial.sss_start_idx = sss_start;
trial.comp_mode = best_comp.mode;
trial.x_freq_corr = best_comp.x_corr;
trial.hex_candidates = best_comp.hex_candidates;
trial.demod_bits_candidates = best_comp.demod_bits_candidates;
trial.best_ber = best_comp.score_ber;
trial.best_rot_idx = best_comp.score_rot_idx;
trial.best_ref_idx = best_comp.score_ref_idx;
trial.best_err = best_comp.score_err;
trial.best_n = best_comp.score_n;
trial.best_rot_deg = rot_angles_deg(max(1, best_comp.score_rot_idx));
trial.carrier_shift = carrier_info.selected_shift;
trial.carrier_shift_raw_ber = carrier_info.raw_ber;
trial.carrier_shift_ber = carrier_info.shifted_ber;
trial.valid_count = numel(valid_idxs);
trial.valid_rel_freq_indices = rel_freq(valid_idxs);
end

function comp = make_blind_pow4_trial_local(x_freq, rel_freq, valid_idxs, freq_indices, adjacent_mask, opts)
syms_valid = x_freq(valid_idxs);
syms_pow4 = syms_valid.^4;
syms_pow4_unit = syms_pow4 ./ (abs(syms_pow4) + eps);

phase_steps = conj(syms_pow4_unit(1:end-1)) .* syms_pow4_unit(2:end);
phase_steps = phase_steps(adjacent_mask);
if numel(phase_steps) < 8
    error('连续有效载波不足，无法完成四次方相位斜率估计。');
end

slope4 = angle(sum(phase_steps));
phase0_4 = angle(sum(syms_pow4_unit .* exp(-1j * slope4 * freq_indices)));
phase_corr = (slope4 * rel_freq + phase0_4) / 4;
x_corr = x_freq(:) .* exp(-1j * phase_corr);
comp = score_comp_local('blind_pow4', x_corr, valid_idxs, opts);
end

function comp = make_ref_linear_trial_local( ...
    x_freq, ref_syms, rel_freq, valid_idxs, freq_indices, adjacent_mask, opts, ref_idx, rot_deg)
z = x_freq(valid_idxs) .* conj(ref_syms(valid_idxs));
z = z ./ (abs(z) + eps);

phase_steps = conj(z(1:end-1)) .* z(2:end);
phase_steps = phase_steps(adjacent_mask);
if numel(phase_steps) < 8
    error('连续有效载波不足，无法完成参考线性相位估计。');
end

slope = angle(sum(phase_steps));
phase0 = angle(sum(z .* exp(-1j * slope * freq_indices)));
phase_corr = slope * rel_freq + phase0;
x_corr = x_freq(:) .* exp(-1j * phase_corr);
mode = sprintf('ref_linear_ref%d_rot%d', ref_idx, rot_deg);
comp = score_comp_local(mode, x_corr, valid_idxs, opts);
end

function comp = make_ref_guided_fit_trial_local(x_freq, ref_syms, rel_freq, valid_idxs, opts, ref_idx, rot_deg)
[phase_corr, ~] = fit_ref_guided_phase_local( ...
    x_freq, ref_syms, rel_freq, valid_idxs, opts.ref_guided_phase_fit_order, ...
    opts.ref_guided_fit_max_iters, opts.ref_guided_fit_outlier_sigma);
x_corr = x_freq(:) .* exp(-1j * phase_corr);
mode = sprintf('ref_guided_fit_ref%d_rot%d_ord%d', ref_idx, rot_deg, opts.ref_guided_phase_fit_order);
comp = score_comp_local(mode, x_corr, valid_idxs, opts);
end

function comp = make_ref_anchor_trial_local(x_freq, ref_syms, rel_freq, valid_idxs, opts, ref_idx, rot_deg)
[phase_corr, ~] = fit_ref_segment_anchor_phase_local( ...
    x_freq, ref_syms, rel_freq, valid_idxs, opts.ref_anchor_seg_len, ...
    opts.ref_anchor_min_coh, opts.ref_anchor_fit_order, opts.ref_anchor_min_segments);
x_corr = x_freq(:) .* exp(-1j * phase_corr);
mode = sprintf('ref_anchor_ref%d_rot%d_seg%d_ord%d', ref_idx, rot_deg, ...
    opts.ref_anchor_seg_len, opts.ref_anchor_fit_order);
comp = score_comp_local(mode, x_corr, valid_idxs, opts);
end

function comp = score_comp_local(mode, x_corr, valid_idxs, opts)
% 当前实验解调所有 1024 个载波；valid_idxs 只参与补偿估计。
syms_payload = x_corr(:);
[hex_candidates, bits_candidates] = demod_qpsk_candidates_local(syms_payload, opts.decision_rotate_deg);
[ber, rot_idx, ref_idx, err, n] = score_bits_against_refs_local( ...
    bits_candidates, opts.ref_hex_list, opts.unstable_bit_positions, opts.unstable_nibble_positions);

comp = struct();
comp.mode = mode;
comp.x_corr = x_corr(:);
comp.valid_idxs = valid_idxs;
comp.hex_candidates = hex_candidates;
comp.demod_bits_candidates = bits_candidates;
comp.score_ber = ber;
comp.score_rot_idx = rot_idx;
comp.score_ref_idx = ref_idx;
comp.score_err = err;
comp.score_n = n;
end

function [comp_out, info] = apply_carrier_shift_lock_local(comp_in, valid_idxs, opts)
% 关键步骤：频域整体平移锁定。
% 如果 FFT 起点/频域索引出现整体错位，星座仍可能很好，但 bit 序列会整体错位。
% 这里仅尝试 x_corr 的 circular shift，随后重新判决并算 BER。
comp_out = comp_in;
info = struct('selected_shift', 0, 'raw_ber', comp_in.score_ber, 'shifted_ber', comp_in.score_ber);

if ~opts.enable_ref_carrier_shift_lock
    return;
end

shift_candidates = unique(round(opts.ref_carrier_shift_search(:).'));
if ~any(shift_candidates == 0)
    shift_candidates = unique([0, shift_candidates]);
end

best_shift = 0;
best_score = comp_in.score_ber;
best_comp = comp_in;

for sh = shift_candidates
    x_shift = circshift(comp_in.x_corr(:), sh);
    cand = score_comp_local(sprintf('%s+carrier_shift%+d', comp_in.mode, sh), ...
        x_shift, valid_idxs, opts);

    if isnan(best_score) || (~isnan(cand.score_ber) && cand.score_ber < best_score)
        best_score = cand.score_ber;
        best_shift = sh;
        best_comp = cand;
    end
end

accept_shift = best_shift == 0 || isnan(comp_in.score_ber) || ...
    (~isnan(best_score) && best_score <= comp_in.score_ber - opts.ref_carrier_shift_min_gain);

if accept_shift
    comp_out = best_comp;
    info.selected_shift = best_shift;
    info.shifted_ber = best_score;
end
end

function [phase_full, info] = fit_ref_guided_phase_local( ...
    x_freq, ref_syms, rel_freq, valid_idxs, fit_order, max_iters, outlier_sigma)
z = x_freq(valid_idxs) .* conj(ref_syms(valid_idxs));
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
    res = angle(exp(1j * (ph - ph_fit)));
    med_res = median(res(mask));
    mad_res = median(abs(res(mask) - med_res)) + eps;
    sigma = 1.4826 * mad_res;
    new_mask = abs(res - med_res) <= max(pi/10, outlier_sigma * sigma);
    if nnz(new_mask) <= fit_order + 2 || isequal(new_mask, mask)
        break;
    end
    mask = new_mask;
end

phase_full = polyval(p, rel_freq / f_scale);
phase_full = phase_full(:);
info = struct('n_used', nnz(mask));
end

function [phase_full, info] = fit_ref_segment_anchor_phase_local( ...
    x_freq, ref_syms, rel_freq, valid_idxs, seg_len, min_coh, fit_order, min_segments)
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

    z = x_freq(seg_idx) .* conj(ref_syms(seg_idx));
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

fit_order = min(fit_order, numel(anchors_f)-1);
p = weighted_polyfit_local(anchors_f / f_scale, anchors_phase, fit_order, anchors_coh);
phase_full = polyval(p, rel_freq / f_scale);
phase_full = phase_full(:);
info = struct('n_anchor', numel(anchors_f), 'mean_coh', mean(anchors_coh));
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
        n_use = min(length(bits_est), length(bits_ref));
        if n_use <= 0
            continue;
        end

        keep = true(n_use, 1);
        ub = unstable_bit_positions(:);
        ub = ub(ub >= 1 & ub <= n_use);
        keep(ub) = false;

        un = unstable_nibble_positions(:);
        un = un(un >= 1);
        for k = 1:numel(un)
            b1 = 4*(un(k)-1) + 1;
            b2 = min(4*un(k), n_use);
            if b1 <= n_use
                keep(b1:b2) = false;
            end
        end

        n_keep = nnz(keep);
        if n_keep <= 0
            continue;
        end

        err = nnz((bits_est(1:n_use) ~= bits_ref(1:n_use)) & keep);
        ber = err / n_keep;
        if isnan(best_ber) || ber < best_ber
            best_ber = ber;
            best_rot_idx = i;
            best_ref_idx = j;
            best_err = err;
            best_n = n_keep;
        end
    end
end
end

function [hex_candidates, bits_candidates] = demod_qpsk_candidates_local(syms, decision_rotate_deg)
hex_candidates = cell(1, 4);
bits_candidates = cell(1, 4);
for r = 1:4
    s = syms(:) .* exp(1j * decision_rotate_deg * pi/180) .* exp(-1j * (r-1) * pi/2);
    bits_i = real(s) < 0;
    bits_q = imag(s) < 0;
    bits = zeros(length(s)*2, 1);
    bits(1:2:end) = bits_i;
    bits(2:2:end) = bits_q;
    bits_candidates{r} = bits;
    hex_candidates{r} = bits_to_hex_local(bits);
end
end

function syms = bits_to_qpsk_pre_rotation_local(bits, decision_rotate_deg, rot_idx)
bits = bits(:);
n_sym = floor(length(bits) / 2);
bits = bits(1:2*n_sym);
bits_i = bits(1:2:end);
bits_q = bits(2:2:end);
syms_rot = (1 - 2*double(bits_i)) + 1j * (1 - 2*double(bits_q));
syms_rot = syms_rot / sqrt(2);
syms = syms_rot .* exp(-1j * decision_rotate_deg * pi/180) .* exp(1j * (rot_idx-1) * pi/2);
syms = syms(:);
end

function bits = hex_to_bits_local(hex_str)
hex_str = upper(regexprep(hex_str, '\s+', ''));
hex_str = regexprep(hex_str, '[^0-9A-F]', '');
bits = zeros(4*length(hex_str), 1);
for i = 1:length(hex_str)
    v = hex2dec(hex_str(i));
    bits(4*(i-1)+1) = bitget(v, 4);
    bits(4*(i-1)+2) = bitget(v, 3);
    bits(4*(i-1)+3) = bitget(v, 2);
    bits(4*(i-1)+4) = bitget(v, 1);
end
end

function hex_str = bits_to_hex_local(bits)
bits = bits(:);
hex_str = '';
for k = 1:4:length(bits)
    chunk = bits(k:min(k+3, length(bits)));
    val = 0;
    for b = 1:length(chunk)
        val = val + chunk(b) * 2^(length(chunk)-b);
    end
    hex_str = [hex_str, dec2hex(val)]; %#ok<AGROW>
end
end

function rel_freq = rel_freq_local(N_fft)
rel_freq = zeros(N_fft, 1);
rel_freq(1:N_fft/2) = (0:N_fft/2 - 1).';
rel_freq(N_fft/2+1:N_fft) = (-N_fft/2:-1).';
end

function out = merge_opts_local(defaults, user_opts)
out = defaults;
if nargin < 2 || isempty(user_opts)
    return;
end
names = fieldnames(user_opts);
for k = 1:numel(names)
    out.(names{k}) = user_opts.(names{k});
end
end
