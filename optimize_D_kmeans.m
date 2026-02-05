% optimize_D_kmeans_no_offset.m
% 用 KMeans 聚类评价星座聚拢程度，搜索最佳降采样倍率 D
% 方法：两阶段网格搜索（粗扫 -> 细扫）
% 评价：固定 offset=0，只优化 D（不再遍历 offset）
% 注意：若最佳采样相位不在 0，最优 D 可能会受影响

clear; clc; close all;

%% 1) 基础参数
inFile      = 'sigtest1.iq';
startSample = 15530-874;
readLen     = 874 * 8;

base_D      = 6.398;
range_width = 0.1;
bounds      = [base_D - range_width, base_D + range_width];

% 两阶段搜索
step_coarse = 1e-4;
step_fine   = 1e-5;
refine_half_width = 5e-4;

% 评价参数
K = 4;                         % 你的假设：星座应有 4 簇（QPSK/4QAM）
interpolation_method = 'pchip';% 避免 spline 过冲
offset = 0;                    % 固定 offset=0 （删除 offset 搜索）

% KMeans 参数（稳定/可复现）
rng_seed   = 0;
km_reps    = 10;
km_maxiter = 200;
km_start   = 'plus';

% 可选：固定用于聚类的采样点数，避免不同 D 导致 N 差异影响评分
use_fixed_N = true;
fixed_N = 1200;   % 你可以改，建议 >= 500

fprintf('==== KMeans-based D Search (No Offset) ====\n');
fprintf('Search bounds: [%.6f, %.6f]\n', bounds(1), bounds(2));
fprintf('Coarse step: %.1e, Fine step: %.1e\n', step_coarse, step_fine);
fprintf('Interp: %s, KMeans reps: %d, offset=%.2f\n', interpolation_method, km_reps, offset);

%% 2) 读取数据
[x_raw, ~] = iq_read_int16_le(inFile, startSample, readLen);
x_raw = double(x_raw);
x_raw = x_raw - mean(x_raw);
x_raw = x_raw / mean(abs(x_raw));
t_raw = 0:(length(x_raw)-1);

%% 3) 粗搜索
Ds_coarse = bounds(1):step_coarse:bounds(2);
nC = numel(Ds_coarse);

score_coarse = inf(nC,1);

fprintf('\n[Stage 1] Coarse grid search: %d points...\n', nC);
for i = 1:nC
    D = Ds_coarse(i);
    score_coarse(i) = evaluate_score_no_offset( ...
        D, offset, t_raw, x_raw, interpolation_method, ...
        K, rng_seed, km_reps, km_maxiter, km_start, ...
        use_fixed_N, fixed_N);

    fprintf('  %4d / %4d, D=%.6f, score=%.6f\n', i, nC, D, score_coarse(i));
end

[best_score_coarse, idxC] = min(score_coarse);
bestD_coarse = Ds_coarse(idxC);
fprintf('Coarse best: D=%.8f, score=%.8f\n', bestD_coarse, best_score_coarse);

%% 4) 细搜索（围绕粗最优）
a = max(bounds(1), bestD_coarse - refine_half_width);
b = min(bounds(2), bestD_coarse + refine_half_width);
Ds_fine = a:step_fine:b;
nF = numel(Ds_fine);

score_fine = inf(nF,1);

fprintf('\n[Stage 2] Fine grid search: %d points in [%.6f, %.6f]...\n', nF, a, b);
for i = 1:nF
    D = Ds_fine(i);
    score_fine(i) = evaluate_score_no_offset( ...
        D, offset, t_raw, x_raw, interpolation_method, ...
        K, rng_seed, km_reps, km_maxiter, km_start, ...
        use_fixed_N, fixed_N);

    if mod(i, max(1, floor(nF/20))) == 0
        fprintf('  %4d / %4d, D=%.6f, score=%.6f\n', i, nF, D, score_fine(i));
    end
end

[best_score, idxF] = min(score_fine);
bestD = Ds_fine(idxF);

fprintf('\n==== DONE ====\n');
fprintf('Best D: %.8f\n', bestD);
fprintf('Best score: %.8f\n', best_score);

%% 5) 画图验证
[ok, x_best] = resample_with_offset(bestD, offset, t_raw, x_raw, interpolation_method, use_fixed_N, fixed_N);
if ~ok
    warning('Resampling failed at bestD. Check data length / bounds.');
end

figure('Position',[100,100,1100,450], 'Name','D Optimization (No Offset)');
subplot(1,2,1);
plot(Ds_fine, score_fine, 'k-'); grid on;
xlabel('D'); ylabel('Score (cluster compactness)');
title('Fine-grid score curve');
xline(bestD, 'r--', sprintf('Best D=%.8f', bestD));

subplot(1,2,2);
plot(x_best, '.', 'MarkerSize', 4);
axis square; grid on;
title(sprintf('Constellation @ Best D=%.8f (offset=0)', bestD));
xlabel('I'); ylabel('Q');
xlim([-2 2]); ylim([-2 2]);

%% =================== 辅助函数 ===================

function score = evaluate_score_no_offset( ...
    D, offset, t_raw, x_raw, interp_method, ...
    K, rng_seed, km_reps, km_maxiter, km_start, ...
    use_fixed_N, fixed_N)

    [ok, x_res] = resample_with_offset(D, offset, t_raw, x_raw, interp_method, use_fixed_N, fixed_N);
    if ~ok
        score = 1e9;
        return;
    end

    X = [real(x_res(:)), imag(x_res(:))];

    if size(X,1) < 100
        score = 1e9;
        return;
    end

    rng(rng_seed);
    try
        [~, ~, sumd] = kmeans(X, K, ...
            'Start', km_start, ...
            'MaxIter', km_maxiter, ...
            'Replicates', km_reps, ...
            'Display', 'off');

        score = sum(sumd) / size(X,1);
    catch
        score = 1e9;
    end
end

function [ok, x_res] = resample_with_offset(D, offset, t_raw, x_raw, interp_method, use_fixed_N, fixed_N)
    t_new = offset : D : (length(x_raw)-1);

    if numel(t_new) < 50 || any(t_new < 0) || any(t_new > t_raw(end))
        ok = false; x_res = []; return;
    end

    if use_fixed_N && numel(t_new) > fixed_N
        t_new = t_new(1:fixed_N);
    end

    x_res = interp1(t_raw, x_raw, t_new, interp_method);

    if any(isnan(x_res))
        ok = false; x_res = []; return;
    end

    ok = true;
end
