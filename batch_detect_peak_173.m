% batch_detect_peak_173.m
% 遍历 sigtest1.iq ~ sigtest173.iq，按 autocorr_test2_slide128.m 的方法
% 计算滑动自相关实部曲线 M_real，并提取“正峰值”对应的 x 轴坐标。

clear; clc;

%% 参数区（按需修改）
file_prefix = 'sigtest';
file_suffix = '.iq';
file_range = 1:173;

startSample0 = 0;      % 0-based 起始读取点
numSamples = 50000;    % 每个文件读取点数；设 [] 表示读到文件尾

normalizeToUnit = true;
removeMean = true;

W = 874;               % 滑动积分窗口
D = 874;               % 延迟

peak_min = 0.5;        % 仅把 >=peak_min 的局部极大值视作“目标峰”候选

out_csv = 'peak_xcoord_173.csv';

%% 主循环
results = cell(numel(file_range), 5);
row = 0;

fprintf('Start batch peak detection...\n');
for k = 1:numel(file_range)
    idx = file_range(k);
    inFile = sprintf('%s%d%s', file_prefix, idx, file_suffix);

    if ~exist(inFile, 'file')
        fprintf('[%3d] %-16s | missing\n', idx, inFile);
        continue;
    end

    try
        nRead = local_num_samples(inFile, numSamples, startSample0);
        [x, meta] = iq_read_int16_le(inFile, startSample0, nRead);
        L = meta.numSamplesRead;

        if L <= D || L < W
            fprintf('[%3d] %-16s | too short (L=%d)\n', idx, inFile, L);
            continue;
        end

        x = double(x);
        if normalizeToUnit
            x = x / 32768;
        end
        if removeMean
            x = x - mean(x);
        end

        % === 与 autocorr_test2_slide128.m 一致的滑动自相关 ===
        rx_delayed = x(1+D:end);
        rx_base    = x(1:end-D);

        conj_prod = conj(rx_base) .* rx_delayed;
        P_metric = filter(ones(1, W), 1, conj_prod);

        R_base = filter(ones(1, W), 1, abs(rx_base).^2);
        R_delayed = filter(ones(1, W), 1, abs(rx_delayed).^2);

        M_complex = P_metric ./ (sqrt(R_base .* R_delayed) + 1e-10);
        M_real = real(M_complex);
        M_real_full = [M_real; zeros(D,1)];

        t_axis = (0:length(x)-1).' + startSample0;

        % === 检测“正峰”x坐标：优先找 >=peak_min 的局部极大值 ===
        [peak_x, peak_y, peak_idx] = local_pick_target_peak(t_axis, M_real_full, peak_min);

        row = row + 1;
        results{row,1} = idx;
        results{row,2} = inFile;
        results{row,3} = peak_x;
        results{row,4} = peak_y;
        results{row,5} = peak_idx;

        fprintf('[%3d] %-16s | peak_x=%8d | peak_y=%+.4f\n', idx, inFile, peak_x, peak_y);

    catch ME
        fprintf('[%3d] %-16s | error: %s\n', idx, inFile, ME.message);
    end
end

%% 导出结果
results = results(1:row, :);
if row == 0
    warning('没有得到任何有效结果。');
    return;
end

T = cell2table(results, 'VariableNames', ...
    {'file_index','file_name','peak_x_sample','peak_value','peak_idx_in_vector'});

writetable(T, out_csv);
fprintf('\nDone. Results saved to %s\n', out_csv);

% 可选：在命令窗显示前 20 行
disp(T(1:min(20,height(T)), :));

%% ===== local helpers =====
function n = local_num_samples(inFile, numSamples, startSample0)
if ~isempty(numSamples)
    n = numSamples;
    return;
end
bytesPerComplexSample = 4;
info = dir(inFile);
if isempty(info)
    error('找不到文件: %s', inFile);
end
Ntotal = floor(info.bytes / bytesPerComplexSample);
if startSample0 >= Ntotal
    n = 0;
else
    n = Ntotal - startSample0;
end
end

function [peak_x, peak_y, peak_idx] = local_pick_target_peak(t_axis, y, peak_min)
% 手写局部极大值（避免依赖特殊工具箱）
n = length(y);
if n < 3
    [peak_y, peak_idx] = max(y);
    peak_x = t_axis(peak_idx);
    return;
end

cand = find(y(2:end-1) > y(1:end-2) & y(2:end-1) >= y(3:end)) + 1;

% 优先目标：大于阈值的正峰
cand_strong = cand(y(cand) >= peak_min);
if ~isempty(cand_strong)
    [~, iBest] = max(y(cand_strong));
    peak_idx = cand_strong(iBest);
else
    % 回退：全局最大值
    [~, peak_idx] = max(y);
end

peak_y = y(peak_idx);
peak_x = t_axis(peak_idx);
end
