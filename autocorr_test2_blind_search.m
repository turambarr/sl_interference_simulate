% 独立脚本：对 iq 文件做“盲”滑动自相关结构检测
% 目标：遍历不同的 'window size' (同时假定 Delay = Window)，寻找显著的“反相-正相”特征
% 原理：PSS结构是 [..., -B, B, ...] 或类似，周期如果是 Unknown，就得扫。
% 这里我们扫描范围为 64 ~ 256 (或其他合理范围)，步进1。

clear; clc; close all;

%% 1. 参数设置
inFile = 'test2.iq';

startSample0 = 0;   
numSamples = 500000; % 为了速度，可以先读一部分 (比如 50万点)
                     % 如果是 []，可以读全文件，但遍历会慢

normalizeToUnit = true; 
removeMean = true;      

% 待扫描的长度范围 Range of Block Length (Blind Estimation)
scanRange = 870:1:880; % 修改为步长 1，范围可根据需要调整（例如 64~300 或 64~1024）
                     % 步长设为1以精确搜索未知长度

feature_threshold_neg = -0.6; % 负峰门限 (越负越像)
feature_threshold_pos = 0.6;  % 正峰门限 (越正越像)

%% 2. 读取数据
[x, meta] = iq_read_int16_le(inFile, startSample0, local_num_samples(inFile, numSamples, startSample0));
L = meta.numSamplesRead;
if L == 0, error('未读取到数据。'); end
x = double(x);
if normalizeToUnit, x = x / 32768; end
if removeMean, x = x - mean(x); end

fprintf('已读取信号长度: %d\n正在进行盲扫描 (W range: %d ~ %d)...\n', length(x), min(scanRange), max(scanRange));

%% 3. 遍历扫描 (Blind Search)
% 我们要寻找一个 W，使得 Metric = |Real(M)| 最大，或者使得 "负峰后紧接正峰" 最显著。
% 这里定义一个 Score： score(W) = max( 某个特征指标 )

best_W = 0;
best_score = 0;
best_M_curve = [];

scores = zeros(size(scanRange));

% 我们可以简单定义 Score 为 "全局最深的那个负峰" 的幅度
% 或者更严格："负峰深度 + 随后正峰的高度"
% 因为不知道具体结构，先找“存在显著负峰”的那个 W

for i = 1:length(scanRange)
    W_curr = scanRange(i);
    D_curr = W_curr; % 假设 Delay = Window Size (这是 Schmidl & Cox 的前提)
    
    if length(x) <= D_curr + W_curr
        continue; 
    end
    
    % ---- Core Calculation ----
    rx_delayed = x(1+D_curr:end);
    rx_base    = x(1:end-D_curr);
    
    conj_prod = conj(rx_base) .* rx_delayed;
    
    % 为了速度，这里可以用 conv 而不是 filter，或者手动优化
    % 但 filter 已经够快了
    P_metric = filter(ones(1, W_curr), 1, conj_prod);
    
    rx_power = abs(rx_base).^2;
    R_energy = filter(ones(1, W_curr), 1, rx_power);
    
    M_val = real( P_metric ./ (R_energy + 1e-10) );
    
    % ---- Scoring Strategy ----
    % 策略：寻找“最像反相结构”的 W
    % 如果结构是 [... -B, B, B ...]
    % 那么 M(n) 应该有一个接近 -1 的点，后面(隔W点)跟着接近 +1 的点
    
    % 1. 找到所有负峰 (低于 -0.5)
    min_val = min(M_val);
    
    % 2. 如果存在很深的负峰，我们再看它后面是不是有正峰
    % 为了简单稳健，我们直接用 abs(min_val) 作为基础分
    % 只有当 min_val < -0.4 时才认为可能是这种结构
    
    current_score = 0;
    if min_val < -0.4
         current_score = abs(min_val);
         % 可选：增加对 正峰 的检查
         % [min_v, min_idx] = min(M_val);
         % check_idx = min_idx + W_curr; % 理论上下一个块应该是正的
         % if check_idx < length(M_val) && M_val(check_idx) > 0.4
         %     current_score = current_score + M_val(check_idx);
         % end
    end
    
    scores(i) = current_score;
    
    if current_score > best_score
        best_score = current_score;
        best_W = W_curr;
        best_M_curve = M_val; % 暂存最佳曲线
        % best_M_curve 还需要补齐长度以便后续画图
        best_M_curve = [best_M_curve; zeros(D_curr, 1)]; 
    end
    
    if mod(i, 1) == 0 % 修改为每1次循环都输出，或者按需调整，例如每10次
        fprintf('Scanned W=%d, CurScore=%.3f, BestW=%d, BestScore=%.3f\n', W_curr, current_score, best_W, best_score);
    end
end

fprintf('扫描完成。Best W = %d, Score = %.3f\n', best_W, best_score);

if best_W == 0
    warning('未找到显著的反相结构（Score过低）。可能是阈值太高或信号中无此结构。显示得分最高的那个结果。');
    [~, max_idx] = max(scores);
    best_W = scanRange(max_idx);
    
    % Re-calculate for plotting
    W_curr = best_W; D_curr = best_W;
    rx_delayed = x(1+D_curr:end);
    rx_base    = x(1:end-D_curr);
    P_metric = filter(ones(1, W_curr), 1, conj_base .* rx_delayed); % Typo here? conj_prod logic
    conj_prod = conj(rx_base) .* rx_delayed;
    P_metric = filter(ones(1, W_curr), 1, conj_prod);
    rx_power = abs(rx_base).^2;
    R_energy = filter(ones(1, W_curr), 1, rx_power);
    best_M_curve = real( P_metric ./ (R_energy + 1e-10) );
    best_M_curve = [best_M_curve; zeros(D_curr, 1)];
end

%% 4. 结果展示
figure('Position', [100, 100, 1200, 800], 'Name', 'Blind Autocorrelation Search');

% 1. 评分分布图
subplot(3,1,1);
plot(scanRange, scores, '.-', 'LineWidth', 1);
xlabel('Window Size (W)'); ylabel('Inverted Structure Score');
title('盲扫描评分 (寻找最佳 W)');
grid on;
xline(best_W, 'r--', ['Best W=' num2str(best_W)]);

% 2. 信号时域图
t_axis = (0:length(x)-1) + startSample0;
subplot(3,1,2);
plot(t_axis, abs(x), 'Color', [0.6 0.6 0.6]);
ylabel('|x[n]|'); title(['原始信号 - ' inFile]);
grid on; axis tight;

% 3. 最佳 W 下的特征曲线
subplot(3,1,3);
plot(t_axis, best_M_curve, 'b', 'LineWidth', 1.5);
hold on;
yline(0, 'k-');
yline(-0.5, 'r:', 'Neg Thresh');
yline(0.5, 'r:', 'Pos Thresh');
ylabel('Real( M[n] )'); 
xlabel('Sample Index');
title(sprintf('最佳匹配 W=%d 时的自相关曲线 (负峰指示反相结构)', best_W));
grid on; axis tight;
ylim([-1.1, 1.1]);

% 自动缩放到负峰附近
[~, global_min_idx] = min(best_M_curve);
zoom_range = 2000; % 显示负峰左右各2000点
xlim([max(t_axis(1), t_axis(global_min_idx)-zoom_range), min(t_axis(end), t_axis(global_min_idx)+zoom_range)]);


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
