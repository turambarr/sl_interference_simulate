% 独立脚本：标注 sigtest1.iq 中指定区域的“突起”结构
% 目标：在 36000-44000 样本区间内，寻找幅度高于背景噪声但低于主信号的突起
% 并绘图标注

clear; clc; close all;

%% 1. 参数设置
inFile = 'sigtest1.iq';
startSample0 = 0;
numSamples = []; % 读全部

% 关注的搜索区间 (大致范围，用于限定搜索)
searchRegion = [36000, 45000];

% 平滑窗口长度 (用于提取包络)
smoothWin = 256;

%% 2. 读取数据
[x, meta] = iq_read_int16_le(inFile, startSample0, local_num_samples(inFile, numSamples, startSample0));
L = meta.numSamplesRead;
if L == 0, error('未读取到数据。'); end
x = double(x);
x = x / 32768;      % 归一化
x = x - mean(x);    % 去直流

% 计算幅度包络
amp = abs(x);
t_axis = (0:length(x)-1) + startSample0;

% 平滑包络 (Moving Average)
amp_smooth = filter(ones(1, smoothWin)/smoothWin, 1, amp);
% 修正 filter 带来的相位延迟 (向左平移 win/2)
shift = round(smoothWin/2);
amp_smooth = [amp_smooth(shift+1:end); zeros(shift,1)];

%% 3. 区域检测
% 策略：
% 1. 取 searchRegion 内的数据
% 2. 取 searchRegion 右侧的数据作为 "Noise Reference" (例如 50000-60000)
% 3. 设定阈值

idxROI = find(t_axis >= searchRegion(1) & t_axis <= searchRegion(2));
ampROI = amp_smooth(idxROI);

noiseRegion = [50000, 60000]; % 假设这里是纯背景噪声
idxNoise = find(t_axis >= noiseRegion(1) & t_axis <= noiseRegion(2));
if isempty(idxNoise)
    % 如果文件没那么长，就取最后 1000 点
    idxNoise = (length(x)-1000):length(x);
end
meanNoise = mean(amp_smooth(idxNoise));
stdNoise  = std(amp_smooth(idxNoise));

% 设定检测阈值 (Lower Bound)
lowerThreshold = meanNoise + 2.0 * stdNoise;

% 设定过滤阈值 (Upper Bound) - 用于剔除左侧过高的主信号尾部
% 策略：参考 ROI 左侧一段信号的电平，取其与噪声底之间的某个比例位置作为上限
% 例如：认为是 "protrusion" 的信号不应超过主信号动态范围的 20%
leftRefIdx = max(1, idxROI(1)-1000) : idxROI(1);
meanHigh = mean(amp_smooth(leftRefIdx));
upperThreshold = meanNoise + 0.2 * (meanHigh - meanNoise);

fprintf('噪声基底: %.4f, 左侧高电平参考: %.4f\n', meanNoise, meanHigh);
fprintf('检测门限: Lower=%.4f, Upper=%.4f\n', lowerThreshold, upperThreshold);

% 在 ROI 中寻找: 高于底噪 & 低于主信号尾部
isBurst = ampROI > lowerThreshold & ampROI < upperThreshold;

% 简单的逻辑：找到最长的连续段，或者包围所有 true 的范围
if any(isBurst)
    % 寻找 ROI 内满足条件的起止点
    % 这里简单处理：取 ROI 内第一个和最后一个超过阈值的点作为边界
    % 也可以用 find_bursts 逻辑做更细致的分割
    
    burstIdxRel = find(isBurst);
    detectedStart = t_axis(idxROI(burstIdxRel(1)));
    detectedEnd   = t_axis(idxROI(burstIdxRel(end)));
    
    fprintf('检测到目标结构：\n');
    fprintf('  范围: %d -- %d (长度 %d)\n', detectedStart, detectedEnd, detectedEnd - detectedStart + 1);
    fprintf('  平均幅度: %.4f (vs 噪声 %.4f)\n', mean(ampROI(burstIdxRel)), meanNoise);
else
    fprintf('未在指定区间 %d-%d 内检测到明显突起 (阈值 %.4f)。\n', searchRegion(1), searchRegion(2), threshold);
    % 如果没检测到，就默认标注整个搜索区间，由用户判断
    detectedStart = searchRegion(1);
    detectedEnd   = searchRegion(2);
end

%% 4. 绘图展示
figure('Position', [100, 100, 1200, 600], 'Name', 'Highlight Protrusion');

plot(t_axis, amp, 'Color', [0.8 0.8 0.8], 'DisplayName', 'Raw Amplitude'); hold on;
plot(t_axis, amp_smooth, 'b', 'LineWidth', 1.5, 'DisplayName', 'Smoothed Envelope');

% 画出搜索区间边界
xline(searchRegion(1), 'g--', 'Search Start');
xline(searchRegion(2), 'g--', 'Search End');

% 画出噪声参考线和阈值
yline(meanNoise, 'k--', 'Noise Floor');
yline(lowerThreshold, 'r--', 'Lower Thresh');
yline(upperThreshold, 'r-.', 'Upper Thresh');

% 标注检测到的区域 (高亮背景)
yl = ylim;
patch([detectedStart detectedEnd detectedEnd detectedStart], ...
      [yl(1) yl(1) yl(2) yl(2)], ...
      'y', 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName', 'Detected Region');

text(detectedStart, yl(2)*0.9, 'Target Burst', 'Color', [0.8 0.5 0], 'FontWeight', 'bold');

title(['sigtest1.iq - 目标结构标注 (' num2str(detectedStart) '-' num2str(detectedEnd) ')']);
xlabel('Sample Index'); ylabel('Amplitude');
legend;
grid on;
xlim([max(0, searchRegion(1)-5000), min(length(x), searchRegion(2)+5000)]);

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
