% 独立脚本：查看IQ信号时域图（仅 |IQ|）
% 已知：Fs=409.6MHz，数据格式=int16小端序，I/Q交织（纯数据文件，无文件头）

clear; clc;

%% 参数区（按需修改）
inFile = '20250912222305_part1.iq';
Fs = 409.6e6;      % 采样率 Hz

startSample = 0;   % 从第几个“复采样点”开始（0-based）
Nplot = 10000000;    % 绘图点数（复采样点数）

% ===== 时域“突起(burst)”检测（可选） =====
enableBurstDetect = true;

% 基于短时能量的包络检测（单位：复采样点）
detWinLen = 4096;
detHop = 1024;
detSmoothLen = 9;      % 以“bin”为单位
detThrMethod = 'mad';  % 'mad' 或 'quantile'
detThrK = 6;           % thr = median(E) + K*MAD(E)
detThrQuantile = 0.90; % quantile 方法

minBurstLen = 20000;   % 最短突起段长度（复采样点）
minGapLen = 20000;     % 合并相近突起段：gap 小于该值则合并（复采样点）

%% 绘图
[x, meta] = iq_read_int16_le(inFile, startSample, Nplot);
Nread = meta.numSamplesRead;
if Nread == 0
	error('未读取到数据。');
end

sampleIndex = startSample + (0:Nread-1);

figure('Name', sprintf('|IQ| Time Domain: %s', inFile));
hLine = plot(sampleIndex, abs(x));
grid on;
xlabel('Sample Index');
ylabel('|IQ| (int16 magnitude)');
title(sprintf('|IQ| (start=%d, N=%d)', startSample, Nread));

% 可选：检测并标注“突起(burst)”区间
if enableBurstDetect
	mag2 = abs(double(x)).^2;
	if numel(mag2) >= detWinLen
		acc = filter(ones(detWinLen, 1), 1, mag2);
		Evalid = acc(detWinLen:end) / detWinLen; % 窗末端从 (winLen-1) 到 end
		Ebins = Evalid(1:detHop:end);
		Ebins = local_smooth(Ebins, detSmoothLen);

		thr = local_threshold(Ebins, detThrMethod, detThrK, detThrQuantile);
		isActive = Ebins > thr;

		[runs1, ~] = local_runs(isActive);

		% 最短长度过滤（bin域）
		minBurstBins = ceil(minBurstLen / detHop);
		runs1 = runs1((runs1(:,2) - runs1(:,1) + 1) >= minBurstBins, :);

		% 合并间隔太短的突起段（bin域）
		if ~isempty(runs1)
			minGapBins = ceil(minGapLen / detHop);
			merged = runs1(1,:);
			for i = 2:size(runs1,1)
				gap = runs1(i,1) - merged(end,2) - 1;
				if gap <= minGapBins
					merged(end,2) = runs1(i,2);
				else
					merged = [merged; runs1(i,:)]; %#ok<AGROW>
				end
			end
			runs1 = merged;
		end

		% bin -> sampleIndex（0-based，全局索引）
		% Evalid(1) 对应窗末端 = startSample + (detWinLen-1)
		k0 = (0:numel(Ebins)-1).';
		binsEnd0 = startSample + (detWinLen - 1) + k0 * detHop;

		burstStart0 = zeros(0,1);
		burstEnd0 = zeros(0,1);
		if ~isempty(runs1)
			burstEnd0 = binsEnd0(runs1(:,2));
			burstStart0 = binsEnd0(runs1(:,1)) - (detWinLen - 1);

			% 裁剪到当前绘图区间
			burstStart0 = max(burstStart0, sampleIndex(1));
			burstEnd0 = min(burstEnd0, sampleIndex(end));
			keep = burstEnd0 >= burstStart0;
			burstStart0 = burstStart0(keep);
			burstEnd0 = burstEnd0(keep);
		end

		fprintf('=== Time-domain Burst Detection ===\n');
		fprintf('range=[%d..%d], winLen=%d, hop=%d, bins=%d\n', sampleIndex(1), sampleIndex(end), detWinLen, detHop, numel(Ebins));
		fprintf('thr=%.6g, bursts=%d\n', thr, numel(burstStart0));
		if ~isempty(burstStart0)
			for i = 1:numel(burstStart0)
				fprintf('burst #%d: [%d..%d] (len=%d)\n', i, burstStart0(i), burstEnd0(i), burstEnd0(i) - burstStart0(i) + 1);
			end
		end

		% 画突起段背景色块
		hold on;
		yl = ylim;
		for i = 1:numel(burstStart0)
			x1 = burstStart0(i);
			x2 = burstEnd0(i);
			patch([x1 x2 x2 x1], [yl(1) yl(1) yl(2) yl(2)], [0.85 0.90 1.0], ...
				'FaceAlpha', 0.18, 'EdgeColor', 'none');
		end
		uistack(hLine, 'top');
	else
		warning('数据长度(%d) < detWinLen(%d)，跳过 burst 检测。', numel(mag2), detWinLen);
	end
end

% 交互放大/查看
zoom on;
pan on;

% 点击查看坐标（数据光标）
datacursormode on;

function y = local_smooth(x, smoothLen)
if smoothLen <= 1
	y = x(:);
	return;
end
x = x(:);
w = ones(smoothLen,1) / smoothLen;
y = filter(w, 1, x);
if numel(y) >= smoothLen
	y(1:smoothLen-1) = x(1:smoothLen-1);
end
end

function thr = local_threshold(E, method, K, q)
E = E(:);
switch lower(method)
	case 'mad'
		med = median(E);
		m = median(abs(E - med));
		thr = med + K * m;
	case 'quantile'
		q = max(0, min(1, q));
		thr = local_quantile(E, q);
	otherwise
		error('detThrMethod 仅支持 mad / quantile');
end
thr = max(thr, eps);
end

function v = local_quantile(x, q)
x = sort(x(:));
n = numel(x);
if n == 0
	v = NaN;
	return;
end
if q <= 0
	v = x(1);
	return;
end
if q >= 1
	v = x(end);
	return;
end
pos = 1 + q * (n - 1);
lo = floor(pos);
hi = ceil(pos);
if lo == hi
	v = x(lo);
else
	frac = pos - lo;
	v = x(lo) * (1 - frac) + x(hi) * frac;
end
end

function [runs1, runs0] = local_runs(mask)
mask = logical(mask(:));
if isempty(mask)
	runs1 = zeros(0,2);
	runs0 = zeros(0,2);
	return;
end
d = diff([false; mask; false]);
starts1 = find(d == 1);
ends1 = find(d == -1) - 1;
runs1 = [starts1 ends1];

mask0 = ~mask;
d0 = diff([false; mask0; false]);
starts0 = find(d0 == 1);
ends0 = find(d0 == -1) - 1;
runs0 = [starts0 ends0];
end
