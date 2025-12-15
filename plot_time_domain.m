% 独立脚本：查看IQ信号时域图（仅 |IQ|）
% 已知：Fs=409.6MHz，数据格式=int16小端序，I/Q交织（纯数据文件，无文件头）

clear; clc;

%% 参数区（按需修改）
inFile = '20250912222305_part1.iq';
Fs = 409.6e6;      % 采样率 Hz

startSample = 0;   % 从第几个“复采样点”开始（0-based）
Nplot = 20000000;    % 绘图点数（复采样点数）

%% 绘图
[x, meta] = iq_read_int16_le(inFile, startSample, Nplot);
Nread = meta.numSamplesRead;
if Nread == 0
	error('未读取到数据。');
end

sampleIndex = startSample + (0:Nread-1);

figure('Name', sprintf('|IQ| Time Domain: %s', inFile));
plot(sampleIndex, abs(x));
grid on;
xlabel('Sample Index');
ylabel('|IQ| (int16 magnitude)');
title(sprintf('|IQ| (start=%d, N=%d)', startSample, Nread));

% 交互放大/查看
zoom on;
pan on;

% 点击查看坐标（数据光标）
datacursormode on;
