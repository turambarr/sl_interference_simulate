% check_ofdm_constellation.m
% 功能：读取指定位置的一个 OFDM 符号长度的数据，做 FFT，然后画出频域星座图
% 用途：检查子载波上的调制方式 (QPSK, 16QAM, etc.)
%
% 【重要提示】
% 为了看到清晰的星座图，必须满足以下苛刻条件：
% 1. start_sample 必须精确对应去掉了CP（循环前缀）后的 OFDM 符号起始位置。
%    (如果包含了CP或者位置偏了，星座图会旋转变成圆环)
% 2. n_fft 必须等于发送端的 FFT 长度。
% 3. 必须没有显著的频率偏移 (CFO)，否则也会转圈。

clear; clc; close all;

%% 1. 参数设置
filename = 'sigtest1.iq'; 

% 您需要手动尝试滑动这个 offset 来找到最佳对齐点
% 如果图是圆环，说明没对齐，请尝试 +/- 少量样点
start_sample = 15515-874+6992+162+324;    

n_fft = 1024;            % 假设的 FFT 长度 (OFDM子载波数)
                         % 常见值: 1024, 2048, 4096
                         
fs = 409600000;          % 采样率 (用于标记频率轴，可选)

%% 2. 读取数据 (仅读取一个 FFT 窗口的长度)
fprintf('正在读取 %s (Start=%d, Len=%d)...\n', filename, start_sample, n_fft);
[x, ~] = iq_read_int16_le(filename, start_sample, n_fft);

if length(x) < n_fft
    error('数据不足，无法进行FFT');
end

x = double(x);
% x = x - mean(x); % 去直流 (即去除 DC 子载波分量，可选)

%% 3. 核心处理: FFT
% 变换到频域
Y = fft(x, n_fft);

% 通常我们会做 fftshift，把 DC (0Hz) 移到中间
Y_shifted = fftshift(Y);

% 归一化 (为了画图好看，缩放到[-1, 1]附近)
Y_shifted = Y_shifted / mean(abs(Y_shifted)); 

%% 4. 绘图
figure('Position', [100, 100, 1200, 600], 'Name', 'OFDM Constellation Check');

% 子图1: 频域幅度谱 (辅助看哪些子载波是激活的)
subplot(1, 2, 1);
f_axis = linspace(-fs/2, fs/2, n_fft) / 1e6; % MHz
plot(f_axis, 20*log10(abs(Y_shifted)), 'b');
title('频域幅度谱');
xlabel('频率 (MHz)'); ylabel('幅度 (dB)');
grid on; axis tight;
subtitle('波峰处代表有效子载波，波谷是空闲子载波');

% 子图2: 频域星座图
subplot(1, 2, 2);
plot(real(Y_shifted), imag(Y_shifted), '.', 'MarkerSize', 8);
axis square; grid on;
title(sprintf('频域星座图 (Start=%d, FFT=%d)', start_sample, n_fft));
xlabel('In-Phase (I)'); ylabel('Quadrature (Q)');

% 画个单位圆作参考
hold on; 
th = 0:pi/50:2*pi;
plot(cos(th), sin(th), 'r:', 'LineWidth', 0.5); 
xlim([-2 2]); ylim([-2 2]);

fprintf('\n>>> 调试指南 <<<\n');
fprintf('1. 如果看到乱糟糟的一团：说明 start_sample 没对准 OFDM 符号的起始位置（还在 CP 里或跨符号了）。\n');
fprintf('2. 如果看到清晰的圆环：说明位置对得差不多了，但有【频偏】。\n');
fprintf('3. 如果看到 4个/16个/64个 清晰的簇：说明对齐完美！\n');
