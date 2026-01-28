% check_ofdm_psd.m
% 功能：读取指定文件的特定样点区间，绘制功率谱密度(PSD)
% 用途：通过观察频谱形状判断是否为 OFDM 信号
%
% OFDM 信号的典型频谱特征：
% 1. "矩形谱"：带内功率比较平坦（像一个门或者是平台）。
% 2.带外衰减快：边缘比较陡峭。
% 3. 如果是空闲时隙，则可能是底噪。
% 4. 如果是单载波(SC-FDMA/DFT-s-OFDM)，频谱形状可能会接近高斯或升余弦形状，与纯OFDM有所不同。

clear; clc; close all;

%% 1. 参数设置
filename = 'sigtest1.iq'; % 目标文件

% 选取样点区间
% Start=xxxx, Len=yyyy
start_sample = 15515-874+6992+162+324;      % 起始样点 (0-based)
read_len     = 6992;       % 读取长度 (建议2的幂次，如2048, 4096，便于FFT)

n_fft = 1024;              % FFT点数/Welch分段点数
fs = 409.6e6;              % 采样率 409.6 MHz

%% 2. 读取数据
% 假设 iq_read_int16_le 在路径中可用
fprintf('正在读取 %s (Start=%d, Len=%d)...\n', filename, start_sample, read_len);
[x, meta] = iq_read_int16_le(filename, start_sample, read_len);

if isempty(x)
    error('未读取到数据，请检查起始位置是否超出文件范围。');
end

x = double(x);

% 简单预处理
x = x - mean(x);          % 去直流
x = x / rms(x);           % 归一化功率

%% 3. 计算功率谱 (Welch法)
% 使用 Hamming 窗
window = hamming(n_fft);
n_overlap = n_fft / 2;

[pxx, f] = pwelch(x, window, n_overlap, n_fft, fs, 'centered');

% 转换为 dB
pxx_db = 10*log10(pxx);

%% 4. 绘图
figure('Position', [300, 300, 800, 600], 'Name', 'PSD Analysis');

% 时域图 (辅助确认是不是全是噪声)
subplot(2,1,1);
t_axis = (0:length(x)-1) / fs * 1e6; % 时间轴 (us)
plot(t_axis, abs(x));
title(sprintf('时域幅度 (Start=%d, Len=%d)', start_sample, read_len));
xlabel('时间 (\mus)'); ylabel('幅度');
grid on; axis tight;

% 频域图
subplot(2,1,2);
f_axis_mhz = f / 1e6; % 频率轴 (MHz)
plot(f_axis_mhz, pxx_db, 'LineWidth', 1.5, 'Color', 'b');
title('功率谱密度 (PSD Estimate)');
xlabel('频率 (MHz)');
ylabel('功率谱密度 (dBm/Hz 归一化)');
grid on; axis tight;


% 辅助线
y_max = max(pxx_db);
yline(y_max - 3, 'r--', 'Label', '-3dB');
yline(y_max - 10, 'g--', 'Label', '-10dB');

fprintf('\n>>> 分析指南 <<<\n');
fprintf('1. 观察下面子图的蓝色曲线形状。\n');
fprintf('2. 【OFDM特征】：呈现明显的“平台”状（矩形），顶部平坦，边缘陡峭。\n');
fprintf('3. 【单载波特征】：通常呈现中间高、两边低的山峰状。\n');
fprintf('4. 【噪声】：全频带杂乱无章，没有明显的凸起结构。\n');
