% plot_full_signal_energy.m
% 读整段信号的能量并绘制时域包络和频域能量分布

clear; clc; close all;

inFile = 'sigtest1.iq';
fs_source = 409.6e6;

try
    d = dir(inFile);
    full_len = d.bytes / 4; % int16 complex，每个IQ点占4字节
    
    fprintf('Loading whole file: %s (length: %d)...\n', inFile, full_len);
    [x_full, ~] = iq_read_int16_le(inFile, 0, full_len); % 0-based start
    x_full = double(x_full);
    
    % 计算每一时刻的能量 (复数模的平方)
    signal_power_time = abs(x_full).^2; 
    
    % 1. 时域包络图 (幅度/能量)
    figure('Position', [100, 500, 1200, 400], 'Name', 'Whole File Time Domain Power');
    plot(signal_power_time, 'LineWidth', 1);
    title('全文件时域能量包络 (|x|^2)');
    xlabel('Sample Index (MATLAB 1-based)'); 
    ylabel('Power');
    grid on;
    
    % 可选：高亮你在 sss_demodulation 里面的读取区域
    startSample = 16386+874*6+218; 
    readLen = 7000;
    hold on;
    y_lim = ylim;
    high_x = [startSample+1, startSample+readLen, startSample+readLen, startSample+1];
    high_y = [y_lim(1), y_lim(1), y_lim(2), y_lim(2)];
    patch(high_x, high_y, 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    legend('Signal Power', 'SSS Processing Window (sss\_demodulation)');
    
    % 2. 频域能量频谱图 (Power Spectral Density)
    fprintf('Calculating Frequency Domain Power...\n');
    figure('Position', [100, 50, 1200, 400], 'Name', 'Whole File Frequency Spectrum');
    
    % 使用 Welch 方法计算全局功率谱密度，使图像更平滑
    nfft = 8192;
    window_pts = hamming(nfft);
    [pxx, f] = pwelch(x_full, window_pts, round(nfft/2), nfft, fs_source, 'centered');
    
    % 防止出现 -Inf 导致无法画图
    power_dB = 10 * log10(pxx + 1e-15);
    
    plot(f / 1e6, power_dB, 'b-', 'LineWidth', 1.5);
    grid on;
    title('全文件频域能量谱分布 (Welch PSD)');
    xlabel('Frequency (MHz)');
    ylabel('Power Spectral Density (dB/Hz)');
    
    % 根据纵坐标情况限制一下范围，排除过低的极点
    max_pow = max(power_dB);
    ylim([max_pow - 60, max_pow + 5]);
    
    fprintf('绘制完成。\n');
catch ME
    warning('Failed to process: %s', ME.message);
end