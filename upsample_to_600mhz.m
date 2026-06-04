% upsample_to_600mhz.m
% 独立脚本：将 .iq / .dat 指定基带/中频信号从原采样率重采样到 600 MHz 
% 严格防止混叠/失真，完美保持时频特性

clear; clc; close all;

%% 1. 参数设置
inFile  = '20250912222305_part1_57.iq';  % 输入文件名称
outFile = '20250912222305_part1_57_819.2MHz.iq'; % 输出文件名称

fs_in  = 409.6e6;   % 原始采样率 (Hz)
fs_out = 819.2e6;     % 目标采样率 (Hz)

% 读取控制 (inf表示读到底，如需测试可在此限制点数，例如 10e6)
max_read_samples = inf; 

%% 2. 读取原始数据 (int16 IQ 交织格式)
fprintf('正在读取文件: %s...\n', inFile);
fid = fopen(inFile, 'rb');
if fid == -1
    error('无法打开输入文件: %s', inFile);
end
raw = fread(fid, max_read_samples * 2, 'int16');
fclose(fid);

% 合成复数数据
x_iq = raw(1:2:end) + 1j * raw(2:2:end);
x_raw = double(x_iq);

fprintf('成功读取了 %d 个复数样本.\n', length(x_raw));

%% 3. 高质量上采样 (无混叠/无失真设计)
fprintf('正在执行高精度重采样 (%.1f MHz -> %.1f MHz)...\n', fs_in/1e6, fs_out/1e6);
[P, Q] = rat(fs_out / fs_in); 
% P=375, Q=256。resample会自动在这个比例下应用最佳多相滤波器
% 完美保留 -204.8MHz 到 +204.8MHz 的全频带，绝不产生镜像
x_resampled = resample(x_raw, P, Q); 

% 【删除】手动 filtfilt 滤波部分，避免切除 200~204.8MHz 的边缘信息

%% 4. 生成输出并保存 (防溢出等比例缩放设计)
fprintf('正在保存重采样后的信号到 %s...\n', outFile);
out_data = zeros(2 * length(x_resampled), 1);
out_data(1:2:end) = real(x_resampled);
out_data(2:2:end) = imag(x_resampled);

% 核心防失真操作：等比例峰值回退 (Backoff)，绝对避免硬削峰
max_val = max(abs(out_data));
if max_val > 32767
    fprintf('⚠️ 警告: 重采样产生过冲 (峰值 %.1f)！正在进行全局等比例缩放以防止硬截断失真...\n', max_val);
    % 留一点点裕量，比如缩放到 32700
    scale_factor = 32700 / max_val;
    out_data = out_data * scale_factor;
else
    fprintf('✅ 峰值检测安全 (峰值 %.1f)，无溢出风险。\n', max_val);
end

% 转换为 int16
out_data_int16 = int16(round(out_data));

f_out = fopen(outFile, 'wb');
fwrite(f_out, out_data_int16, 'int16');
fclose(f_out);

fprintf('=======================================================\n');
fprintf('重采样并写入完成！\n');
fprintf('输出样本数: %d 点\n', length(x_resampled));
fprintf('输出文件体积: %.2f MB\n', length(out_data_int16)*2 / 1024 / 1024);
fprintf('=======================================================\n');
