```matlab
% convert_iq_to_wv.m
% 独立脚本：将 .iq 或 .dat 文件（纯 int16 I/Q 交织裸数据）转换为 R&S .wv 波形文件
% 专供 Rohde & Schwarz 波形发生器使用

clear; clc;

%% 1. 参数设置
inFile  = '20250912222305_part1.iq';     % 输入的纯 IQ 数据文件
outFile = '20250912222305_part1_conv.wv';  % 导出的 .wv 波形文件名
fs_hz   = 409.6e6;                       % 采样率 (Hz)

%% 2. 读取没有头的 IQ 裸数据 (int16交织)
fprintf('正在读取输入文件: %s...\n', inFile);
f_in = fopen(inFile, 'rb');
if f_in == -1
    error('无法打开输入文件: %s', inFile);
end
data = fread(f_in, inf, 'int16');
fclose(f_in);

if isempty(data) || mod(length(data), 2) ~= 0
    error('数据为空或点数不是偶数，请检查文件格式是否为交织的 I/Q (int16)。');
end

%% 3. 规一化并求功率参数
% 对于有符号16位整数(int16)，全量程(Full Scale)为 32767
fprintf('正在计算信号的 RMS 和 PEAK 功率特性...\n');
x = data(1:2:end) + 1j * data(2:2:end);
x_norm = double(x) / 32767.0; 

rms_val = sqrt(mean(abs(x_norm).^2));
peak_val = max(abs(x_norm));

% R&S 仪表需要以 dB 形式记录 RMS 与最大偏离，方便其进行 DAC 功率校准
if rms_val > 0
    rms_db = 20 * log10(rms_val);
else
    rms_db = -100;
end

if peak_val > 0
    peak_db = 20 * log10(peak_val);
else
    peak_db = -100;
end

%% 4. 生成 R&S 合法的 ASCII 头信息
fprintf('正在构造 .wv 文件头...\n');
f_out = fopen(outFile, 'wb');
if f_out == -1
    error('无法创建输出文件: %s', outFile);
end

date_str = datestr(now, 'yyyy-mm-dd;HH:MM:SS');
num_complex_samples = length(x);

% 组装 Header 字符串
header = sprintf('{TYPE: SMU-WV, 0}{COMMENT: Converted from IQ by MATLAB}{DATE: %s}{CLOCK: %d}{SAMPLES: %d}{RMS: %.5f}{PEAK: %.5f}', ...
                 date_str, fs_hz, num_complex_samples, rms_db, peak_db);

% WAVEFORM 指示标记：后面的数字是 "二进制总字节数 + 1（因为带个'#'）"
data_bytes = length(data) * 2; % 每个 int16 是 2 个字节
waveform_tag = sprintf('{WAVEFORM-%d: #', data_bytes + 1);

%% 5. 写入文件
% 写入文本头
fwrite(f_out, [header, waveform_tag], 'char');

% 写入交织的无头二进制数据
fprintf('正在写入波形负载数据 (约 %.2f MB)...\n', data_bytes / 1024 / 1024);
fwrite(f_out, data, 'int16');

% 写入结尾符号 '}'
fwrite(f_out, '}', 'char');
fclose(f_out);

fprintf('=======================================================\n');
fprintf('转换完毕！!\n');
fprintf('已生成文件: %s\n', outFile);
fprintf('配置参数 -> 采样率: %.1f MHz, RMS电平: %.2f dBFS, 峰值电平: %.2f dBFS\n', ...
        fs_hz/1e6, rms_db, peak_db);
fprintf('=======================================================\n');

```