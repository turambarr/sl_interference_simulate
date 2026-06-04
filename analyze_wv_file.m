% analyze_wv_file.m
% 专门用于分析 .wv (通常是 Rohde & Schwarz 波形文件) 结构的脚本
% 并包含一个将 .iq 转换为 .wv 的示例功能

clear; clc;

wv_filename = '20250912222305_part1_57.wv';

if ~exist(wv_filename, 'file')
    fprintf('文件 %s 不存在，请检查路径。\n', wv_filename);
    return;
end

fprintf('===== 开始分析 .wv 文件: %s =====\n', wv_filename);

% 1. 分析文件头 (ASCII Header)
fid = fopen(wv_filename, 'rb');
header_str = '';
in_header = true;
binary_start_pos = 0;

% 逐字符读取，直到遇到特定的波形标签或二进数据标志
% R&S .wv 文件的 tag 形如 {TAG: VALUE}
while in_header && ~feof(fid)
    c = fread(fid, 1, 'uint8=>char');
    header_str = [header_str, c]; %#ok<AGROW>
    
    % 检查是否读到了波形数据的起始标志
    % 通常是 {WAVEFORM-xxxxx: #  其中 xxxxx 是字节数，# 后面紧跟纯二进制数据
    if endsWith(header_str, '#')
        % 检查前面是不是有 WAVEFORM 等字样
        if contains(header_str, 'WAVEFORM')
            in_header = false;
            binary_start_pos = ftell(fid);
        end
    end
    
    % 限制header读取的范围防止读完整个巨型文件
    if length(header_str) > 10000
        fprintf('警告: 超过10KB未找到 WAVERFORM 标签，可能不是标准 R&S wv 格式或无特殊头。\n');
        break;
    end
end

% 处理并打印头信息
fprintf('\n[检测到的 Headers / Metadata]:\n');
tags = regexp(header_str, '\{(.*?)\}', 'match');
for i = 1:length(tags)
    fprintf('  %s\n', tags{i});
end

fprintf('\n[文件结构参数]:\n');
fseek(fid, 0, 'eof');
file_size = ftell(fid);
fprintf('  总文件大小: %.2f MB\n', file_size / 1024 / 1024);

if binary_start_pos > 0
    binary_size = file_size - binary_start_pos;
    num_samples = binary_size / 4; % 假设 int16 I + int16 Q = 4 bytes
    fprintf('  Header 占用: %d bytes\n', binary_start_pos);
    fprintf('  Binary Data 占用: %d bytes\n', binary_size);
    fprintf('  推测复采样点数: %d\n', num_samples);
else
    fprintf('  未能明确切割 Header 和 Binary Data。\n');
end
fclose(fid);

%% 提供 .iq 转 .wv 的附加函数
fprintf('\n===================\n');
fprintf('附加说明：.iq 到 .wv 转换演示完毕（详见脚本源码内的 iq_to_wv 函数）。\n');

function iq_to_wv_example(iqFile, wvFile, fs_hz)
    % 这是一个通用的转换函数示例，将无头的 I/Q (int16) 转换为 R&S .wv 格式
    
    % 1. 读出所有 IQ 数据
    f_in = fopen(iqFile, 'rb');
    data = fread(f_in, inf, 'int16');
    fclose(f_in);
    
    % 组装复数以计算功率特性
    x = data(1:2:end) + 1j * data(2:2:end);
    
    % R&S .wv 通常使用有符号16位全量程 [-32768, 32767]
    % 归一化为小数以计算 RMS 和 PEAK
    % 如果原始iq已经是刚好不过载的int16，这里直接计算
    x_norm = double(x) / 32767.0; 
    
    rms_val = sqrt(mean(abs(x_norm).^2));
    peak_val = max(abs(x_norm));
    
    % R&S 要求写入时记录 RMS 和 PEAK 偏离（峰值到满量程的距离）等，
    % 以便信号源在重播时校准DAC增益。
    % 这里通常要求 RMS = 20*log10(rms_val) (参考全量程1) -> 负数
    rms_db = 20 * log10(rms_val);
    peak_db = 20 * log10(peak_val);
    
    % 2. 写文件
    f_out = fopen(wvFile, 'wb');
    
    % 创建标准的头信息
    date_str = datestr(now, 'yyyy-mm-dd;HH:MM:SS');
    
    header = '';
    header = [header, '{TYPE: SMU-WV, 0}'];               % 波形类型
    header = [header, '{COMMENT: Converted from IQ}'];    % 注释
    header = [header, sprintf('{DATE: %s}', date_str)];   % 日期
    header = [header, sprintf('{CLOCK: %d}', fs_hz)];     % 给定采样率
    header = [header, sprintf('{SAMPLES: %d}', length(x))]; % 样本数
    header = [header, sprintf('{RMS: %.5f}', rms_db)];    % RMS功率
    header = [header, sprintf('{PEAK: %.5f}', peak_db)];  % 偏置/峰值
    
    % 写入头的二进制长度标识
    data_bytes = length(data) * 2; % 每个 int16 是 2 字节
    waveform_tag = sprintf('{WAVEFORM-%d: #', data_bytes + 1); % +1是因为#算一个字节
    header = [header, waveform_tag];
    
    % 写入 ASCII Header
    fwrite(f_out, header, 'char');
    
    % 写入 二进制数据 (I/Q 必须交织，可能需要注意大小端，R&S通常接受Little Endian)
    fwrite(f_out, data, 'int16');
    
    % 写入 结尾符
    fwrite(f_out, '}', 'char');
    
    fclose(f_out);
    fprintf('已将 %s 转换为 %s (采样率: %d Hz)\n', iqFile, wvFile, fs_hz);
end
