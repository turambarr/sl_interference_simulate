% 读取IQ文件的MATLAB脚本
% 文件格式：100字节文件头 + IQ采样数据
%
% IQ数据类型识别方法：
% 1. 魔数：文件头前2字节 = 0x6971 ('iq') 表示IQ文件
% 2. 文件头结构：通常包含采样率、时间戳、数据类型标志等
% 3. 文件大小反推：(总大小 - 100) / 字节数 = 样本数
%    - 如果能被8整除但不被16整除 → float32 (4+4字节)
%    - 如果能被16整除 → float64 (8+8字节)
%    - 如果能被4整除但不被8整除 → int16/uint16 (2+2字节)
% 4. 检查文件头中的类型标志位（通常在字节4-7或其他位置）

clear; clc;

% 文件路径
filename = '20250912222305_part1.iq';

% 打开文件
fid = fopen(filename, 'rb');
if fid == -1
    error(['无法打开文件: ' filename]);
end
cleanupObj = onCleanup(@() fclose(fid));

%% 读取文件头 (100字节)
fprintf('===== 读取文件头信息 =====\n\n');

% 重新定位到文件开始
frewind(fid);

% 读取整个文件头
header = fread(fid, 100, 'uint8=>uint8');
if numel(header) ~= 100
    error('文件头读取失败：期望100字节，实际读取%d字节。', numel(header));
end

fprintf('===== 文件头原始内容（100字节，HEX+ASCII）=====\n');
disp(hexdump_with_offsets(header));

fprintf('\n===== 按小端解析（Little-Endian）=====\n\n');

magic_bytes = header(1:2);
magic_str = char(magic_bytes.');
magic_u16 = u16le(header, 1);

flags_u16 = u16le(header, 3);
header_len_u32 = u32le(header, 5);
field_u64_08 = u64le(header, 9);

u32_16 = u32le(header, 17);
u32_20 = u32le(header, 21);
u32_24 = u32le(header, 25);
u32_28 = u32le(header, 29);

fprintf('偏移 0x00-0x01: Magic        = 0x%04X ("%s")\n', magic_u16, magic_str);
fprintf('偏移 0x02-0x03: Flags/Resv   = 0x%04X (%u)\n', flags_u16, flags_u16);
fprintf('偏移 0x04-0x07: HeaderLen    = 0x%08X (%u)\n', header_len_u32, header_len_u32);
if header_len_u32 ~= 100
    fprintf('  注意：HeaderLen与约定的100字节不一致，当前脚本仍按100字节读取。\n');
end

fprintf('偏移 0x08-0x0F: Field_u64_08 = 0x%016lX (%lu)\n', field_u64_08, field_u64_08);
fprintf('偏移 0x10-0x13: Field_u32_16 = 0x%08X (%u)  [推测：采样率 Hz]\n', u32_16, u32_16);
fprintf('偏移 0x14-0x17: Field_u32_20 = 0x%08X (%u)  [推测：频率/中频 Hz]\n', u32_20, u32_20);
fprintf('偏移 0x18-0x1B: Field_u32_24 = 0x%08X (%u)  [推测：频率/LO Hz]\n', u32_24, u32_24);
fprintf('偏移 0x1C-0x1F: Field_u32_28 = 0x%08X (%u)  [多为0：保留]\n', u32_28, u32_28);

nonzero_offsets = find(header ~= 0) - 1; % 转为0-based偏移
fprintf('\n非零字节位置（0-based偏移）: ');
if isempty(nonzero_offsets)
    fprintf('无（全为0）\n');
else
    fprintf('%s\n', sprintf('0x%02X ', nonzero_offsets));
end

if all(header(33:100) == 0)
    fprintf('偏移 0x20-0x63: 其余大部分为0（保留区）。\n');
else
    fprintf('偏移 0x20-0x63: 存在非零字节（请根据厂商格式进一步解释）。\n');
end

fprintf('\n===== 读取IQ数据 =====\n\n');

% 获取文件大小
fseek(fid, 0, 'eof');
fileSize = ftell(fid);
fprintf('文件总大小: %d 字节 (%.2f GB)\n', fileSize, fileSize/(1024^3));

% 计算数据部分大小（除去100字节文件头）
dataBytes = fileSize - 100;
fprintf('IQ数据部分大小: %d 字节\n', dataBytes);

fprintf('\n===== 数据格式自动检测 =====\n\n');

% 常见IQ数据格式及其字节大小
dataFormats = {
    'float32',  8,   '32位浮点 I/Q对 (4字节I + 4字节Q)';
    'float64',  16,  '64位浮点 I/Q对 (8字节I + 8字节Q)';
    'int16',    4,   '16位整数 I/Q对 (2字节I + 2字节Q)';
    'int32',    8,   '32位整数 I/Q对 (4字节I + 4字节Q)';
    'uint16',   4,   '16位无符号整数 I/Q对';
    'uint32',   8,   '32位无符号整数 I/Q对';
    'uint8',    2,   '8位无符号整数 I/Q对 (1字节I + 1字节Q)';
};

fprintf('可能的数据格式及对应样本数:\n');
fprintf('%-10s %-20s %s\n', '格式', '样本数', '字节大小');
fprintf('%-10s %-20s %s\n', '----', '------', '------');

for i = 1:size(dataFormats, 1)
    format_name = dataFormats{i, 1};
    bytes_per_sample = dataFormats{i, 2};
    format_desc = dataFormats{i, 3};
    
    numSamples_candidate = dataBytes / bytes_per_sample;
    
    % 检查样本数是否为整数
    if mod(numSamples_candidate, 1) == 0
        fprintf('%-10s %-20.0f %s\n', format_name, numSamples_candidate, format_desc);
    end
end

fprintf('根据数据大小，很多格式都可整除（本文件 dataBytes=%d 可被2/4/8/16整除）。\n', dataBytes);
fprintf('因此改用“探测数据分布”的方式推断（更可靠）。\n\n');

% 读取一小段数据用于探测（不依赖文件头里的类型标志位）
fseek(fid, 100, 'bof');
probe = fread(fid, 4096, 'uint8=>uint8');
[detectedPrecision, bytesPerSample, detectNote] = guess_iq_datatype(probe);
numSamples = dataBytes / bytesPerSample;

fprintf('探测结论: %s（%s）\n', detectedPrecision, detectNote);
fprintf('按该格式估计样本数: %.0f\n', numSamples);

% 首先定位到数据开始位置
fseek(fid, 100, 'bof');

% 读取一些样本来验证数据格式
fprintf('\n===== 前5个样本的IQ数据（按探测格式读取）=====\n');
samplesPerRead = 5;

switch detectedPrecision
    case 'int16'
        iqRaw = fread(fid, samplesPerRead * 2, 'int16=>double');
    case 'uint16'
        iqRaw = fread(fid, samplesPerRead * 2, 'uint16=>double');
    case 'int32'
        iqRaw = fread(fid, samplesPerRead * 2, 'int32=>double');
    case 'float32'
        iqRaw = fread(fid, samplesPerRead * 2, 'single=>double');
    case 'float64'
        iqRaw = fread(fid, samplesPerRead * 2, 'double=>double');
    otherwise
        error('不支持的探测格式: %s', detectedPrecision);
end

if numel(iqRaw) < samplesPerRead * 2
    error('数据读取不足：期望%d个数值，实际读取%d个。', samplesPerRead * 2, numel(iqRaw));
end

% 重新整理为I-Q对
iq_complex = complex(iqRaw(1:2:end), iqRaw(2:2:end));
for i = 1:samplesPerRead
    fprintf('样本 %d: I=%.6f, Q=%.6f, 幅度=%.6f\n', i, real(iq_complex(i)), imag(iq_complex(i)), abs(iq_complex(i)));
end

% 关闭文件（由onCleanup负责，这里不再手动关闭）

fprintf('\n===== 文件头信息总结 =====\n');
fprintf('文件类型: IQ采样数据\n');
fprintf('文件头大小: 100 字节\n');
fprintf('估计样本数: %.2e\n', numSamples);
fprintf('数据格式(探测): %s I/Q 对\n', detectedPrecision);


function v = u16le(b, oneBasedIndex)
v = typecast(uint8(b(oneBasedIndex:oneBasedIndex+1)), 'uint16');
end


function v = u32le(b, oneBasedIndex)
v = typecast(uint8(b(oneBasedIndex:oneBasedIndex+3)), 'uint32');
end


function v = u64le(b, oneBasedIndex)
v = typecast(uint8(b(oneBasedIndex:oneBasedIndex+7)), 'uint64');
end


function out = hex_dump(bytes)
bytes = uint8(bytes(:)).';
out = sprintf('%02X ', bytes);
end


function out = ascii_dump(bytes)
bytes = uint8(bytes(:)).';
chars = char(bytes);
mask = bytes < 32 | bytes > 126;
chars(mask) = '.';
out = chars;
end


function out = hexdump_with_offsets(bytes)
% 以 16字节/行 输出：偏移 + HEX + ASCII
bytes = uint8(bytes(:)).';
n = numel(bytes);
lines = strings(0);
for off = 0:16:(n-1)
    chunk = bytes(off+1:min(off+16, n));
    hexPart = sprintf('%02X ', chunk);
    if numel(chunk) < 16
        hexPart = [hexPart repmat('   ', 1, 16-numel(chunk))];
    end
    asciiPart = ascii_dump(chunk);
    lines(end+1) = sprintf('0x%04X: %s |%s|', off, hexPart, asciiPart);
end
out = strjoin(lines, newline);
end


function [precision, bytesPerSample, note] = guess_iq_datatype(probe)
% 通过数据分布探测IQ样本的基础类型（little-endian假设）。
% 返回precision为 'int16'/'float32'/...；bytesPerSample为每个复样本(I+Q)字节数。

probe = uint8(probe(:));

% float32探测：若大量NaN/Inf或出现接近3e38的极值，通常说明不是float32
nf4 = floor(numel(probe) / 4);
f4 = typecast(probe(1:4*nf4), 'single');
nanRatio = sum(isnan(f4)) / numel(f4);
infRatio = sum(isinf(f4)) / numel(f4);
finiteF4 = f4(isfinite(f4));
hugeRatio = 0;
if ~isempty(finiteF4)
    hugeRatio = sum(abs(double(finiteF4)) > 1e6) / numel(finiteF4);
end

% int16探测：统计幅度范围（若多数落在合理ADC范围内，倾向int16）
ni2 = floor(numel(probe) / 2);
i2 = typecast(probe(1:2*ni2), 'int16');
i2abs = abs(double(i2));
adcLikeRatio = sum(i2abs <= 4096) / numel(i2abs); % 经验阈值：很多采集卡有效位数<=12

if (nanRatio > 0.01) || (infRatio > 0) || (hugeRatio > 0.01)
    precision = 'int16';
    bytesPerSample = 4;
    note = sprintf('float32解释出现异常（NaN比例=%.3f, huge比例=%.3f），int16更像ADC计数(<=4096比例=%.3f)', nanRatio, hugeRatio, adcLikeRatio);
else
    % float32看起来“可能”合理，但仍优先看int16是否更像
    if adcLikeRatio > 0.95
        precision = 'int16';
        bytesPerSample = 4;
        note = sprintf('两者都可能，但int16更像ADC计数(<=4096比例=%.3f)', adcLikeRatio);
    else
        precision = 'float32';
        bytesPerSample = 8;
        note = sprintf('float32解释较正常（NaN比例=%.3f, huge比例=%.3f）', nanRatio, hugeRatio);
    end
end

end
