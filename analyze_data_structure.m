% analyze_data_structure.m
% 专门用于探测IQ数据中是否存在异常的零填充或未使用的位
% 帮助判断是否存在 I,Q,0,0 或 12-bit packed in 16-bit 等情况

clear; clc;

filename = '20250912222305_part1.iq';
fileHeaderSize = 100;

% 读取 10MB 数据用于统计
bytesToRead = 1000 * 1024 * 1024; 

fid = fopen(filename, 'rb');
if fid == -1
    error('无法打开文件');
end

fseek(fid, fileHeaderSize, 'bof');
raw_bytes = fread(fid, bytesToRead, 'uint8=>uint8');
fclose(fid);

fprintf('===== 数据结构深度分析 =====\n');
fprintf('分析数据量: %.2f MB\n', numel(raw_bytes)/1024/1024);

%% 1. 转换为 int16 观察
% 假设是 Little Endian int16
raw_int16 = typecast(raw_bytes, 'int16');
numSamples = numel(raw_int16);

% 统计数值为0的比例
zeroCount = sum(raw_int16 == 0);
zeroRatio = zeroCount / numSamples;

fprintf('\n>>> 零值统计 (Zero Value Statistics)\n');
fprintf('总int16样本数: %d\n', numSamples);
fprintf('数值完全为0的样本数: %d\n', zeroCount);
fprintf('零值比例: %.2f%%\n', zeroRatio * 100);

if zeroRatio > 0.1
    fprintf('警告：零值比例很高，可能存在填充零（Padding）或信号间歇发射。\n');
else
    fprintf('提示：零值比例正常，不太可能是 I,Q,0,0 这种简单的填充结构。\n');
end

%% 2. 周期性零检测 (Periodic Zero Check)
% 检查是否每隔N个点出现一次零（针对 I,Q,0 结构的探测）
max_period = 8;
fprintf('\n>>> 周期性零检测 (Periodic Zero Check, Period 2-%d)\n', max_period);
found_period = false;
for p = 2:max_period
    % Reshape数据为 p 行
    % 截断末尾以整除
    n_cols = floor(numSamples / p);
    matrix = reshape(raw_int16(1:n_cols*p), p, n_cols);
    
    % 统计每一行全是0的比例
    for row = 1:p
        row_zeros = sum(matrix(row, :) == 0);
        row_zero_ratio = row_zeros / n_cols;
        if row_zero_ratio > 0.9  % 如果某一位置超过90%是0
            fprintf('  发现周期规律：周期 P=%d，第 %d 个位置几乎全为0 (%.2f%%)\n', p, row, row_zero_ratio*100);
            found_period = true;
        end
    end
end
if ~found_period
    fprintf('  未发现明显的周期性零填充结构。\n');
end

%% 3. 有效位数检测 (Effective Bit Depth)
% 对所有正数、所有负数分别做 bitor，看哪些位用到了
% 注意：负数在补码中高位是1，所以主要看正数的分布或绝对值的分布
% 这里简单做所有数据的 bitwise OR
accumulated_bits = raw_int16(1);
for i = 2:numel(raw_int16)
    accumulated_bits = bitor(accumulated_bits, raw_int16(i));
end

fprintf('\n>>> 位深度检测 (Bitwise Analysis)\n');
fprintf('累积按位或结果 (HEX): 0x%04X\n', uint16(accumulated_bits));
fprintf('二进制: %s\n', dec2bin(typecast(accumulated_bits, 'uint16'), 16));

% 检查低位是否全0
if bitand(accumulated_bits, 15) == 0 % 0x000F
    fprintf('  -> 低4位似乎未被使用（可能是12-bit左对齐）\n');
elseif bitand(accumulated_bits, 3) == 0
    fprintf('  -> 低2位似乎未被使用（可能是14-bit左对齐）\n');
else
    fprintf('  -> 低位有翻转，看起来利用了完整的精度。\n');
end

% 检查高位（对于正数，看最大值）
max_val = max(abs(double(raw_int16)));
fprintf('最大绝对值幅值: %d (2^%.2f)\n', max_val, log2(max_val));
if max_val < 2048 % 2^11
    fprintf('  -> 幅值较小，未充满12位动态范围。\n');
elseif max_val < 32768
    fprintf('  -> 幅值正常，使用了16位符号整数范围。\n');
end
