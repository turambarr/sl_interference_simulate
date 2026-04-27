% 把 20250912222305_part1.iq 的前 1/10 裁切出来
% 默认按：100字节文件头 + int16 I/Q交织（每复采样点4字节）处理

clear; clc;

inFile  = '20250912222305_part1.iq';
outFile = '20250912222305_part1_first10pct.iq';

headerBytes = 100;          % 你的 .iq 文件头长度
bytesPerSample = 4;         % int16 I + int16 Q
ratio = 0.1;                % 前十分之一

fin = fopen(inFile, 'rb');
if fin == -1
    error('无法打开输入文件: %s', inFile);
end
cleanupIn = onCleanup(@() fclose(fin));

fout = fopen(outFile, 'wb');
if fout == -1
    error('无法创建输出文件: %s', outFile);
end
cleanupOut = onCleanup(@() fclose(fout));

% 获取总文件大小
fseek(fin, 0, 'eof');
fileSize = ftell(fin);

if fileSize <= headerBytes
    error('文件大小异常：总大小(%d) <= headerBytes(%d)', fileSize, headerBytes);
end

dataBytes = fileSize - headerBytes;

% 仅保留完整复采样点
totalSamples = floor(dataBytes / bytesPerSample);
cutSamples = floor(totalSamples * ratio);
cutDataBytes = cutSamples * bytesPerSample;

if cutSamples <= 0
    error('可裁切样本数为0，请检查文件或参数。');
end

% 先拷贝文件头
fseek(fin, 0, 'bof');
header = fread(fin, headerBytes, 'uint8=>uint8');
if numel(header) ~= headerBytes
    error('文件头读取失败：期望 %d 字节，实际 %d 字节', headerBytes, numel(header));
end
fwrite(fout, header, 'uint8');

% 再拷贝前 1/10 数据区（分块拷贝，避免占内存）
fseek(fin, headerBytes, 'bof');
remaining = cutDataBytes;
chunkBytes = 16 * 1024 * 1024; % 16MB

while remaining > 0
    nRead = min(remaining, chunkBytes);
    buf = fread(fin, nRead, 'uint8=>uint8');
    if isempty(buf)
        warning('提前到达文件尾，实际写入可能少于预期。');
        break;
    end
    fwrite(fout, buf, 'uint8');
    remaining = remaining - numel(buf);
end

fprintf('裁切完成: %s\n', outFile);
fprintf('总复采样点: %d\n', totalSamples);
fprintf('裁切复采样点: %d (%.2f%%)\n', cutSamples, cutSamples / totalSamples * 100);
fprintf('输出文件大小: %d 字节\n', headerBytes + (cutDataBytes - remaining));
