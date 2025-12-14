function [x, meta] = iq_read_int16_le(filename, startSample, numSamples, headerBytes)
%IQ_READ_INT16_LE 读取int16小端序IQ文件（I/Q交织），返回复数向量。
%
% [x, meta] = iq_read_int16_le(filename, startSample, numSamples, headerBytes)
%
% - filename: .iq文件路径
% - startSample: 从第几个“复采样点”开始读（0-based）
% - numSamples: 读取多少个“复采样点”
% - headerBytes: 头大小（默认100）
%
% 数据格式：int16 little-endian，按 I0,Q0,I1,Q1,... 交织。

if nargin < 4 || isempty(headerBytes)
    headerBytes = 100;
end
if nargin < 3
    error('需要参数: filename, startSample, numSamples');
end
if startSample < 0 || numSamples < 0
    error('startSample/numSamples 必须为非负');
end

bytesPerComplexSample = 4; % int16 I + int16 Q

fid = fopen(filename, 'rb');
if fid == -1
    error('无法打开文件: %s', filename);
end
cleanupObj = onCleanup(@() fclose(fid));

% 读取文件头
fseek(fid, 0, 'bof');
header = fread(fid, headerBytes, 'uint8=>uint8');
if numel(header) ~= headerBytes
    error('文件头读取失败：期望%d字节，实际读取%d字节。', headerBytes, numel(header));
end

% 定位到数据区
dataOffset = headerBytes + startSample * bytesPerComplexSample;
status = fseek(fid, dataOffset, 'bof');
if status ~= 0
    error('fseek失败：dataOffset=%d', dataOffset);
end

raw = fread(fid, numSamples * 2, 'int16=>double');
if numel(raw) < numSamples * 2
    warning('实际读取到%d个int16值，少于期望%d（文件可能到尾）。', numel(raw), numSamples*2);
end

I = raw(1:2:end);
Q = raw(2:2:end);
N = min(numel(I), numel(Q));
I = I(1:N);
Q = Q(1:N);

x = complex(I, Q);

meta = struct();
meta.headerBytes = headerBytes;
meta.bytesPerComplexSample = bytesPerComplexSample;
meta.startSample = startSample;
meta.numSamplesRequested = numSamples;
meta.numSamplesRead = N;
meta.header = header;
end
