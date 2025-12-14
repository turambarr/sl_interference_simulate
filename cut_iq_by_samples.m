% 独立脚本：按复采样点数裁切IQ文件并导出（输出为纯数据，无文件头）
% 已知：输入数据格式=int16小端序，I/Q交织，输入文件头=100字节

clear; clc;

%% 参数区（按需修改）
inFile = '20250912222305_part1.iq';
outFile = '20250912222305_part1_cut1.iq';

headerBytes = 100;
writeHeader = false; % false=输出纯数据（无文件头）

startSample = 500000;   % 从第几个“复采样点”开始裁切（0-based）
numSamples = 1700000;  % 裁切多少个“复采样点”

%% 裁切
iq_cut_by_samples(inFile, outFile, startSample, numSamples, headerBytes, [], writeHeader);

fprintf('裁切完成: %s\n', outFile);


function iq_cut_by_samples(inFile, outFile, startSample, numSamples, headerBytes, chunkSamples, writeHeader)
%IQ_CUT_BY_SAMPLES 按复采样点数裁切IQ文件。
%
% - 数据格式：int16 little-endian，I/Q交织
% - 输入文件带文件头（headerBytes，默认100），裁切时会跳过该文件头
% - writeHeader=true：复制输入文件头到输出文件；false：输出纯数据（无文件头）

if nargin < 6 || isempty(chunkSamples)
	chunkSamples = 2e6;
end
if nargin < 5 || isempty(headerBytes)
	headerBytes = 100;
end
if nargin < 7 || isempty(writeHeader)
	writeHeader = false;
end
if startSample < 0 || numSamples < 0
	error('startSample/numSamples 必须为非负');
end

bytesPerComplexSample = 4; % int16 I + int16 Q

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

% 读取输入文件头
fseek(fin, 0, 'bof');
header = fread(fin, headerBytes, 'uint8=>uint8');
if numel(header) ~= headerBytes
	error('文件头读取失败：期望%d字节，实际读取%d字节。', headerBytes, numel(header));
end

% 可选：写输出文件头
if writeHeader
	fwrite(fout, header, 'uint8');
end

% 定位输入数据起点（跳过文件头）
startOffset = headerBytes + startSample * bytesPerComplexSample;
status = fseek(fin, startOffset, 'bof');
if status ~= 0
	error('fseek失败：startOffset=%d', startOffset);
end

remaining = numSamples;
while remaining > 0
	nThis = min(remaining, chunkSamples);
	raw = fread(fin, nThis * 2, 'int16=>int16');
	if isempty(raw)
		warning('已到达文件尾，提前结束裁切。');
		break;
	end
	fwrite(fout, raw, 'int16');
	remaining = remaining - floor(numel(raw) / 2);

	if numel(raw) < nThis * 2
		warning('读取不足（文件尾），提前结束裁切。');
		break;
	end
end

end
