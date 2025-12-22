% 独立脚本：按复采样点数裁切IQ文件并导出（输入/输出均为纯数据，无文件头）
% 已知：数据格式=int16小端序，I/Q交织

clear; clc;

%% 参数区（按需修改）
inFile = '20250912222305_part1.iq';
outFile = 'test4.iq';

startSample = 5218304;   % 从第几个“复采样点”开始裁切（0-based）
endSample = 5262335;   % 裁切到第几个“复采样点”（0-based，包含端点）

numSamples = endSample - startSample + 1;
if numSamples <= 0
	error('endSample 必须 >= startSample（且为包含端点的索引）。');
end

%% 裁切
iq_cut_by_samples(inFile, outFile, startSample, numSamples);

fprintf('裁切完成: %s\n', outFile);


function iq_cut_by_samples(inFile, outFile, startSample, numSamples, chunkSamples)
%IQ_CUT_BY_SAMPLES 按复采样点数裁切IQ文件。
%
% - 数据格式：int16 little-endian，I/Q交织
% - 输入/输出均为纯数据（无文件头）

if nargin < 5 || isempty(chunkSamples)
	chunkSamples = 2e6;
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

% 定位输入数据起点
startOffset = startSample * bytesPerComplexSample;
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
