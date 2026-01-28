% batch_cut_bursts.m
% 批量裁切 burst #13 到 #100，跳过前 12 个
% 读取 FILES.md 解析 burst 信息

clear; clc;

inputFile = '20250912222305_part1.iq';
headerBytes = 100;

% 检查输入文件是否存在
if ~exist(inputFile, 'file')
    error('Input file %s not found. Make sure you have the file with header.', inputFile);
end

% 由于 FILES.md 内容不完整，改读 burst_info.txt
fid = fopen('burst_info.txt', 'r');
if fid == -1
    error('Could not open burst_info.txt');
end

fprintf('Start processing bursts...\n');

while ~feof(fid)
    line = fgetl(fid);
    if ~ischar(line), break; end
    
    % Example line: burst #1: [303104..348159] (len=45056)
    tokens = regexp(line, 'burst #(\d+): \[(\d+)\.\.(\d+)\]', 'tokens');
    
    if ~isempty(tokens)
        idx = str2double(tokens{1}{1});
        startSample = str2double(tokens{1}{2});
        endSample = str2double(tokens{1}{3});
        
        % 用户要求截取 13-100 (1-12已截取)
        if idx >= 101 && idx <= 173
             outFile = sprintf('sigtest%d.iq', idx);
             
             numSamples = endSample - startSample + 1;
             
             % 必须考虑文件头偏移，否则会错位 25 个样本 (100 bytes / 4 = 25 samples)
             % 这里的逻辑严格等效于：在无头文件上运行 cut_iq_by_samples
             fprintf('Cutting burst #%d: [%d..%d] len=%d -> %s (skip header=%d)\n', ...
                 idx, startSample, endSample, numSamples, outFile, headerBytes);
             
             try
                do_cut(inputFile, outFile, startSample, numSamples, headerBytes);
             catch ME
                 fprintf('Error cutting burst #%d: %s\n', idx, ME.message);
             end
        elseif idx > 100
             % 超过100就不处理了
             break; 
        end
    end
end
fclose(fid);
fprintf('All done.\n');

function do_cut(inFile, outFile, startSample, numSamples, headerBytes)
    % 每个 sample 4 字节 (int16 I, int16 Q)
    bytesPerSample = 4;
    
    fin = fopen(inFile, 'rb');
    if fin == -1
        error('Cannot open input file');
    end
    cleanupIn = onCleanup(@() fclose(fin));
    
    % 跳过文件头 + 跳到 startSample
    seekPos = headerBytes + startSample * bytesPerSample;
    fseek(fin, seekPos, 'bof');
    
    % 读取数据 (numSamples 个复数 -> 2 * numSamples 个 int16)
    data = fread(fin, numSamples * 2, 'int16');
    
    if length(data) ~= numSamples * 2
        warning('Read fewer samples than expected. Wanted %d, got %d', numSamples, length(data)/2);
    end
    
    % 写入输出文件
    fout = fopen(outFile, 'wb');
    if fout == -1
        error('Cannot open output file');
    end
    cleanupOut = onCleanup(@() fclose(fout));
    
    fwrite(fout, data, 'int16');
end
