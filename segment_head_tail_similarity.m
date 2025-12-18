function result = segment_head_tail_similarity(inFile, d0, offset, L, opts)
%SEGMENT_HEAD_TAIL_SIMILARITY 量化“任意段首尾相似程度”（适用于OFDM CP自证）
%
% result = segment_head_tail_similarity(inFile, d0, offset, L, opts)
%
% 输入（索引均为 0-based 复采样点索引）：
% - inFile : 纯IQ数据文件（int16 little-endian，I/Q交织）
% - d0     : 第一段起点（比如OFDM CP起点）
% - offset : 第二段相对第一段的偏移（OFDM时 offset = N）
% - L      : 对比长度（OFDM时 L = Ng）
% - opts   : 可选结构体
%            .normalizeToUnit (默认 true)  -> int16/32768
%            .removeMean      (默认 false) -> 是否去均值
%            .makePlot        (默认 true)  -> 是否画图
%
% 输出 result 结构体：
% - rho          : 复相关系数（归一化）
% - rhoMag       : |rho|（越接近1越相似）
% - rhoPhaseRad  : angle(rho)
% - ampRatio     : rms(tail)/rms(head)
% - nmse         : sum(|head-tail|^2)/sum(|tail|^2)
% - nmseDb       : 10log10(nmse)
% - d0, offset, L: 输入回显
% - head, tail   : 两段数据（复数向量）

if nargin < 5 || isempty(opts)
    opts = struct();
end
if ~isfield(opts, 'normalizeToUnit'); opts.normalizeToUnit = true; end
if ~isfield(opts, 'removeMean');      opts.removeMean = false; end
if ~isfield(opts, 'makePlot');        opts.makePlot = true; end
if ~isfield(opts, 'verbose');         opts.verbose = true; end

if d0 < 0 || offset < 0 || L <= 0
    error('d0/offset 必须非负，L 必须为正');
end

% 只读取需要的最小数据窗：[d0 .. d0+offset+L-1]
Nneed = offset + L;
[x, meta] = iq_read_int16_le(inFile, d0, Nneed);
if meta.numSamplesRead < Nneed
    error('读取不足：需要%d点，实际读取%d点（可能到文件尾）。', Nneed, meta.numSamplesRead);
end

x = double(x);
if opts.normalizeToUnit
    x = x / 32768;
end
if opts.removeMean
    x = x - mean(x);
end

head = x(1:L);
tail = x(offset+1:offset+L);

% 归一化复相关系数
num = head' * tail;
den = (norm(head) * norm(tail) + eps);
rho = num / den;

% NMSE
err = head - tail;
nmse = (sum(abs(err).^2)) / (sum(abs(tail).^2) + eps);

% 避免依赖工具箱的 rms()
tailMagRms = sqrt(mean(abs(tail).^2));
headMagRms = sqrt(mean(abs(head).^2));
ampRatio = tailMagRms / (headMagRms + eps);

result = struct();
result.d0 = d0;
result.offset = offset;
result.L = L;
result.rho = rho;
result.rhoMag = abs(rho);
result.rhoPhaseRad = angle(rho);
result.ampRatio = ampRatio;
result.nmse = nmse;
result.nmseDb = 10*log10(nmse + eps);
result.head = head;
result.tail = tail;

% 打印摘要
if opts.verbose
    fprintf('=== Head/Tail Similarity (0-based) ===\n');
    fprintf('d0=%d, offset=%d, L=%d\n', d0, offset, L);
    fprintf('|rho|=%.6f, angle(rho)=%.6f rad\n', result.rhoMag, result.rhoPhaseRad);
    fprintf('ampRatio=RMS(|tail|)/RMS(|head|)=%.6f\n', result.ampRatio);
    fprintf('NMSE=%.6e (%.2f dB)\n', result.nmse, result.nmseDb);
end

% 可视化：首尾对齐对比
if opts.makePlot
    n = 0:L-1;
    figure('Name', sprintf('Head/Tail Similarity @ d0=%d', d0));

    subplot(3,1,1);
    plot(n, real(head), 'b'); hold on;
    plot(n, real(tail), 'r--');
    grid on;
    xlabel('n'); ylabel('Real');
    title('Real(head) vs Real(tail)');
    legend('head', 'tail');

    subplot(3,1,2);
    plot(n, imag(head), 'b'); hold on;
    plot(n, imag(tail), 'r--');
    grid on;
    xlabel('n'); ylabel('Imag');
    title('Imag(head) vs Imag(tail)');
    legend('head', 'tail');

    subplot(3,1,3);
    plot(n, abs(head), 'b'); hold on;
    plot(n, abs(tail), 'r--');
    grid on;
    xlabel('n'); ylabel('Magnitude');
    title(sprintf('|head| vs |tail|   |rho|=%.4f, NMSE(dB)=%.2f', result.rhoMag, result.nmseDb));
    legend('head', 'tail');

    zoom on;
    pan on;
    datacursormode on;
end

end
