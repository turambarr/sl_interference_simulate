% 参数脚本：任意段“首尾相似程度”检测（CP自证）
% 说明：此脚本只放参数与调用，算法在 segment_head_tail_similarity.m

% 注意：如果你希望直接复用 locate_ofdm_symbols_cp.m 的输出 symCpStart0，
% 不要把它 clear 掉。
clearvars -except symCpStart0 symDataStart0;

%% 参数区（按需修改）
% 输入为纯数据（无文件头）：int16 小端序，I/Q交织
inFile = '20250912222305_part1_cut2.iq';

% OFDM已知参数
N = 1024;
Ng = 192;

% 选择一个起点 d0（0-based复采样点索引）
% 方式A：手工指定
% d0 = 123456;
clc;
% 方式B：如果你在同一个MATLAB会话里先运行了 locate_ofdm_symbols_cp.m
% 里面会生成 symCpStart0，则可以直接取其中一个
if exist('symCpStart0', 'var') && ~isempty(symCpStart0)
    idx = 150;
    d0 = symCpStart0(idx);
else
    d0 = 0; % 默认从0开始（你通常需要改成检测到的符号起点）
end

% 这里的“首尾”定义：对比 [d0 .. d0+Ng-1] 与 [d0+N .. d0+N+Ng-1]
offset = N;
L = Ng;

% 可选：快速排查参数是否匹配（不画图，只打印）
% - 如果你怀疑 N 或 Ng 不是 1024/192，把候选填进来
% - 留空 [] 表示不扫描
offsetCandidates = []; % 例如：[512 768 1024 1536 2048]
LCandidates = [];      % 例如：[128 160 192 224 256]

opts = struct();
opts.normalizeToUnit = true;
opts.removeMean = true;
opts.makePlot = true;
opts.verbose = true;

% 自检：在 d0 附近扫描，判断是否存在小偏移（索引/峰值微调问题）
% - 设为 0 表示不扫描
scanRadius = 500;   % 例如：50 或 200
scanStep = 1;     % 一般用1
scanMetric = 'rho'; % 'rho' 或 'nmseDb'

%% 调用
if scanRadius > 0
    % 先看中心点本身表现（不画图，避免干扰）
    centerOpts = opts;
    centerOpts.makePlot = false;
    centerOpts.verbose = false;
    rCenter = segment_head_tail_similarity(inFile, d0, offset, L, centerOpts);
    fprintf('=== Center Check ===\n');
    fprintf('d0=%d, |rho|=%.6f, NMSE(dB)=%.2f\n', d0, rCenter.rhoMag, rCenter.nmseDb);

    ds = (-scanRadius:scanStep:scanRadius);
    rhoMags = nan(size(ds));
    nmseDbs = nan(size(ds));

    scanOpts = opts;
    scanOpts.makePlot = false;
    scanOpts.verbose = false;

    for ii = 1:numel(ds)
        dTry = d0 + ds(ii);
        if dTry < 0
            continue;
        end
        r = segment_head_tail_similarity(inFile, dTry, offset, L, scanOpts);
        rhoMags(ii) = r.rhoMag;
        nmseDbs(ii) = r.nmseDb;
    end

    switch lower(scanMetric)
        case 'rho'
            [~, bestIdx] = max(rhoMags);
        case 'nmsedb'
            [~, bestIdx] = min(nmseDbs);
        otherwise
            error('scanMetric 只能是 ''rho'' 或 ''nmseDb''');
    end

    bestD0 = d0 + ds(bestIdx);
    fprintf('=== Local Scan ===\n');
    fprintf('center d0=%d, scanRadius=%d, step=%d\n', d0, scanRadius, scanStep);
    fprintf('best d0=%d (delta %+d), |rho|=%.6f, NMSE(dB)=%.2f\n', bestD0, ds(bestIdx), rhoMags(bestIdx), nmseDbs(bestIdx));

    % 打印 top 候选（便于看峰是否尖、是否多峰）
    topK = 10;
    switch lower(scanMetric)
        case 'rho'
            [~, order] = sort(rhoMags, 'descend');
        case 'nmsedb'
            [~, order] = sort(nmseDbs, 'ascend');
    end
    order = order(~isnan(rhoMags(order)) & ~isnan(nmseDbs(order)));
    order = order(1:min(topK, numel(order)));
    fprintf('=== Top Candidates ===\n');
    for jj = 1:numel(order)
        dCand = d0 + ds(order(jj));
        fprintf('%2d) d0=%d (delta %+d), |rho|=%.6f, NMSE(dB)=%.2f\n', jj, dCand, ds(order(jj)), rhoMags(order(jj)), nmseDbs(order(jj)));
    end

    % 用最佳点再跑一遍（带打印/画图）
    result = segment_head_tail_similarity(inFile, bestD0, offset, L, opts);

    % 额外：对 bestD0 做参数候选扫描
    if ~isempty(offsetCandidates) || ~isempty(LCandidates)
        if isempty(offsetCandidates); offsetCandidates = offset; end
        if isempty(LCandidates);      LCandidates = L; end

        fprintf('=== Param Sweep @ bestD0=%d ===\n', bestD0);
        sweepOpts = opts;
        sweepOpts.makePlot = false;
        sweepOpts.verbose = false;
        for oi = 1:numel(offsetCandidates)
            for li = 1:numel(LCandidates)
                off = offsetCandidates(oi);
                LL = LCandidates(li);
                r = segment_head_tail_similarity(inFile, bestD0, off, LL, sweepOpts);
                fprintf('offset=%-5d, L=%-4d  |rho|=%.6f  NMSE(dB)=%.2f\n', off, LL, r.rhoMag, r.nmseDb);
            end
        end
    end
else
    result = segment_head_tail_similarity(inFile, d0, offset, L, opts);
end
