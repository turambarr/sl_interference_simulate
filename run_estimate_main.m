% 参数脚本：估计功能主入口（点击运行版）
% 说明：
% - 只包含“估计/检测”功能（OFDM CP 定位、帧间保护间隔检测）。
% - 不包含验证类（segment_head_tail_similarity / run_segment_similarity）。
% - 直接在本文件修改参数，然后点击运行。

clear; clc;

%% 选择功能
% 可选：'ofdm_cp' 或 'guard_intervals'
% action = 'guard_intervals';
action = 'ofdm_cp';
%% 通用输入
inFile = 'test4.iq';
startSample0 = 0;
endSample0 = []; % [] 表示到文件尾

%% OFDM CP 定位参数（action='ofdm_cp'）
Fs = 409.6e6;
N = 1024;
Ng = 192;

% 是否画出 M(d) 与符号起点标记（仅对 ofdm_cp 生效）
makePlotOfdmMetric = true;

%% 帧间保护间隔参数（action='guard_intervals'）
winLen = 4096;
hop = 1024;
smoothLen = 9;
% 保护间隔检测：如果“帧/空隙”主要体现为均值/直流分量变化，removeMean=true 可能会抹平对比度
removeMeanGuard = false;
thrMethod = 'mad'; % 'mad' 或 'quantile'
thrK = 8;
thrQuantile = 0.90;
minBurstLen = 20000;
minGuardLen = 20000;
chunkSamples = 2e6;
makePlotGuard = true;

%% 运行
switch lower(action)
    case 'ofdm_cp'
        out = estimate_main('ofdm_cp', ...
            'inFile', inFile, ...
            'N', N, 'Ng', Ng, ...
            'startSample0', startSample0, 'endSample0', endSample0, ...
            'Fs', Fs, ...
            'normalizeToUnit', true, 'removeMean', true, ...
            'verbose', true);

        % 便于后续脚本直接用
        symCpStart0 = out.symCpStart0; %
        symDataStart0 = out.symDataStart0; %

        if makePlotOfdmMetric
            M = out.diag.M;
            d0_rel = out.diag.d0_rel;
            symRel = out.diag.symRel;
            Lsym = out.diag.Lsym;

            figure('Name', 'OFDM CP Correlation Metric');
            plot(startSample0 + d0_rel, M, 'LineWidth', 1);
            hold on;
            plot(symCpStart0, M(symRel+1), 'r.', 'MarkerSize', 12);
            grid on;
            xlabel('Sample Index (0-based)');
            ylabel('M(d)');
            title(sprintf('CP metric (N=%d, Ng=%d, Lsym=%d)', N, Ng, Lsym));
            zoom on; pan on; datacursormode on;
        end

    case 'guard_intervals'
        out = estimate_main('guard_intervals', ...
            'inFile', inFile, ...
            'startSample0', startSample0, 'endSample0', endSample0, ...
            'normalizeToUnit', true, 'removeMean', removeMeanGuard, ...
            'winLen', winLen, 'hop', hop, 'smoothLen', smoothLen, ...
            'thrMethod', thrMethod, 'thrK', thrK, 'thrQuantile', thrQuantile, ...
            'minBurstLen', minBurstLen, 'minGuardLen', minGuardLen, ...
            'chunkSamples', chunkSamples, ...
            'makePlot', makePlotGuard, ...
            'verbose', true);

        % 便于后续脚本直接用
        framesStart0 = out.framesStart0; %#ok<NASGU>
        framesEnd0   = out.framesEnd0; %#ok<NASGU>
        guardsStart0 = out.guardsStart0; %#ok<NASGU>
        guardsEnd0   = out.guardsEnd0; %#ok<NASGU>

    otherwise
        error('未知 action=%s（仅支持 ofdm_cp / guard_intervals）', action);
end
