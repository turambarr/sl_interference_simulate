% 独立脚本：用循环前缀(CP)相关定位OFDM符号边界
% 已知：输入为纯数据（无文件头），int16小端序，I/Q交织
% 参数：子载波数 N=1024，循环前缀长度 Ng=192
%
% 输出：
% - symCpStart0 : 每个OFDM符号CP起点（0-based复采样点索引）
% - symDataStart0: 每个OFDM符号有效数据起点（0-based）= symCpStart0 + Ng

clear; clc;

%% 参数区（按需修改）
inFile = '20250912222305_part1_cut2.iq';

Fs = 409.6e6;   % Hz（仅用于可选显示）
N = 1024;       % 有效符号长度（FFT点数）
Ng = 192;       % CP长度
Lsym = N + Ng;  % 一个OFDM符号总长度

% 要搜索的范围（0-based复采样点索引，包含端点）；留空则用全文件
startSample0 = 0;
endSample0 = []; % 例如：startSample0+Lsym*200-1

% 归一化/去直流（通常有利）
normalizeToUnit = true; % int16/32768
removeMean = true;

% 从全局最大峰开始建立“符号栅格”，每隔Lsym在附近微调
refineRadius = round(Lsym/4); % 每个期望位置±refineRadius内找局部最大
maxSymbolsToReport = 5000;    % 防止打印太多

% 鲁棒性：真实信号里全局最大峰可能是杂峰。
% 这里尝试多个 top 峰作为“锚点”，用生成的 symRel 处 M 值均值打分，选最优锚点。
numAnchorsToTry = 20;
anchorMinSeparation = round(Lsym/2);

% 诊断：如果你想看某个符号的度量值（比如与 run_segment_similarity.m 的 idx 对齐）
debugSymIndex = 140;

%% 调用可复用估计函数
opts = struct();
opts.normalizeToUnit = normalizeToUnit;
opts.removeMean = removeMean;
opts.refineRadius = refineRadius;
opts.maxSymbolsToReport = maxSymbolsToReport;
opts.numAnchorsToTry = numAnchorsToTry;
opts.anchorMinSeparation = anchorMinSeparation;
opts.debugSymIndex = debugSymIndex;
opts.topK = 20;
opts.verbose = true;

out = estimate_ofdm_cp_locations(inFile, N, Ng, startSample0, endSample0, Fs, opts);
symCpStart0 = out.symCpStart0;
symDataStart0 = out.symDataStart0;

M = out.diag.M;
d0_rel = out.diag.d0_rel;
symRel = out.diag.symRel;

%% 可视化：度量曲线 + 标记符号起点
figure('Name', 'OFDM CP Correlation Metric');
plot(startSample0 + d0_rel, M, 'LineWidth', 1);
hold on;
plot(symCpStart0, M(symRel+1), 'r.', 'MarkerSize', 12);
grid on;
xlabel('Sample Index (0-based)');
ylabel('M(d)');
title(sprintf('CP metric (N=%d, Ng=%d, Lsym=%d)', N, Ng, Lsym));

zoom on;
pan on;
datacursormode on;
