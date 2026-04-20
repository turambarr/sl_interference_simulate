% match_voted_consensus.m
% 专门用于将 2-Pass 投票提炼出的唯一“众数”序列的 4 种相位变体
% 与 导频.txt (及其 16 种内部隐射变体) 进行终极比对。

clear; clc;

%% 1. 解析 导频.txt 以构建 16 种映射变体的庞大字典
pilot_file = '导频.txt';
if ~isfile(pilot_file), error('缺失 %s 文件', pilot_file); end
lines = readlines(pilot_file);
pilot_hex = cell(16, 1);
pilot_idx = 1;
for i = 1:length(lines)
    str = strtrim(lines(i));
    if startsWith(str, 'qp')
        parts = split(str, '=');
        if length(parts) == 2
            hex_str = regexprep(char(parts(2)), '[^0-9a-fA-F]', '');
            pilot_hex{pilot_idx} = hex_str;
            pilot_idx = pilot_idx + 1;
        end
    end
end

hex_dict = containers.Map(...
    {'0','1','2','3','4','5','6','7','8','9','A','B','C','D','E','F',...
     'a','b','c','d','e','f'}, ...
    {[0 0 0 0],[0 0 0 1],[0 0 1 0],[0 0 1 1],[0 1 0 0],[0 1 0 1],[0 1 1 0],[0 1 1 1],...
     [1 0 0 0],[1 0 0 1],[1 0 1 0],[1 0 1 1],[1 1 0 0],[1 1 0 1],[1 1 1 0],[1 1 1 1],...
     [1 0 1 0],[1 0 1 1],[1 1 0 0],[1 1 0 1],[1 1 1 0],[1 1 1 1]});

num_hex_chars = length(pilot_hex{1});
num_bits = num_hex_chars * 4;
num_symbols = floor(num_bits / 2);
ref_bits_matrix = zeros(32, num_symbols);
for k = 1:16
    hex_str = pilot_hex{k};  bits_k = zeros(1, length(hex_str)*4);
    for c = 1:length(hex_str)
        bits_k((c-1)*4+1 : c*4) = hex_dict(hex_str(c));
    end
    use_sym = min(floor(length(bits_k) / 2), num_symbols);
    
    % 改为单纯的 Interleaved [I,Q] 排列：
    % 即所有映射假设的基准数据构建时，奇数行存 I，偶数行存 Q
    ref_bits_matrix(2*k - 1, 1:use_sym) = bits_k(1:2:2*use_sym); % 奇数行为 I
    ref_bits_matrix(2*k, 1:use_sym)     = bits_k(2:2:2*use_sym); % 偶数行为 Q
end
num_symbols = use_sym;
ref_bits_matrix = ref_bits_matrix(:, 1:num_symbols);

% 只使用单纯的 [I,Q] 映射规则，不再自动搜索其他假说
% 因为前面的 ref_bits_matrix 已经严格按照 奇数行=I，偶数行=Q 构建好了
% 所以 idx 直接设为 1:32 得到的就是单纯匹配 iq 的基准数据！
maps = {};
maps{1} = struct('name', 'Pure [I,Q]', 'idx', 1:32);

num_maps = length(maps);
ref_matrices = cell(num_maps, 1);
for m = 1:num_maps
    idx = maps{m}.idx;
    ref_matrices{m} = ref_bits_matrix(idx, :); 
end

%% 2. 解析 voted_consensus_4phases.txt
voted_file = 'voted_consensus_4phases.txt';
if ~isfile(voted_file), error('缺失文件 %s', voted_file); end

voted_bits = zeros(4, 32);
lines = readlines(voted_file);
for i = 1:length(lines)
    str = char(strtrim(lines(i)));
    tokens = regexp(str, 'Phase (\d+).*:\s+([01]{32})', 'tokens');
    if ~isempty(tokens)
        p_idx = str2double(tokens{1}{1});
        bit_str = tokens{1}{2};
        voted_bits(p_idx, :) = bit_str - '0';
    end
end

%% 3. 终极比对
global_best_err = inf;
global_best_p = 0;
global_best_t = 0;
global_best_m = 0;
global_best_hyp = zeros(32, 1);

% 遍历记录排名前三的解，防止多种假设同样接近
results_list = [];

for p = 1:4
    hyp_bits = voted_bits(p, :).';
    for m = 1:num_maps
        M = ref_matrices{m};
        % 计算该相位假设在不同映射规则下，与时间轴上所有刻度的 Hamming距离
        err_vector = sum(M ~= hyp_bits, 1);
        [min_e, min_t] = min(err_vector);
        
        results_list = [results_list; p, m, min_t, min_e]; %#ok<AGROW>
        
        if min_e < global_best_err
            global_best_err = min_e;
            global_best_p = p;
            global_best_t = min_t;
            global_best_m = m;
            global_best_hyp = hyp_bits;
        end
    end
end

%% 4. 输出结论
fprintf('================ 投票合成序列的终极比对报告 ================\n\n');

if global_best_err > 12
    fprintf('【警告】绝对匹配误差非常高 (最佳也要错 %d/32)，这意味着 voted 序列未能汇聚到有效值。\n', global_best_err);
    fprintf('这可能是由于多个时间帧混杂导致投票结果变为无意义纯噪声。\n');
end

% 排序所有的检测结果
[~, sort_idx] = sort(results_list(:, 4), 'ascend');
top_results = results_list(sort_idx(1:min(5, size(results_list,1))), :);

fprintf('【全局最优匹配】：\n');
fprintf('最优匹配误差: %02d / 32 位\n', global_best_err);
fprintf('命中时间帧 T: %d\n', global_best_t);
fprintf('对应相位假说: Phase %d (%.0f°)\n', global_best_p, (global_best_p-1)*90);
fprintf('命中映射规则: %s\n\n', maps{global_best_m}.name);

fprintf('【最优匹配的逐位透视图】 (误差为 %d):\n', global_best_err);
M_best = ref_matrices{global_best_m};
ref_best_bits = M_best(:, global_best_t);
vis_str = char(zeros(1, 32));
for b = 1:32
    if ref_best_bits(b) == global_best_hyp(b)
        vis_str(b) = num2str(global_best_hyp(b)); % 匹配正确
    else
        vis_str(b) = 'X'; % 匹配错误
    end
end
fprintf('理想序列(T=%03d) : %s\n', global_best_t, sprintf('%d', ref_best_bits));
fprintf('投票序列透视图   : %s\n\n', vis_str);

fprintf('【排名前五的各种尝试】 (可能存在互相对称/等效关系):\n');
fprintf('排名 | 相位  | 时间 T | 误差/32 |  推断映射规则\n');
fprintf('------------------------------------------------------------------------\n');
for i = 1:size(top_results, 1)
    pr = top_results(i, 1);
    mr = top_results(i, 2);
    tr = top_results(i, 3);
    er = top_results(i, 4);
    fprintf(' #%d  | Ph %d | T=%03d  | %02d/32   |  %s\n', i, pr, tr, er, maps{mr}.name);
end
fprintf('\n');