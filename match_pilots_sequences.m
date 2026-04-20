% match_pilots_sequences.m
% 用于验证在 all_decoded_bits_4qam_4phases.txt 中解析出的 173 个文件的序列，
% 是否能和 导频.txt 中提供的序列相匹配。

clear; clc;

%% 1. 解析 导频.txt
pilot_file = '导频.txt';
if ~isfile(pilot_file)
    error('缺失 %s 文件', pilot_file);
end

lines = readlines(pilot_file);
pilot_hex = cell(16, 1);
pilot_idx = 1;

for i = 1:length(lines)
    str = strtrim(lines(i));
    if startsWith(str, 'qp')
        % 提取 = 后面的十六进制字符串并去掉空格
        parts = split(str, '=');
        if length(parts) == 2
            % 转为 char 数组并只保留有效的十六进制字符，消除 string 标量下标和隐藏换行符问题
            hex_str = regexprep(char(parts(2)), '[^0-9a-fA-F]', '');
            pilot_hex{pilot_idx} = hex_str;
            pilot_idx = pilot_idx + 1;
        end
    end
end

if pilot_idx <= 16
    error('导频.txt 中未找到完整的 16 个(488~495, 528~535)的导频子载波定义。');
end

% 将十六进制字符映射为二进制比特 (每个字符4个比特，左侧为MSB)
hex_dict = containers.Map(...
    {'0','1','2','3','4','5','6','7','8','9','A','B','C','D','E','F',...
     'a','b','c','d','e','f'}, ...
    {[0 0 0 0],[0 0 0 1],[0 0 1 0],[0 0 1 1],...
     [0 1 0 0],[0 1 0 1],[0 1 1 0],[0 1 1 1],...
     [1 0 0 0],[1 0 0 1],[1 0 1 0],[1 0 1 1],...
     [1 1 0 0],[1 1 0 1],[1 1 1 0],[1 1 1 1],...
     [1 0 1 0],[1 0 1 1],[1 1 0 0],[1 1 0 1],[1 1 1 0],[1 1 1 1]});

% 计算可以构成的 QPSK 符号个数 (每 2 个 bit 算作 1 个符号)
num_hex_chars = length(pilot_hex{1});
num_bits = num_hex_chars * 4;
num_symbols = floor(num_bits / 2);

% 构建所有时间点上理想的 32bit 观测矩阵
% Size: 32 (bits/symbol) x num_symbols
ref_bits_matrix = zeros(32, num_symbols);

for k = 1:16
    hex_str = pilot_hex{k};
    bits_k = zeros(1, length(hex_str)*4);
    for c = 1:length(hex_str)
        bits_k((c-1)*4+1 : c*4) = hex_dict(hex_str(c));
    end
    
    cur_symbols = floor(length(bits_k) / 2);
    use_sym = min(cur_symbols, num_symbols);
    
    % 分配给当前子载波 k: 行 2k-1 放 I比特，行 2k 放 Q比特
    ref_bits_matrix(2*k - 1, 1:use_sym) = bits_k(1:2:2*use_sym);
    ref_bits_matrix(2*k, 1:use_sym)     = bits_k(2:2:2*use_sym);
end

% 截断多余长度并保证一致
num_symbols = use_sym;
ref_bits_matrix = ref_bits_matrix(:, 1:num_symbols);
fprintf('成功从 导频.txt 加载 16 个子载波序列，等效时间长度 T=%d 个 QPSK 符号。\n\n', num_symbols);

%% 1.5 扩展各种映射假说
% 实际通信系统中可能存在各种 I/Q 和子载波的映射顺序。
% 这里我们构建涵盖 I/Q交换、Block块排列、子载波反序、共轭反相（Q取反）等的 16 种基准矩阵。
I_idx = 1:2:31;
Q_idx = 2:2:32;

maps = {};
maps{end+1} = struct('name', 'Interleaved [I,Q]',               'idx', 1:32, 'flip_Q', false);
maps{end+1} = struct('name', 'Interleaved [Q,I]',               'idx', reshape([Q_idx; I_idx], 1, 32), 'flip_Q', false);
maps{end+1} = struct('name', 'Interleaved RevSub [I,Q]',        'idx', reshape([fliplr(I_idx); fliplr(Q_idx)], 1, 32), 'flip_Q', false);
maps{end+1} = struct('name', 'Interleaved RevSub [Q,I]',        'idx', reshape([fliplr(Q_idx); fliplr(I_idx)], 1, 32), 'flip_Q', false);
maps{end+1} = struct('name', 'Block [I1..I16, Q1..Q16]',        'idx', [I_idx, Q_idx], 'flip_Q', false);
maps{end+1} = struct('name', 'Block [Q1..Q16, I1..I16]',        'idx', [Q_idx, I_idx], 'flip_Q', false);
maps{end+1} = struct('name', 'Block RevSub [I16..I1, Q16..Q1]', 'idx', [fliplr(I_idx), fliplr(Q_idx)], 'flip_Q', false);
maps{end+1} = struct('name', 'Block RevSub [Q16..Q1, I16..I1]', 'idx', [fliplr(Q_idx), fliplr(I_idx)], 'flip_Q', false);

% 添加共轭版本 (即 Q 反相)，这覆盖了硬判决取反或发射端共轭的情况
num_base = length(maps);
for m = 1:num_base
    new_map = maps{m};
    new_map.name = [new_map.name, ' (+ Conj/Flip Q)'];
    new_map.flip_Q = true;
    maps{end+1} = new_map;
end

num_maps = length(maps);
ref_matrices = cell(num_maps, 1);
for m = 1:num_maps
    idx = maps{m}.idx;
    M = ref_bits_matrix(idx, :); % 按索引重排矩阵行 (映射重排)
    if maps{m}.flip_Q
        % 将对应的 Q比特位翻转 (0变1，1变0)
        for r = 1:32
            if ismember(idx(r), Q_idx)
                M(r, :) = 1 - M(r, :);
            end
        end
    end
    ref_matrices{m} = M;
end

fprintf('已生成 %d 种扩展映射假说准备进行全面比对。\n', num_maps);

%% 2. 解析提取出的四个相位的假说文本
res_file = 'all_decoded_bits_4qam_4phases.txt';
if ~isfile(res_file)
    error('缺失 %s 文件', res_file);
end

fid = fopen(res_file, 'r');
file_matches = [];
current_file_idx = 0;
current_file_hyp = zeros(4, 32);

while ~feof(fid)
    line_str = strtrim(fgetl(fid));
    
    if startsWith(line_str, '--- File ')
        tokens = regexp(line_str, '\d+', 'match');
        if ~isempty(tokens)
            current_file_idx = str2double(tokens{1});
        end
    elseif startsWith(line_str, 'Phase ')
        % 提取: Phase x (xxx deg): 1100100...
        tokens = regexp(line_str, 'Phase (\d+).*:\s+([01]+)', 'tokens');
        if ~isempty(tokens)
            p_idx = str2double(tokens{1}{1});
            bit_str = tokens{1}{2};
            if length(bit_str) == 32
                current_file_hyp(p_idx, :) = bit_str - '0';
            end
            
            % 读取完第 4 个相位后，与所有的 ref 矩阵组合做匹配
            if p_idx == 4 && current_file_idx > 0
                best_err = inf;
                best_p = 0;
                best_t = 0;
                best_m = 0;
                
                for p = 1:4
                    hyp_bits = current_file_hyp(p, :).';
                    
                    for m = 1:num_maps
                        M = ref_matrices{m};
                        % 按列计算汉明距离
                        err_vector = sum(M ~= hyp_bits, 1);
                        
                        [min_e, min_t] = min(err_vector);
                        if min_e < best_err
                            best_err = min_e;
                            best_p = p;
                            best_t = min_t;
                            best_m = m;
                        end
                    end
                end
                
                % 记录: [文件ID, 相位ID, 时间T, 误差Err, 映射方式ID]
                file_matches = [file_matches; current_file_idx, best_p, best_t, best_err, best_m]; %#ok<AGROW>
                current_file_idx = 0;
            end
        end
    end
end
fclose(fid);

%% 3. 统计展示结果
fprintf('\n================ 扩展映射检验匹配验证结果 ================\n');
if isempty(file_matches)
    fprintf('警告：未能从文件中成功解析任何序列，请检查文件格式。\n');
    return;
end

% 匹配容错阈值（32 个 Bit 中，误码小于等于这个数的认为是真正的基准序列中的部分）
% 在做全面比对后，误差如果归零，阈值设小一点也能完全匹配上
match_threshold = 4;
valid_matches = file_matches(file_matches(:, 4) <= match_threshold, :);

fprintf('有效完美/近完美匹配数： %d / %d 个被测试文件（容许最大错误位 <= %d 位）\n', ...
    size(valid_matches, 1), size(file_matches, 1), match_threshold);

mean_err = mean(file_matches(:, 4));
fprintf('\n平均最匹配错误比特数: %.2f / 32\n', mean_err);

% 看看出现最多的映射方式
best_maps = file_matches(:, 5);
map_modes = mode(best_maps);
fprintf('\n>> 当下判定系统最有可能采用的内部映射策略是:【%s】\n', maps{map_modes}.name);

%% 4. 输出完整结果到新 txt 文件
out_res_file = 'pilot_match_extended_results.txt';
fid_out = fopen(out_res_file, 'w');
fprintf(fid_out, '================ 扩展映射验证完整结果 ================\n');
fprintf(fid_out, '测试总文件数：%d\n', size(file_matches, 1));
fprintf(fid_out, '推测系统最可能的映射策略：%s\n', maps{map_modes}.name);
fprintf(fid_out, '平均匹配误差：%.2f / 32 bits\n\n', mean_err);

fprintf(fid_out, '【所有被测试文件的最佳匹配详情】\n');
fprintf(fid_out, '文件编号  | 最佳相位假说     | 时间 T  | 误差/32 | 映射假说\n');
fprintf(fid_out, '------------------------------------------------------------------------------------\n');
for i = 1:size(file_matches, 1)
    f_idx = file_matches(i, 1);
    p_idx = file_matches(i, 2);
    t_idx = file_matches(i, 3);
    err_bits = file_matches(i, 4);
    m_idx = file_matches(i, 5);
    
    fprintf(fid_out, 'File %03d | Phase %d (%.0f°) | T=%03d | %02d/32  | %s\n', ...
        f_idx, p_idx, (p_idx-1)*90, t_idx, err_bits, maps{m_idx}.name);
end

edges = 0:2:32;
err_counts = histcounts(file_matches(:, 4), edges);
fprintf(fid_out, '\n【整体误差分布】\n');
for e = 1:length(err_counts)
    if err_counts(e) > 0
        fprintf(fid_out, ' 错 %2d~%2d 位: %d 个文件\n', edges(e), edges(e+1)-1, err_counts(e));
    end
end

fclose(fid_out);
fprintf('\n---> 扩展映射比对验证结果已经保存至新文件: %s\n', out_res_file);
