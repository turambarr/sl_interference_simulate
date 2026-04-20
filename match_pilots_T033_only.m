% match_pilots_T033_only.m
% 定向查错脚本：强制要求采用 Interleaved [I,Q] 策略，且紧盯时间索引 T=33
% 旨在输出直观的“按位比对图”，帮助肉眼排查到底哪几个位反反复复出错。

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
        parts = split(str, '=');
        if length(parts) == 2
            % 转为 char 数组并只保留有效的十六进制字符
            hex_str = regexprep(char(parts(2)), '[^0-9a-fA-F]', '');
            pilot_hex{pilot_idx} = hex_str;
            pilot_idx = pilot_idx + 1;
        end
    end
end

if pilot_idx <= 16
    error('导频.txt 中未找到完整的 16 个子载波定义。');
end

hex_dict = containers.Map(...
    {'0','1','2','3','4','5','6','7','8','9','A','B','C','D','E','F',...
     'a','b','c','d','e','f'}, ...
    {[0 0 0 0],[0 0 0 1],[0 0 1 0],[0 0 1 1],...
     [0 1 0 0],[0 1 0 1],[0 1 1 0],[0 1 1 1],...
     [1 0 0 0],[1 0 0 1],[1 0 1 0],[1 0 1 1],...
     [1 1 0 0],[1 1 0 1],[1 1 1 0],[1 1 1 1],...
     [1 0 1 0],[1 0 1 1],[1 1 0 0],[1 1 0 1],[1 1 1 0],[1 1 1 1]});

num_hex_chars = length(pilot_hex{1});
num_bits = num_hex_chars * 4;
num_symbols = floor(num_bits / 2);

ref_bits_matrix = zeros(32, num_symbols);

for k = 1:16
    hex_str = pilot_hex{k};
    bits_k = zeros(1, length(hex_str)*4);
    for c = 1:length(hex_str)
        bits_k((c-1)*4+1 : c*4) = hex_dict(hex_str(c));
    end
    
    cur_symbols = floor(length(bits_k) / 2);
    use_sym = min(cur_symbols, num_symbols);
    
    % 改为单纯的 Interleaved [Q,I] 排列：
    % 即 16 个子载波，每个子载波先排 Q 比特，再排 I 比特
    ref_bits_matrix(2*k - 1, 1:use_sym) = bits_k(2:2:2*use_sym); % 奇数行为 Q
    ref_bits_matrix(2*k, 1:use_sym)     = bits_k(1:2:2*use_sym); % 偶数行为 I
end
num_symbols = use_sym;
ref_bits_matrix = ref_bits_matrix(:, 1:num_symbols);

%% 2. 定向固定比对 (T_target = 33, 映射 = 单纯的 Interleaved [Q,I])
target_T = 33;
if target_T > num_symbols
    error('参考矩阵长度不足 T=%d', target_T);
end
ref_target_bits = ref_bits_matrix(:, target_T);

res_file = 'all_decoded_bits_4qam_4phases.txt';
if ~isfile(res_file)
    error('缺失 %s 文件', res_file);
end

fid = fopen(res_file, 'r');
current_file_idx = 0;
current_file_hyp = zeros(4, 32);

out_res_file = 'pilot_match_T033_only.txt';
fid_out = fopen(out_res_file, 'w');

fprintf(fid_out, '================ 定向侦查：单纯的 Interleaved [Q,I], 强制固定 T=033 ================\n');
fprintf(fid_out, '说明：这里将所有文件强制与目标时刻 T=33 的理想32位导频进行逐位对齐。\n');
fprintf(fid_out, '如果解出的序列匹配该位，显示为对应的 0 或 1；如果不匹配，将在该位用字母 X 高亮标出。\n\n');

% 输出标尺，方便数刻度
fprintf(fid_out, '位索引标尺:     12345678901234567890123456789012\n');
fprintf(fid_out, '理想目标(T=33): %s\n\n', sprintf('%d', ref_target_bits));
fprintf(fid_out, '文件编号  | 最优相位 | 误差/32 |  按位重合情况透视图 (X为报错位)\n');
fprintf(fid_out, '-----------------------------------------------------------------------------------\n');

total_err_count = 0;
file_cnt = 0;
err_arr = [];
all_best_bits = zeros(173, 32);

while ~feof(fid)
    line_str = strtrim(fgetl(fid));
    
    if startsWith(line_str, '--- File ')
        tokens = regexp(line_str, '\d+', 'match');
        if ~isempty(tokens)
            current_file_idx = str2double(tokens{1});
        end
    elseif startsWith(line_str, 'Phase ')
        tokens = regexp(line_str, 'Phase (\d+).*:\s+([01]+)', 'tokens');
        if ~isempty(tokens)
            p_idx = str2double(tokens{1}{1});
            bit_str = tokens{1}{2};
            if length(bit_str) == 32
                current_file_hyp(p_idx, :) = bit_str - '0';
            end
            
            % 读取完第 4 个相位，强制进行 T=33 对比
            if p_idx == 4 && current_file_idx > 0
                best_err = inf;
                best_p = 0;
                best_hyp_bits = zeros(32, 1);
                
                % 在4个假设中挑出一个跟它最像的（容错最好的）
                for p = 1:4
                    hyp_bits = current_file_hyp(p, :).';
                    err = sum(ref_target_bits ~= hyp_bits);
                    if err < best_err
                        best_err = err;
                        best_p = p;
                        best_hyp_bits = hyp_bits;
                    end
                end
                
                % 生成透视字符串
                vis_str = char(zeros(1, 32));
                for b = 1:32
                    if ref_target_bits(b) == best_hyp_bits(b)
                        vis_str(b) = num2str(best_hyp_bits(b));
                    else
                        vis_str(b) = 'X';
                    end
                end
                
                all_best_bits(current_file_idx, :) = best_hyp_bits.';
                
                fprintf(fid_out, 'File %03d | Ph %d (%03d°) | %02d/32  |  %s\n', ...
                    current_file_idx, best_p, (best_p-1)*90, best_err, vis_str);
                
                total_err_count = total_err_count + best_err;
                file_cnt = file_cnt + 1;
                err_arr = [err_arr; best_err]; %#ok<AGROW>
                current_file_idx = 0;
            end
        end
    end
end
fclose(fid);

% 统计常犯错的比特坑位 (针对该 T=33 假设下)
error_per_bit = sum(ref_target_bits.' ~= all_best_bits(1:file_cnt, :), 1);
[sorted_err, sorted_idx] = sort(error_per_bit, 'descend');

fprintf(fid_out, '\n================ 强制约束统计信息 ================\n');
fprintf(fid_out, '1. 总体平均误差: %.2f / 32 bits\n', total_err_count / file_cnt);

fprintf(fid_out, '\n2. 【常报错比特位置侦测排行】 (在所有173个文件中，累计错误率最高的位):\n');
for i = 1:length(sorted_idx)
    if sorted_err(i) > 0
        bit_idx = sorted_idx(i);
        sub_c_idx = ceil(bit_idx / 2);
        
        % 因为现在是 [Q, I] 映射，所以奇数索引(1,3,5...)是 Q 通道，偶数索引是 I 通道
        iq_label = 'Q'; 
        if mod(bit_idx, 2) == 0
            iq_label = 'I'; 
        end
        
        fprintf(fid_out, '   比特索引 %02d (即第 %d 根载波的 %s 通道) -> 共报错 %03d 次 (错误率 %05.2f%%)\n', ...
            bit_idx, sub_c_idx, iq_label, sorted_err(i), sorted_err(i)/file_cnt*100);
    end
end

fclose(fid_out);
fprintf('专门针对单纯的 [ Interleaved Q,I 且 T=033 ] 的强制逐位对比透视图已生成。\n请查阅文件：%s\n', out_res_file);