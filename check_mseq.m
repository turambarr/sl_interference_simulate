function result = check_mseq(seq, assumeFullPeriod)
% CHECK_MSEQ  检验某段序列是否为 m序列（或像LFSR/m序列片段）
%
% 用法:
%   r = check_mseq([1 0 1 1 0 ...]);          % 默认按“完整周期”尝试严格判断
%   r = check_mseq(hexStr);                    % hex字符串输入
%   r = check_mseq(seq, false);                % 不假设完整周期，只做结构性判断
%
% 输入:
%   seq             : 二进制向量(0/1) 或 十六进制字符串(如 '592D7A79...')
%   assumeFullPeriod: true/false (默认 true)
%
% 输出:
%   result: 结构体，含线性复杂度、BM反馈多项式、自相关、严格/非严格判断等
%
% 说明:
%   - 严格 m序列判断需要“完整周期”序列，长度必须满足 N = 2^m - 1
%   - 截断片段只能做“像不像LFSR/m序列”的结构性判断（BM线性复杂度等）
%


    if nargin < 2
        assumeFullPeriod = true;
    end

    % 1) 解析输入为 0/1 行向量
    s = parse_input_to_bits(seq);
    N = numel(s);

    % 2) Berlekamp-Massey (GF(2))，求最短LFSR阶数和连接多项式
    [L, C] = bm_gf2(s);     % C(1)=1, 长度为 L+1
    % 对应递推关系（GF(2)）:
    % s(n) = C(2)*s(n-1) + C(3)*s(n-2) + ... + C(L+1)*s(n-L)  (mod 2)

    % 3) 周期循环自相关（±1映射）
    %    x = +1/-1，m序列完整周期时应满足 R(0)=N, 其余 R(t)=-1
    R = cyclic_autocorr_pm1(s);

    % 4) 平衡性（完整周期m序列中两符号个数只差1）
    n1 = sum(s == 1);
    n0 = sum(s == 0);
    balance_abs_diff = abs(n1 - n0);

    % 5) 长度是否可能是完整m序列周期
    [isLenMseq, mCandidate] = is_mseq_period_length(N);

    % 6) 严格判定（只有假设完整周期时有意义）
    strict_len_ok   = isLenMseq;
    strict_bal_ok   = (balance_abs_diff == 1);
    strict_corr_ok  = (R(1) == N) && all(R(2:end) == -1);
    strict_bm_ok    = isLenMseq && (L == mCandidate);

    is_mseq_strict = false;
    if assumeFullPeriod
        is_mseq_strict = strict_len_ok && strict_bal_ok && strict_corr_ok && strict_bm_ok;
    end

    % 7) 结构性判断（截断片段适用）
    %    对随机序列，长度128时线性复杂度通常接近64左右；
    %    若L显著偏小，说明更像LFSR生成。
    is_lfsr_like = (L <= max(8, floor(N/4)));   % 经验阈值，可自行调
    lfsr_like_score = 1 - min(L / max(1, N/2), 1);  % 粗略评分, 越大越像简单LFSR

    % 8) 如果有通信工具箱，尝试检查多项式是否本原（可选）
    primitiveCheckAvailable = (exist('gfprimck', 'file') == 2);
    primitive_ok = NaN;
    primitive_msg = 'gfprimck 不可用（需要 Communications Toolbox），跳过本原多项式检查';
    if primitiveCheckAvailable && L >= 2
        % gfprimck 一般需要降幂排列，如 [1 ... 1] 表示 x^L + ... + 1
        % 这里 C 本身就是 [1 c1 c2 ... cL]，对应 x^L + c1*x^(L-1)+...+cL
        % 但BM得到的C与具体LFSR实现方向有关，这里仅作尝试性检查。
        try
            primitive_ok = gfprimck(C, 2); %#ok<GFLD>
            if primitive_ok
                primitive_msg = 'gfprimck: 该连接多项式判定为本原（方向匹配前提下）';
            else
                primitive_msg = 'gfprimck: 该连接多项式不是本原（或方向定义不匹配）';
            end
        catch ME
            primitive_ok = NaN;
            primitive_msg = ['gfprimck 调用失败: ' ME.message];
        end
    end

    % 9) 打包结果
    result = struct();
    result.N = N;
    result.bits = s;
    result.ones_count = n1;
    result.zeros_count = n0;
    result.balance_abs_diff = balance_abs_diff;

    result.linear_complexity_L = L;
    result.connection_poly_C = C;  % [1 c1 ... cL]
    result.connection_poly_desc = poly_to_string(C);

    result.cyclic_autocorr_pm1 = R;
    result.cyclic_autocorr_two_level = (R(1) == N) && all(R(2:end) == -1);

    result.length_is_2m_minus_1 = isLenMseq;
    result.m_candidate_if_full_period = mCandidate;

    result.assumeFullPeriod = assumeFullPeriod;
    result.strict_checks = struct( ...
        'length_ok', strict_len_ok, ...
        'balance_ok', strict_bal_ok, ...
        'corr_ok', strict_corr_ok, ...
        'bm_ok', strict_bm_ok);

    result.is_mseq_strict = is_mseq_strict;

    result.is_lfsr_like = is_lfsr_like;
    result.lfsr_like_score = lfsr_like_score;

    result.primitive_check_available = primitiveCheckAvailable;
    result.primitive_ok = primitive_ok;
    result.primitive_msg = primitive_msg;

    % 10) 命令行打印摘要
    fprintf('\n===== m序列检验结果 =====\n');
    fprintf('长度 N = %d\n', N);
    if isLenMseq
        fprintf('N 满足 2^m-1, 候选 m = %d\n', mCandidate);
    else
        fprintf('N 不满足 2^m-1（因此不可能是“完整周期”m序列）\n');
    end
    fprintf('0/1 个数: zeros=%d, ones=%d, |diff|=%d\n', n0, n1, balance_abs_diff);
    fprintf('BM线性复杂度 L = %d\n', L);
    fprintf('BM连接多项式: %s\n', result.connection_poly_desc);

    if assumeFullPeriod
        fprintf('严格判定 is_mseq_strict = %d\n', is_mseq_strict);
        fprintf('  - length_ok = %d\n', strict_len_ok);
        fprintf('  - balance_ok = %d\n', strict_bal_ok);
        fprintf('  - corr_ok    = %d\n', strict_corr_ok);
        fprintf('  - bm_ok      = %d\n', strict_bm_ok);
    else
        fprintf('未假设完整周期，跳过严格m序列判定。\n');
    end

    fprintf('结构性判断 is_lfsr_like = %d (score=%.3f)\n', is_lfsr_like, lfsr_like_score);
    fprintf('%s\n', primitive_msg);
    fprintf('========================\n\n');
end


%% ========= 子函数：输入解析 =========
function s = parse_input_to_bits(seq)
    if ischar(seq) || isstring(seq)
        hexStr = upper(strtrim(char(seq)));
        hexStr = regexprep(hexStr, '^0X', '');
        hexStr = regexprep(hexStr, '\s+', '');
        if isempty(hexStr) || ~all(ismember(hexStr, ['0':'9' 'A':'F']))
            error('字符串输入必须是十六进制字符串（可含空格）。');
        end
        % 每个hex字符 -> 4bit（MSB->LSB）
        s = zeros(1, 4*numel(hexStr));
        idx = 1;
        for k = 1:numel(hexStr)
            v = hex2dec(hexStr(k));
            b = bitget(uint8(v), 4:-1:1); % MSB到LSB
            s(idx:idx+3) = double(b);
            idx = idx + 4;
        end
    elseif isnumeric(seq) || islogical(seq)
        s = seq(:).';
        if ~all(ismember(s, [0 1]))
            error('数值输入必须是0/1序列。');
        end
        s = double(s);
    else
        error('不支持的输入类型。请输入0/1向量或十六进制字符串。');
    end
end


%% ========= 子函数：BM算法（GF(2)） =========
function [L, C] = bm_gf2(s)
% Berlekamp-Massey over GF(2)
% 输入 s: 0/1行向量
% 输出:
%   L: 线性复杂度
%   C: 连接多项式系数 [1 c1 c2 ... cL]
%
% 递推形式（GF(2)）:
%   s(n) = c1*s(n-1) + c2*s(n-2) + ... + cL*s(n-L)

    s = s(:).';
    N = numel(s);

    Cfull = zeros(1, N); Cfull(1) = 1;  % 当前多项式
    Bfull = zeros(1, N); Bfull(1) = 1;  % 上次更新时的多项式

    L = 0;
    m = -1;  % 上次发生L更新的位置（0-based）

    for n = 0:(N-1)  % 0-based
        % discrepancy d = s(n) + sum_{i=1..L} C(i+1)*s(n-i) mod 2
        d = s(n+1);
        for i = 1:L
            d = xor(d, Cfull(i+1) & s(n-i+1));
        end

        if d == 1
            T = Cfull;
            shift = n - m;  % 要把B左移shift位并异或到C上

            % C = C + x^(shift) * B (mod 2)
            for j = 1:(N-shift)
                Cfull(j+shift) = xor(Cfull(j+shift), Bfull(j));
            end

            if 2*L <= n
                L = n + 1 - L;
                Bfull = T;
                m = n;
            end
        end
    end

    C = Cfull(1:L+1);
    C = double(C);
end


%% ========= 子函数：循环自相关（±1映射） =========
function R = cyclic_autocorr_pm1(s)
% s为0/1，映射到x in {+1,-1}
% 这里用 0->+1, 1->-1（取反也不影响两值自相关特性）
    s = s(:).';
    N = numel(s);
    x = 1 - 2*s;  % 0->+1, 1->-1

    R = zeros(1, N);
    for tau = 0:(N-1)
        idx = mod((0:N-1) + tau, N) + 1;
        R(tau+1) = sum(x .* x(idx));
    end
end


%% ========= 子函数：长度是否为2^m-1 =========
function [tf, m] = is_mseq_period_length(N)
    val = N + 1;
    m = NaN;
    if val <= 0
        tf = false;
        return;
    end
    % 检查 val 是否为2的幂
    tf = (bitand(uint64(val), uint64(val - 1)) == 0);
    if tf
        m = round(log2(double(val)));
    end
end


%% ========= 子函数：多项式字符串显示 =========
function str = poly_to_string(C)
% C = [1 c1 c2 ... cL], 表示 x^L + c1*x^(L-1) + ... + cL
    L = numel(C) - 1;
    if L < 0
        str = 'N/A';
        return;
    elseif L == 0
        str = '1';
        return;
    end

    terms = {};
    terms{end+1} = sprintf('x^%d', L);
    for i = 1:L
        if C(i+1) == 1
            p = L - i;
            if p > 1
                terms{end+1} = sprintf('x^%d', p);
            elseif p == 1
                terms{end+1} = 'x';
            else
                terms{end+1} = '1';
            end
        end
    end
    str = strjoin(terms, ' + ');
end