% test_mseq_qpsk_correlation.m
% 仿真脚本：生成两段不同的m序列进行QPSK调制
% 构造组合信号：[S1, S1, S1, S2, S1, S1]
% 并验证互相关检测效果（分别用 S1 和 S2 作为模板）

clear; clc; close all;

%% 1. 参数设置
m_order = 8;             % m序列阶数 (长度 = 2^8 - 1 = 255 bits)
num_samples = 255;       % 截断或补齐用于QPSK调制的比特数
snr_db = 20;             % 混合信号的信噪比 (dB)，模拟真实环境

%% 2. 生成两段不同的 m序列 (Bits)
% 使用两组不同的多项式反馈抽头
% 阶数8: 255长度
% Poly1: x^8 + x^4 + x^3 + x^2 + 1 (Taps: [8 4 3 2])
% Poly2: x^8 + x^6 + x^5 + x + 1   (Taps: [8 6 5 1])

bits1 = generate_m_sequence(m_order, [8 4 3 2]);
bits2 = generate_m_sequence(m_order, [8 6 5 1]);

% 确保长度是偶数以便QPSK调制 (丢弃最后一个或补0)
bits1 = bits1(1:254);
bits2 = bits2(1:254);

%% 3. QPSK 调制
% 映射规则: 2 bits -> 1 symbol
s1 = qpsk_modulate(bits1);
s2 = qpsk_modulate(bits2);

% 归一化幅度
s1 = s1 / mean(abs(s1));
s2 = s2 / mean(abs(s2));

fprintf('信号生成完成：\n');
fprintf('  S1 长度: %d 样点 (基于 m序列1)\n', length(s1));
fprintf('  S2 长度: %d 样点 (基于 m序列2)\n', length(s2));

%% 4. 构造混合信号 S3
% 模式: S1, S1, S1, S2, S1, S1
s3_clean = [s1; s1; s1; s2; s1; s1];
L_segment = length(s1);

% 添加高斯白噪声
s3 = awgn(s3_clean, snr_db, 'measured');

fprintf('  S3 构造完成 (模式: S1x3 + S2 + S1x2), 总长度: %d\n', length(s3));

%% 5. 执行互相关 (仿照 cross_corr_sig1_sig2.m)

% --- Case A: 用 S1 去搜 S3 ---
[xc1, lags1, max1, best_lag1] = perform_normalized_correlation(s3, s1);

% --- Case B: 用 S2 去搜 S3 ---
[xc2, lags2, max2, best_lag2] = perform_normalized_correlation(s3, s2);

%% 6. 绘图展示
figure('Position', [100, 100, 1000, 800], 'Name', 'M-Seq QPSK Correlation Test');

% 子图1: S3 时域波形
subplot(3, 1, 1);
plot(abs(s3), 'Color', [0.6 0.6 0.6]); 
hold on; 
% 简单的标注一下真实位置区域
x_ticks = 1:L_segment:length(s3);
for k = 1:length(x_ticks)
    xline(x_ticks(k), 'k--');
end
title(['信号 S3 幅度 (结构: S1 - S1 - S1 - S2 - S1 - S1)  [SNR=' num2str(snr_db) 'dB]']);
xlim([1 length(s3)]); grid on;
% 标注块名称
text(L_segment*0.5, 1.5, 'S1', 'Color', 'b', 'Ho', 'center');
text(L_segment*1.5, 1.5, 'S1', 'Color', 'b', 'Ho', 'center');
text(L_segment*2.5, 1.5, 'S1', 'Color', 'b', 'Ho', 'center');
text(L_segment*3.5, 1.5, 'S2', 'Color', 'r', 'Ho', 'center', 'FontWeight', 'bold');
text(L_segment*4.5, 1.5, 'S1', 'Color', 'b', 'Ho', 'center');
text(L_segment*5.5, 1.5, 'S1', 'Color', 'b', 'Ho', 'center');

% 子图2: S1 模板匹配结果
subplot(3, 1, 2);
plot(lags1, xc1, 'b', 'LineWidth', 1);
title('用 S1 作为模板进行互相关 (预期: 5个强峰, 1个弱区)');
ylabel('模值 (Abs)');
grid on; axis tight;
ylim([0 1.2]);
xline(best_lag1, 'r--'); % 标记最大峰

% 子图3: S2 模板匹配结果
subplot(3, 1, 3);
plot(lags2, xc2, 'r', 'LineWidth', 1);
title('用 S2 作为模板进行互相关 (预期: 仅 1 个强峰)');
xlabel('Lag (延迟采样数)');
ylabel('模值 (Abs)');
grid on; axis tight;
ylim([0 1.2]);
xline(best_lag2, 'b--'); % 标记最大峰

fprintf('\n=== 测试结果 ===\n');
fprintf('Case 1 (S1 vs S3) 最大匹配度: %.4f @ Lag %d\n', max1, best_lag1);
fprintf('Case 2 (S2 vs S3) 最大匹配度: %.4f @ Lag %d\n', max2, best_lag2);


%% === 本地函数 ===

function [xc_abs, lags, max_val, best_lag] = perform_normalized_correlation(long_sig, template_sig)
    % 核心算法与 cross_corr_sig1_sig2.m 保持一致
    % 1. 计算原始 xcorr
    [xc, lags] = xcorr(long_sig, template_sig);
    
    % 2. 归一化 (除以模板能量)
    template_energy = sum(abs(template_sig).^2);
    if template_energy > 0
        xc = xc / template_energy;
    end
    
    % 3. 取模
    xc_abs = abs(xc);
    
    [max_val, idx] = max(xc_abs);
    best_lag = lags(idx);
end

function syms = qpsk_modulate(bits)
    % 简单的 QPSK 映射
    % 00 -> 1+1j
    % 01 -> -1+1j
    % 11 -> -1-1j
    % 10 -> 1-1j (Gray code wise or straightforward)
    % 这里我们用简单的直接映射:
    % Odd bit -> I, Even bit -> Q (0->1, 1->-1)
    
    i_bits = bits(1:2:end);
    q_bits = bits(2:2:end);
    
    % BPSK map: 0->1, 1->-1
    i_val = 1 - 2*i_bits;
    q_val = 1 - 2*q_bits;
    
    syms = i_val + 1j * q_val;
end

function seq = generate_m_sequence(order, taps)
    % 简单的 LFSR 实现
    % state 初始化为非零
    state = ones(1, order);
    len = 2^order - 1;
    seq = zeros(len, 1);
    
    for k = 1:len
        seq(k) = state(end);
        % 反馈计算 (XOR sum of taps)
        % taps是指数位置，对应 state 索引需要转换
        % 比如 x^8 + x^4... 对应的寄存器位
        feedback = 0;
        for t = taps
            % taps通常指由高到低，这里简化：
            % 假设 state 是 [x_1, x_2 ... x_n]
            % taps 只是参与 XOR 的位置
            if t <= order
                feedback = xor(feedback, state(t));
            end
        end
        % 移位
        state = [feedback, state(1:end-1)];
    end
end
