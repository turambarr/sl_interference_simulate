% 独立脚本：固定降采样倍率绘制星座图
% 用途：验证特定降采样因子 (D=3.2) 下的星座图效果
% 输入：IQ 文件路径，起始样点，降采样倍率

clear; clc; close all;

%% 1. 参数设置
inFile = 'sigtest1.iq'; 

% 采样范围设置
regionStart = 15530-874+10; 
L_read = 874;      % 读取长度
regionEnd   = regionStart + L_read; 

% 降采样设置
D = 6.828125;            % 降采样倍率 (Downsampling Factor)
                    % 意味着每 3.199 个原始点产出 1 个新点
                    % 新采样率 = 旧采样率 / 3.2

normalizeToUnit = true; 
removeMean = true;      

%% 2. 读取数据
[x, meta] = iq_read_int16_le(inFile, regionStart, local_num_samples_limit(inFile, L_read, regionStart));
L = meta.numSamplesRead;
if L == 0
    error('未读取到数据，请检查起始点是否越界');
end

x = double(x);

if normalizeToUnit
    x = x / 32768;
end
if removeMean
    x = x - mean(x);
end

fprintf('已读取数据: [%d, %d], 长度 %d\n', regionStart, regionStart+L-1, L);
fprintf('目标降采样倍率 D = %.2f\n', D);

% --- 准备两组数据进行对比 ---
% 1. 正常组 (假设文件是对齐的 I, Q, I, Q)
x_normal = x;

% 2. 错位组 (假设文件第一个样点多余/是Q，导致后面全错)
%    操作：丢弃第一个复数里的 I (其实是丢弃原始 int16 流的第1个点)
%    但在已经读成 complex 数组 x 后，这等价于重新构造
%    原始 int16流: [i0, q0, i1, q1, i2, q2...]
%    x 目前是:     (i0+jq0), (i1+jq1), (i2+jq2)...
%    如果真实是:    [垃圾, i0_真, q0_真, i1_真, q1_真...]
%    那我们需要：   real(x)里的q0, imag(x)里的i1 组合... 这很麻烦
%    最简单的办法：重新读，这次从 regionStart+1 (offset 1 sample? No, offset 2 bytes = 1 int16)
%    "把第一个样点q删掉" -> 意思是 "把第一个int16扔掉"

fprintf('正在重新读取以模拟错位 (Skip 1 int16)... \n');
% 偏移 2 字节 = 1 个 int16，从而让原来的 Q 变成 I
[x_shifted_raw, ~] = iq_read_int16_le(inFile, 0, 0, (regionStart*4) + 2); 
% 注意：上面的调用稍微有点 hacky。
% iq_read_int16_le(filename, startSample, numSamples, startBytes)
% 我们的 startSample 是基于复数点(4字节)的。
% 要想偏移半个复数点(2字节)，我们把 startSample 设为 0，把 startBytes 设为绝对偏移量

% 为了保证长度一致，我们只读 L_read 个
% 绝对偏移 = 100(Header) + regionStart*4 + 2(错位)
% 但 iq_read_int16_le 内部有 fseek(fid, startBytes + startSample*4, ...)
% 所以我们可以：
shift_offset_bytes = 100 + regionStart*4 + 2; % 100是文件头，2是我们要跳过的那个 "Q" (int16)
num_to_read = L_read; % 保持长度一致

fid = fopen(inFile, 'rb');
fseek(fid, shift_offset_bytes, 'bof');
raw_shifted = fread(fid, num_to_read * 2, 'int16');
fclose(fid);

i_sh = raw_shifted(1:2:end);
q_sh = raw_shifted(2:2:end);
len_sh = min(length(i_sh), length(q_sh));
x_shifted = complex(i_sh(1:len_sh), q_sh(1:len_sh));
x_shifted = double(x_shifted);

if normalizeToUnit, x_shifted = x_shifted / 32768; end
if removeMean,     x_shifted = x_shifted - mean(x_shifted); end

%% 3. 执行降采样
[P, Q] = rat(1 / D);
fprintf('使用有理数近似: P=%d, Q=%d (Ratio = %.4f)\n', P, Q, P/Q);

% 分别进行重采样
x_res_normal = resample(x_normal, P, Q);
x_res_shifted = resample(x_shifted, P, Q);


% % resample 函数引入的 FIR 滤波器延迟可能导致头尾出现靠近0的过渡点
% trim_len = 500;  % 截掉首尾各 50 个点
% if length(x_res_normal) > 2*trim_len
%     x_res_normal = x_res_normal(trim_len+1 : end-trim_len);
%     x_res_shifted = x_res_shifted(trim_len+1 : end-trim_len);
%     fprintf('已去除首尾各 %d 个样点的滤波过渡段以消除边缘效应。\n', trim_len);
% end

%% 4. 绘图对比
figure('Position', [100, 100, 1200, 600], 'Name', 'Check Sample Shift');

% === 子图1: Normal Reading (Start aligned) ===
subplot(1, 2, 1);
local_plot_const(x_res_normal, D, 'Normal Alignment');

% === 子图2: Shifted Reading (Skip 1 int16) ===
subplot(1, 2, 2);
local_plot_const(x_res_shifted, D, 'Shifted 1 int16 (Skip potential garbage Q)');



%% 内部绘图函数 (Nested function is not supported in simple script unless it's at the end, so we use inline logic or just append to end)
% 为了避免复杂性，这里直接把绘图逻辑放到辅助函数区域，或者在下面直接调用新定义的函数

function local_plot_const(sig, D, titleTag)
    L = length(sig);
    idx_odd = 1:2:L;
    idx_even = 2:2:L;
    
    plot(real(sig(idx_odd)), imag(sig(idx_odd)), 'r.', 'MarkerSize', 6); hold on;
    plot(real(sig(idx_even)), imag(sig(idx_even)), 'b.', 'MarkerSize', 6);
    
    axis square; grid on;
    limit = 0.5; 
    xlim([-limit, limit]);
    ylim([-limit, limit]);
    xlabel('I'); ylabel('Q');
    title(sprintf('%s\n(Red=Odd, Blue=Even)', titleTag), 'FontSize', 12, 'FontWeight', 'bold');
    legend('Odd', 'Even', 'Location', 'southoutside', 'Orientation', 'horizontal');
end

% (Existing helper functions should remain below)



%% 辅助函数
function n = local_num_samples_limit(inFile, reqSamples, startSample0)
    bytesPerComplexSample = 4;
    info = dir(inFile);
    if isempty(info)
        error('找不到文件: %s', inFile);
    end
    Ntotal = floor(info.bytes / bytesPerComplexSample);
    
    if startSample0 >= Ntotal
        n = 0;
    else
        available = Ntotal - startSample0;
        n = min(reqSamples, available);
    end
end

function [x, meta] = iq_read_int16_le(filename, startSample, numSamples, startBytes)
    if nargin < 4, startBytes = 0; end
    if nargin < 3, numSamples = []; end
    
    fid = fopen(filename, 'rb'); % 'rb' = Read Binary
    if fid < 0
        error('无法打开文件: %s', filename);
    end
    
    status = fseek(fid, startBytes + startSample * 4, 'bof');
    if status ~= 0
        fclose(fid);
        error('定位文件失败');
    end
    
   
    % 指定 'int16' => 'ieee-le'。
    
    if isempty(numSamples)
        raw = fread(fid, Inf, 'int16', 0, 'ieee-le'); % 显式指定 'ieee-le' (Little Endian)
    else
        raw = fread(fid, numSamples * 2, 'int16', 0, 'ieee-le'); 
    end
    fclose(fid);
    
    if isempty(raw)
        x = [];
        meta.numSamplesRead = 0;
        return;
    end
    
    i_data = raw(1:2:end);
    q_data = raw(2:2:end);
    len = min(length(i_data), length(q_data));
    x = complex(i_data(1:len), q_data(1:len));
    meta.numSamplesRead = len;
end
