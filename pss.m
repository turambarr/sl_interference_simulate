%生成PSS同步符号
Qpss_hex = 'C1B5D191024D3DC3F8EC52FAA16F3958';           %生成同步符号所用到的十六进制序列
Qpss_bin = [];         %存储二进制序列的数组
for i = 1 : length(Qpss_hex)
    if Qpss_hex(i) == '0'                       %将十六进制转化为二进制并存入数组中
        Qpss_bin = cat(2,Qpss_bin,[0 0 0 0]);
    elseif Qpss_hex(i) == '1'
        Qpss_bin = cat(2,Qpss_bin,[0 0 0 1]);
    elseif Qpss_hex(i) == '2'
        Qpss_bin = cat(2,Qpss_bin,[0 0 1 0]);
    elseif Qpss_hex(i) == '3'
        Qpss_bin = cat(2,Qpss_bin,[0 0 1 1]);
    elseif Qpss_hex(i) == '4'
        Qpss_bin = cat(2,Qpss_bin,[0 1 0 0]);
    elseif Qpss_hex(i) == '5'
        Qpss_bin = cat(2,Qpss_bin,[0 1 0 1]);
    elseif Qpss_hex(i) == '6'
        Qpss_bin = cat(2,Qpss_bin,[0 1 1 0]);
    elseif Qpss_hex(i) == '7'
        Qpss_bin = cat(2,Qpss_bin,[0 1 1 1]);
    elseif Qpss_hex(i) == '8'
        Qpss_bin = cat(2,Qpss_bin,[1 0 0 0]);
    elseif Qpss_hex(i) == '9'
        Qpss_bin = cat(2,Qpss_bin,[1 0 0 1]);
    elseif Qpss_hex(i) == 'A'
        Qpss_bin = cat(2,Qpss_bin,[1 0 1 0]);
    elseif Qpss_hex(i) == 'B'
        Qpss_bin = cat(2,Qpss_bin,[1 0 1 1]);
    elseif Qpss_hex(i) == 'C'
        Qpss_bin = cat(2,Qpss_bin,[1 1 0 0]);
    elseif Qpss_hex(i) == 'D'
        Qpss_bin = cat(2,Qpss_bin,[1 1 0 1]);
    elseif Qpss_hex(i) == 'E'
        Qpss_bin = cat(2,Qpss_bin,[1 1 1 0]);
    elseif Qpss_hex(i) == 'F'
        Qpss_bin = cat(2,Qpss_bin,[1 1 1 1]);
    end
end
b = flip(2 * Qpss_bin - 1);                %得出公式中的b序列,flip函数用于将序列逆向排序
pss = [];                            %存储生成的N/8长度的训练序列
reverse_pss = [];
b_sum_list = [];
for k = 1 : N / 8
    b_sum = 0;                       %存储b序列之和的变量
    reverse_b_sum = 0;
    for n = 1 : k
        b_sum = b_sum + b(n);        %将对应个数的b的元素相加
        reverse_b_sum = reverse_b_sum - b(n);
    end
    pss = cat(2,pss,exp(1j * pi * (- 1 / 4 - 1 / 2 * b_sum)));          %将该符号填入数组中
    reverse_pss = cat(2, reverse_pss, exp(1j * pi * (- 1 / 4 - 1 / 2 * reverse_b_sum)));
    b_sum_list = cat(2, b_sum_list, b_sum);
end
strain_list = [];         %构造pss训练序列
t = 0 : 1 / Fs : (N - 1) / Fs;      %采样时刻
for i = 1 : N
    mysum = 0;        %存储某一时刻信号值的变量
    for k = 0 : N - 1       %计算某一时刻的信号值
        if k < N / 8
            mysum = mysum + sinc(t(i) * Fs - k) * (- pss(k + 1));      %第一个数据块符号与后面7个相反
        else
            mysum = mysum + sinc(t(i) * Fs - k) * pss(mod(k,N / 8) + 1);      %后面7个数据块相同
        end
    end
    strain_list = cat(2,strain_list,mysum);     %将该时刻值存入数组中
end


%生成SSS同步符号
Qsss_hex = 'BD565D5064E9B3A94958F28624DED560946199F5B40F0E4FB5EFCB473B4C24B2D1E0BD01A6A04D5017DE91A8ECC0DA09EBFE57F9F1B44C532F161C583A42490A5C09F2A117F9A28F9B2FD547A74C44BABB4BE85DA6A62B1235E2AD084C00180142A8F7F357DEC4F31316BC58FA404909A3FCA7F88E421902B6A2580AE8030803F65809DB347F590DBC46F010EBE3A25C060D74429FC46BDF9B63719279798D232C5ABA274122FF66AD7E449F44CB40C49C24A1E2629F5BFE82CE531FDC34F8C64A43A963F40D5B71BDE6FB2F13492D6F2E8544B21D449722C635180342CD0026A1E7F7E80E91B175E852F919767E5AF9B6E909AF362F5218E2B908DC005803';           %生成同步符号所用到的十六进制序列
Qsss_Qua = [];         %存储四进制序列的数组
for i = 1 : length(Qsss_hex)
    if Qsss_hex(i) == '0'                       %将十六进制转化为四进制并存入数组中
        Qsss_Qua = cat(2,Qsss_Qua,[0 0]);
    elseif Qsss_hex(i) == '1'
        Qsss_Qua = cat(2,Qsss_Qua,[0 1]);
    elseif Qsss_hex(i) == '2'
        Qsss_Qua = cat(2,Qsss_Qua,[0 2]);
    elseif Qsss_hex(i) == '3'
        Qsss_Qua = cat(2,Qsss_Qua,[0 3]);
    elseif Qsss_hex(i) == '4'
        Qsss_Qua = cat(2,Qsss_Qua,[1 0]);
    elseif Qsss_hex(i) == '5'
        Qsss_Qua = cat(2,Qsss_Qua,[1 1]);
    elseif Qsss_hex(i) == '6'
        Qsss_Qua = cat(2,Qsss_Qua,[1 2]);
    elseif Qsss_hex(i) == '7'
        Qsss_Qua = cat(2,Qsss_Qua,[1 3]);
    elseif Qsss_hex(i) == '8'
        Qsss_Qua = cat(2,Qsss_Qua,[2 0]);
    elseif Qsss_hex(i) == '9'
        Qsss_Qua = cat(2,Qsss_Qua,[2 1]);
    elseif Qsss_hex(i) == 'A'
        Qsss_Qua = cat(2,Qsss_Qua,[2 2]);
    elseif Qsss_hex(i) == 'B'
        Qsss_Qua = cat(2,Qsss_Qua,[2 3]);
    elseif Qsss_hex(i) == 'C'
        Qsss_Qua = cat(2,Qsss_Qua,[3 0]);
    elseif Qsss_hex(i) == 'D'
        Qsss_Qua = cat(2,Qsss_Qua,[3 1]);
    elseif Qsss_hex(i) == 'E'
        Qsss_Qua = cat(2,Qsss_Qua,[3 2]);
    elseif Qsss_hex(i) == 'F'
        Qsss_Qua = cat(2,Qsss_Qua,[3 3]);
    end
end
sss = flip(exp(1j * Qsss_Qua * pi / 2));                %得出sss序列,flip函数用于将序列逆向排序
final_sss = [];                                         %存储重新安置训练序列在子载波上位置的序列
final_sss = cat(2,final_sss,sss((N / 2 - 1) : (N - 4)),zeros(1,4),sss(1 : (N / 2 - 2)));  %将训练序列重新排序，中间留出预防中心频率泄露的四个空子载波
SSSsymbol = ifft(final_sss,N);                  %生成SSS符号