fs_symbol = 60e6;
fs_source = 409.6e6;
fc_sync = 63e6;
hex_str = 'BD3CD0148871751F84CED8C1BE32AC96';
hex_chars = char(hex_str);
bits_dummy = zeros(1, 128);
for i = 1:length(hex_chars)
    bits_dummy((i-1)*4 + 1 : i*4) = dec2bin(hex2dec(hex_chars(i)), 4) - '0';
end
d_phi = (bits_dummy == 1) * (-pi/2) + (bits_dummy == 0) * (pi/2);
m_syms = exp(1j * cumsum(d_phi));
pss_base = [-m_syms, m_syms, m_syms, m_syms, m_syms, m_syms, m_syms, m_syms];
cp_syms = -pss_base(end - 47 : end);
pss_base_with_cp = [cp_syms, pss_base];
[P_res, Q_res] = rat(fs_source / fs_symbol);
pss_up = resample(pss_base_with_cp, P_res, Q_res);
t_local = (0:length(pss_up)-1).' / fs_source;
disp(size(pss_up));
disp(size(t_local));
try
    pss_local = pss_up .* exp(1j * 2 * pi * fc_sync * t_local);
    disp(size(pss_local));
catch me
    disp(me.message);
end
