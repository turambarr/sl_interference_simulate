clear; clc;
inFile1 = '20250912222305_part1.iq';
inFile2 = 'ofdmcp_interference_3dB.iq';

fid = fopen(inFile1, 'rb');
r1 = fread(fid, 1000000, 'int16');
fclose(fid);
x1 = r1(1:2:end) + 1j * r1(2:2:end);

fid = fopen(inFile2, 'rb');
r2 = fread(fid, 1000000, 'int16');
fclose(fid);
x2 = r2(1:2:end) + 1j * r2(2:2:end);

nfft = 8192;
X1 = pwelch(x1, 1024, 512, nfft, 409.6e6);
X2 = pwelch(x2, 1024, 512, nfft, 409.6e6);

diff_db = 10*log10(sum(abs(X1))) - 10*log10(sum(abs(X2)));
corr_val = max(abs(xcorr(x1(1:1000), x2(1:1000)))) / sqrt(sum(abs(x1(1:1000)).^2)*sum(abs(x2(1:1000)).^2));

fprintf('Diff overall power: %f dB\n', diff_db);
fprintf('Correlation coeff of early signal: %f\n', corr_val);

% Find where they differ
diff_idx = find(x1 ~= x2, 1);
fprintf('First diff index: %d\n', diff_idx);
