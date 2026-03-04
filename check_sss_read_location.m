% check_sss_read_location.m
% Visualizes the entire signal and highlights the SSS read region

clear; clc; close all;

inFile = 'sigtest8.iq';

% Current parameters in sss_demodulation.m
pss_start_sample = 15564;
% startSample in script (PSS Start + PSS Len + Offset)
current_read_start = 15564 + 874*8 + 328; 
current_read_len   = 7500;

% Read the WHOLE file
d = dir(inFile);
total_samples = d.bytes / 4; % int16 complex (2+2 bytes)
[x_whole, ~] = iq_read_int16_le(inFile, 1, total_samples);
x_whole = double(x_whole);
x_whole = x_whole - mean(x_whole);

figure('Position', [100, 100, 1200, 600], 'Name', 'Check SSS Read Location');
plot(abs(x_whole), 'b'); hold on;
grid on;
title('Whole File Amplitude with Read Window Highlight');
xlabel('Sample Index'); ylabel('Amplitude');

% Highlight PSS region (Approx)
pss_len = 874*8;
x_pss = [pss_start_sample, pss_start_sample, pss_start_sample+pss_len, pss_start_sample+pss_len];
y_lims = ylim;
y_patch = [y_lims(1), y_lims(2), y_lims(2), y_lims(1)];
patch(x_pss, y_patch, 'g', 'FaceAlpha', 0.2, 'EdgeColor', 'none');

% Highlight SSS Read region
read_end = current_read_start + current_read_len;
if read_end > total_samples
    read_end = total_samples; 
end
x_read = [current_read_start, current_read_start, read_end, read_end];
patch(x_read, y_patch, 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'r', 'LineWidth', 2);

% Add text labels
text(pss_start_sample, y_lims(2)*0.9, 'PSS Start', 'Color', 'g', 'FontWeight', 'bold');
text(current_read_start, y_lims(2)*0.8, 'Read Start', 'Color', 'r', 'FontWeight', 'bold');
legend('Signal Amplitude', 'PSS Region (Green)', 'Current Read Window (Red)');

fprintf('File Total Samples: %d\n', total_samples);
fprintf('PSS Start: %d\n', pss_start_sample);
fprintf('Current Read Start: %d\n', current_read_start);
fprintf('Current Read End: %d\n', read_end);

if current_read_start > total_samples
    warning('Start sample is beyond file length!');
end
if read_end > total_samples
    warning('Read window extends beyond file length!');
end
