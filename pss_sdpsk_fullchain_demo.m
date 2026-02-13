function pss_sdpsk_fullchain_demo()
    close all; clc;

    %% ===================== Fixed design (your scenario) =====================
    Fs   = 409.6e6;      % capture / processing sample rate
    Rs   = 60e6;         % SDPSK symbol/chip rate
    p    = 512; q = 75;  % Fs/Rs = 512/75 exactly
    Lblk = 874;          % samples per block at Fs (after resample+trim)
    pat  = [+1 +1 -1 +1 -1 +1 -1 +1];    % [+A +A -A +A -A +A -A +A]
    Nblk = numel(pat);

    %% ===================== Channel / test knobs =====================
    pre  = 20000;        % padding before PSS (samples @ Fs)
    post = 20000;        % padding after PSS  (samples @ Fs)

    CFO_Hz   = 2000;        % inject CFO here (e.g. 2e5, 1e6)
    SNRdB    = 20;       % finite SNR to see noise in constellation
    tau_samp = 0;        % integer delay (samples @ Fs)
    phi0     = 0;        % constant phase offset (rad)

    %% ===================== TX: Step A - SDPSK @ 60MHz =====================
    s60 = gen_sdpsk_128_downlink_style();
    
    %% ===================== TX: Step B - resample to 409.6MHz =====================
    Afs = resample(s60, p, q);
    Afs = Afs(:);
    Afs = force_len(Afs, Lblk);
    Afs = Afs / rms(Afs);

    %% ===================== TX: Step C - build PSS =====================
    pss_tx = zeros(Nblk*Lblk, 1);
    for k = 1:Nblk
        pss_tx((k-1)*Lblk + (1:Lblk)) = pat(k) * Afs;
    end
    fprintf("TX PSS length = %d samples\n", length(pss_tx));

    %% ===================== Channel =====================
    y = [zeros(pre,1); pss_tx; zeros(post,1)];
    if tau_samp ~= 0, y = [zeros(tau_samp,1); y]; end
    n = (0:length(y)-1).';
    y = y .* exp(1j*(2*pi*CFO_Hz/Fs*n + phi0));
    if isfinite(SNRdB), y = awgn(y, SNRdB, 'measured'); end

    %% ===================== RX Step 1: Find PSS (Correlation) =====================
    % Use correlation to find the start of the signal
    [pk1, idx1, psr1, ~] = corr_peak_valid(y, pss_tx);
    fprintf("[RX Identification] Peak idx=%d, PSR=%.2f dB\n", idx1, psr1);

    %% ===================== RX Step 2: Extract & Downsample =====================
    % Extract the full PSS duration based on found index
    len_pss_fs = length(pss_tx);
    st_idx = idx1;
    ed_idx = st_idx + len_pss_fs - 1;
    
    % Safety check
    if st_idx < 1 || ed_idx > length(y)
        warning("Signal boundary error"); return;
    end
    
    y_raw_fs = y(st_idx:ed_idx);
    
    % Plot 1: 409.6 MHz Constellation (Time Domain Samples)
    figure; plot(real(y_raw_fs(1:2000)), imag(y_raw_fs(1:2000)), '.'); % Plot subset for speed
    axis equal; grid on; title('1. Constellation @ 409.6 MHz (Raw Samples, subset)');
    xlabel('I'); ylabel('Q');

    % Downsample to 60 MHz using Interpolation
    t_fs = 0:(length(y_raw_fs)-1);
    D_ratio = Fs / Rs;
    % For "real" simulation, we might need symbol timing recovery. 
    % Here we assume idx1 gave us good rough alignment and use offset=0.
    t_rs = 0 : D_ratio : (t_fs(end));
    
    x_60mhz = interp1(t_fs, y_raw_fs, t_rs, 'spline').';
    x_60mhz = x_60mhz / mean(abs(x_60mhz)); % Normalize

    % Plot 2: 60 MHz Constellation (Before CFO Comp)
    figure; plot(real(x_60mhz), imag(x_60mhz), '.');
    axis equal; grid on; title('2. Constellation @ 60 MHz (Before CFO Comp)');
    xlabel('I'); ylabel('Q');

    %% ===================== RX Step 3: CFO Compensation (Costas Loop) =====================
    % We process the full signal (all blocks).
    % Costas loop can generally handle 180-deg phase jumps (seen as bit flips) if bandwidth is sufficient.
    
    N_sym = length(x_60mhz);
    x_synced = zeros(size(x_60mhz));
    freq_log = zeros(1, N_sym);
    
    % Costas Parameters
    BL_T = 0.02;    % Bandwidth (tradeoff: fast lock vs noise jitter)
    zeta = 0.707;
    w_n  = BL_T / (zeta + 1/(4*zeta));
    alpha = 2 * zeta * w_n;
    beta  = w_n^2;
    
    phase_est = 0;
    freq_est  = 0;
    
    % --- Initial Frequency Acquisition (Coarse) ---
    % Use differential power spectrum or similar to get close
    % Here use a simple coarse sweep on proper metric for SDPSK (4th power diff)
    test_freqs = linspace(-20e3, 20e3, 41) / Rs * 2*pi;
    best_c0 = -inf;
    for ft = test_freqs
        z_rot = x_60mhz .* exp(-1j * ft * (1:N_sym).');
        % SDPSK metric: adjacent symbols diff should be +/- j
        dz = z_rot(2:end) .* conj(z_rot(1:end-1));
        metric = sum(abs(imag(dz)).^2 - abs(real(dz)).^2);
        if metric > best_c0
            best_c0 = metric;
            freq_est = ft;
        end
    end
    fprintf("Coarse Freq Est: %.2f Hz\n", freq_est * Rs / (2*pi));

    % Value for loop
    for n = 1:N_sym
        % 1. NCO Derotate
        z = x_60mhz(n) * exp(-1j * phase_est);
        x_synced(n) = z;
        
        % 2. Error Detector (QPSK-like for SDPSK)
        % This drives constellation to diagonals (X shape) or axes (+) depending on definition
        % SDPSK symbols are on unit circle axes if aligned? 
        % Actually SDPSK is diff encoded. The SYMBOLS themselves rotate.
        % Wait, if it is SDPSK, the raw symbols 's' are points on the circle.
        % Ideally standard Costas locks to QPSK constellation (+/-1, +/-j).
        
        I = real(z); Q = imag(z);
        
        % Standard Costas Error: Drives to axes
        err = Q * sign(I) - I * sign(Q);
        
        % 3. Filter
        freq_est = freq_est + beta * err;
        phase_est = phase_est + freq_est + alpha * err;
        
        freq_log(n) = freq_est;
    end
    
    final_cfo = mean(freq_log(end-100:end)) * Rs / (2*pi);
    fprintf("Final Loop CFO: %.2f Hz\n", final_cfo);

    % Plot 3: 60 MHz Constellation (After CFO Comp)
    figure; plot(real(x_synced), imag(x_synced), '.');
    axis equal; grid on; title(sprintf('3. Constellation @ 60 MHz (After CFO Comp)\nEst Freq: %.0f Hz', final_cfo));
    xlabel('I'); ylabel('Q');
    
    % Plot 4: Frequency Tracking
    figure; plot(freq_log * Rs / (2*pi)); grid on;
    title('Frequency Tracking'); ylabel('Hz'); xlabel('Symbol');

    fprintf("Done.\n");
end

%% ============================ SDPSK generator ============================
function s = gen_sdpsk_128_downlink_style()
    % m-seq length 127 then append first bit -> 128
    a127 = mseq_127_x7_x3_1();      % {0,1}
    a128 = [a127; a127(1)];

    % map to +/-1
    b = 2*a128 - 1;                % 1->+1, 0->-1

    % symmetric DPSK: +/-pi/2 per bit
    dphi = b * (pi/2);

    % integrate (differential encoding)
    phi = cumsum(dphi);

    % constant envelope
    s = exp(1j*phi);
end

function a = mseq_127_x7_x3_1()
    % Example primitive poly: x^7 + x^3 + 1
    reg = ones(1,7);     % seed
    a = zeros(127,1);

    for n = 1:127
        out = reg(end);
        a(n) = out;
        fb = xor(reg(end), reg(end-4));  % taps at 7 and 3
        reg = [fb, reg(1:end-1)];
    end
end

%% ============================ Correlation + PSR ============================
function [pk, idx, psr_dB, c] = corr_peak_valid(y, ref)
    c = abs(conv(y, flipud(conj(ref)), 'valid'));
    [pk, idx] = max(c);

    W = max(20, round(0.002*length(c)));
    mask = true(size(c));
    mask(max(1,idx-W):min(length(c),idx+W)) = false;

    sidelobe = c(mask);
    if isempty(sidelobe)
        psr_dB = NaN;
    else
        psr_dB = 20*log10(pk / (mean(sidelobe) + eps));
    end
end

%% ============================ Δφ scoring + shift search ============================
function [bestShift, bestScore] = best_shift_by_dphi_score(yseg, maxShift)
    bestScore = -inf;
    bestShift = 0;

    for sh = -maxShift:maxShift
        yy = apply_shift_keep_len(yseg, sh);

        mag = abs(yy);
        keep = mag > 0.1*median(mag);
        yy2 = yy(keep);
        if length(yy2) < 30
            continue;
        end

        dphi = angle(yy2(2:end).*conj(yy2(1:end-1)));
        dphi = wrapToPi(dphi);

        score = phase_concentration_score(dphi, [+pi/2, -pi/2]);

        if score > bestScore
            bestScore = score;
            bestShift = sh;
        end
    end
end

function yy = apply_shift_keep_len(y, sh)
    y = y(:);
    if sh > 0
        if sh >= length(y)
            yy = zeros(size(y));
            return;
        end
        yy = y(1+sh:end);
        yy(end+1:length(y),1) = 0;
    elseif sh < 0
        yy = [zeros(-sh,1); y];
        yy = yy(1:length(y));
    else
        yy = y;
    end
end

function score = phase_concentration_score(dphi, targets)
    dphi = wrapToPi(dphi(:));
    targets = targets(:).';

    D = zeros(length(dphi), numel(targets));
    for k = 1:numel(targets)
        D(:,k) = abs(wrapToPi(dphi - targets(k)));
    end
    dmin = min(D, [], 2);

    score = mean(exp(-(dmin/(pi/8)).^2));
end

%% ============================ Utilities ============================
function x = force_len(x, L)
    x = x(:);
    if length(x) < L
        x(end+1:L) = 0;
    else
        x = x(1:L);
    end
end

function report_sdpsk_phase_steps(tag, s)
    dphi = angle(s(2:end).*conj(s(1:end-1)));
    dphi = wrapToPi(dphi);
    fprintf("\n=== %s: Δφ stats ===\n", tag);
    fprintf("mean(|s|)=%.4f, std(|s|)=%.4g\n", mean(abs(s)), std(abs(s)));
    fprintf("Δφ mean=%.3f rad, std=%.3f rad\n", mean(dphi), std(dphi));
end
