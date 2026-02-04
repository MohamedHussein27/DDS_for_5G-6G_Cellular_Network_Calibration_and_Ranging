%clear; clc; close all;

% =========================================================
% 1. DESIGN PARAMETERS
% =========================================================
Fs       = 500e6;        % System clock (Hz)
T_dur    = 100e-6;        % Signal duration (s)
Nacc     = 32;           % Phase accumulator bits
LUT_bits = 16;           % LUT address bits
nSeeds   = 30;           % Number of random tests

Ns = round(T_dur * Fs);  % Number of samples
t  = (0:Ns-1) / Fs;      % Time vector

% Mode selection (used by your core)
current_mode = 'double';
fprintf('Testing Mode: %s\n', current_mode);

error_log = zeros(1, nSeeds);

% =========================================================
% 2. RANDOMIZED DDS TEST LOOP
% =========================================================
for seed = 1:nSeeds
    rng(seed);

    % -----------------------------------------------------
    % A. RANDOMIZED CHIRP PARAMETERS (PHYSICALLY VALID)
    % -----------------------------------------------------
    % Start frequency and bandwidth are randomized
    % Output frequency is NOT randomized directly.
    
    f0_min = 1e6;
    f0_max = 20e6;
    B_min  = 5e6;
    B_max  = 80e6;

    f0 = f0_min + (f0_max - f0_min) * rand;   % Start frequency
    B  = B_min  + (B_max  - B_min)  * rand;   % Bandwidth

    % Instantaneous frequency law (linear chirp)
    f_inst = f0 + B * (t / T_dur);

    % Safety check (Nyquist)
    assert(max(f_inst) < Fs/2, 'Chirp exceeds Nyquist frequency');

    % DDS tuning word (M)
    M = round(f_inst * (2^Nacc) / Fs);

    % -----------------------------------------------------
    % B. DDS CORE CALL
    % -----------------------------------------------------
    % DDS output is always a sine generated from a phase
    % accumulator and phase-to-amplitude LUT.

    if strcmp(current_mode, 'single') && exist('dds_core_mex', 'file') == 3
        dds_out = dds_core_mex(M, Nacc, LUT_bits, current_mode);
    else
        dds_out = dds_core(M, Nacc, LUT_bits, current_mode);
    end

    % -----------------------------------------------------
    % C. IDEAL DISCRETE-TIME REFERENCE (GROUND TRUTH)
    % -----------------------------------------------------
    % We reconstruct the *ideal* DDS accumulator using
    % infinite precision and no phase truncation.
    %
    % NOTE:
    % A chirp is simply a sinusoid with time-varying
    % instantaneous frequency, so comparing DDS output
    % (sine) to an ideal chirped sine is valid.

    ideal_phase = zeros(1, Ns);
    ideal_phase(2:end) = cumsum(double(M(1:end-1)));

    yExpected = sin(2*pi * ideal_phase / 2^Nacc);

    % -----------------------------------------------------
    % D. ERROR METRIC
    % -----------------------------------------------------
    error_log(seed) = mean(abs(double(dds_out) - yExpected));
end

% =========================================================
% 3. RESULTS VISUALIZATION
% =========================================================
figure;
plot(1:nSeeds, error_log, '.-', 'LineWidth', 1.5);
grid on;
xlabel('Test Seed');
ylabel('Mean Absolute Error');
title(['DDS Chirp Accuracy - ' current_mode ' Mode']);
subtitle('Error due to Phase Truncation and LUT Quantization Only');

fprintf('\n--- RESULTS SUMMARY ---\n');
fprintf('Average Error: %e\n', mean(error_log));
fprintf('Max Error    : %e\n', max(error_log));
fprintf('Min Error    : %e\n', min(error_log));

% DDS instantaneous frequency
f_dds = double(M) * Fs / 2^Nacc;

figure;
plot(t*1e6, f_dds/1e6, 'LineWidth', 1.5); hold on;
plot(t*1e6, f_inst/1e6, '--', 'LineWidth', 1.5);
grid on;
xlabel('Time (\mus)');
ylabel('Frequency (MHz)');
legend('DDS Output Frequency', 'Target Chirp Frequency');
title('DDS Instantaneous Frequency vs Time');
