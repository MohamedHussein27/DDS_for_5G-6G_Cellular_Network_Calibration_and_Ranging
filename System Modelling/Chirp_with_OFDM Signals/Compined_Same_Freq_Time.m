%% Combined OFDM + LFM Signal Transmission with AWGN Only
clear; clc; close all;

%% ---------------- OFDM Parameters ----------------
NFFT = 4096;
SCS = 120e3;                 % Subcarrier spacing [Hz]
NumSCinRB = 12;
NusedFFT = 264 * NumSCinRB;  % active subcarriers
NGuardSC = NFFT - NusedFFT;
NSymbols_data = 100;         % number of data symbols
numPreamble = 3;             % number of identical preamble OFDM symbols
totalSymbols = numPreamble + NSymbols_data;
ModOrder = 256;              % 256-QAM

% Pilot / estimation settings
pilotSpacingSC = 8;
pilotSpacingSym = 2;
pilotBoostdB = 6;
pilotVal = 1+0j;

%% ---------------- LFM Parameters ----------------
c = 3e8;                     % speed of light (m/s)
A_lfm = 1;                   % amplitude of LFM chirp
f0_lfm = 0;                  % start frequency (Hz)
B_lfm = 80e6;                % LFM bandwidth (Hz) - 80 MHz as requested
T_lfm = 50e-6;               % LFM pulse duration (seconds)
K_lfm = B_lfm / T_lfm;       % chirp rate (Hz/s)
theta_lfm = 0;               % initial phase (rad)
R_true = 302;                % true target range in meters
alpha_lfm = 0.3;             % target amplitude attenuation

%% ---------------- Derived Parameters ----------------
fs = NFFT * SCS;             % Sampling frequency (same for both signals)
Ts = 1/fs;
cpLen = 20;                  % Simplified CP length for AWGN only

fprintf('Sampling frequency fs = %.2f MHz\n', fs/1e6);
fprintf('CP Length = %d samples\n', cpLen);

%% ---------------- SNR Sweep ----------------
SNR_dB_vec = 0:2:50;

%% ---------------- OFDM Setup ----------------
ofdmmod = comm.OFDMModulator('FFTLength', NFFT, ...
    'NumSymbols', totalSymbols, ...
    'NumGuardBandCarriers', [NGuardSC/2; NGuardSC/2], ...
    'CyclicPrefixLength', cpLen, ...
    'InsertDCNull', false);
ofdmdemod = comm.OFDMDemodulator('FFTLength', NFFT, ...
    'NumSymbols', totalSymbols, ...
    'NumGuardBandCarriers', [NGuardSC/2; NGuardSC/2], ...
    'CyclicPrefixLength', cpLen, ...
    'RemoveDCCarrier', false);
ofdmInfo = info(ofdmmod);
numSC = ofdmInfo.DataInputSize(1);
numSym = ofdmInfo.DataInputSize(2);

fprintf('OFDM grid: %d subcarriers x %d symbols\n', numSC, numSym);

%% ---------------- OFDM Data Generation ----------------
rng(42);
msg = randi([0 ModOrder-1], numSC, NSymbols_data);
dataQAM = qammod(msg, ModOrder, 'UnitAveragePower', true);

% Pilot mask
pilotMask = false(numSC, NSymbols_data);
pilotMask(1:pilotSpacingSC:end, 1:pilotSpacingSym:end) = true;

% Build OFDM grid
txGrid = zeros(numSC, numSym);
preambleRef = ones(numSC, numPreamble);
txGrid(:, 1:numPreamble) = preambleRef;
txGrid(:, numPreamble+1:end) = dataQAM;

% Insert boosted pilots
[pRows_rel, pCols_rel] = find(pilotMask);
pCols_full = pCols_rel + numPreamble;
pilotBoost = 10^(pilotBoostdB/20);
linearIdxPilots = sub2ind(size(txGrid), pRows_rel, pCols_full);
txGrid(linearIdxPilots) = pilotBoost * pilotVal;

txDataGrid = txGrid(:, numPreamble+1:end);
dataMask = ~pilotMask;
refGrid_data_only = txDataGrid(dataMask);
msg_data_only = msg(dataMask);
k = log2(ModOrder);
txBits = int2bit(msg_data_only, k);

%% ---------------- OFDM Modulation ----------------
txSig_OFDM = ofdmmod(txGrid);
samples_per_symbol = NFFT + cpLen;

fprintf('OFDM signal length: %d samples (%.3f ms)\n', length(txSig_OFDM), length(txSig_OFDM)/fs*1e3);

%% ---------------- LFM Signal Generation ----------------
% Generate LFM chirp with proper pulse duration T_lfm
t_lfm = (0:round(T_lfm*fs)-1) / fs;  % Only generate for pulse duration
s_lfm_complex = A_lfm * exp(1j*(2*pi*(f0_lfm.*t_lfm + 0.5*K_lfm.*t_lfm.^2) + theta_lfm));
s_lfm_complex = s_lfm_complex(:); % Column vector

% Create delayed echo for radar simulation
tau = 2*R_true / c;
delay_samples = round(tau * fs);
fprintf('LFM delay samples: %d (delay = %.3e s)\n', delay_samples, tau);

% Create received LFM with delay - pad with zeros to match OFDM length
s_lfm_received = [zeros(delay_samples, 1); alpha_lfm * s_lfm_complex; zeros(length(txSig_OFDM) - delay_samples - length(s_lfm_complex), 1)];
if length(s_lfm_received) > length(txSig_OFDM)
    s_lfm_received = s_lfm_received(1:length(txSig_OFDM));
end

% Pad transmitted LFM to match OFDM length
s_lfm_complex_padded = [s_lfm_complex; zeros(length(txSig_OFDM) - length(s_lfm_complex), 1)];

fprintf('LFM pulse duration: %.1f µs\n', T_lfm*1e6);
fprintf('LFM chirp samples: %d\n', length(s_lfm_complex));
fprintf('LFM bandwidth: %.1f MHz\n', B_lfm/1e6);
fprintf('LFM chirp rate: %.3f MHz/µs\n', K_lfm/1e12);

%% ---------------- Power Normalization ----------------
% Increase OFDM power by factor of 2 (3 dB)
OFDM_power_boost = 10;  % Adjust this value to increase/decrease OFDM power
txSig_OFDM = txSig_OFDM * sqrt(OFDM_power_boost);

% Normalize signals to have equal power
P_ofdm = mean(abs(txSig_OFDM).^2);
P_lfm = mean(abs(s_lfm_complex_padded).^2);
fprintf('OFDM power: %.4f, LFM power: %.4f\n', P_ofdm, P_lfm);
fprintf('OFDM power boost factor: %.2f (%.2f dB)\n', OFDM_power_boost, 10*log10(OFDM_power_boost));

% Normalize LFM to match OFDM power
s_lfm_complex_padded = s_lfm_complex_padded * sqrt(P_ofdm / P_lfm);
s_lfm_received = s_lfm_received * sqrt(P_ofdm / P_lfm);

%% ---------------- Combined Signal Transmission ----------------
% Transmitted signal = OFDM + LFM (both in same frequency band)
txSig_Combined = txSig_OFDM + s_lfm_complex_padded;

% Received signal = OFDM + LFM_echo (both in same frequency band)
rxSig_Combined_base = txSig_OFDM + s_lfm_received;

fprintf('\nCombined signal created: OFDM + LFM in same frequency band\n');
fprintf('Combined signal length: %d samples\n', length(txSig_Combined));

%% ---------------- Enhanced LFM Analysis ----------------
% Analyze LFM frequency characteristics
N_lfm = length(s_lfm_complex);
f_lfm_axis = (-N_lfm/2:N_lfm/2-1)*(fs/N_lfm);
S_lfm = fftshift(fft(s_lfm_complex));

% Calculate theoretical LFM bandwidth
f_inst = f0_lfm + K_lfm * t_lfm;
actual_B_lfm = max(f_inst) - min(f_inst);
fprintf('\nLFM Bandwidth Analysis:\n');
fprintf('Theoretical bandwidth: %.1f MHz\n', B_lfm/1e6);
fprintf('Actual instantaneous bandwidth: %.1f MHz\n', actual_B_lfm/1e6);

%% ---------------- SNR Sweep Simulation ----------------
fprintf('\n=== Starting SNR Sweep (AWGN Only) ===\n');

% Pre-allocate results
rmsEVM_vec_pilot = zeros(size(SNR_dB_vec));
rmsEVM_vec_preamble = zeros(size(SNR_dB_vec));
rmsEVM_vec_noEQ = zeros(size(SNR_dB_vec));
ber_vec_pilot = zeros(size(SNR_dB_vec));
range_error_vec = zeros(size(SNR_dB_vec));

evmObj = comm.EVM('AveragingDimensions', [1 2]);
berObj = comm.ErrorRate;

for snr_idx = 1:length(SNR_dB_vec)
    current_SNR_dB = SNR_dB_vec(snr_idx);
    
    % Add AWGN to combined signal
    rxSig_Combined = awgn(rxSig_Combined_base, current_SNR_dB, 'measured');
    
    %% -------- OFDM Processing --------
    % Demodulate OFDM (note: LFM acts as interference here)
    rxGrid = ofdmdemod(rxSig_Combined);
    rxPreamble = rxGrid(:, 1:numPreamble);
    rxPreambleAvg = mean(rxPreamble, 2);
    rxDataGrid = rxGrid(:, numPreamble+1:end);
    
    % Estimate noise variance
    signalPower = mean(abs(txDataGrid(dataMask)).^2, 'all');
    noiseVar = signalPower / (10^(current_SNR_dB/10));
    
    % Channel estimation and equalization
    [H_est_pilot_grid, eqPilotVec, eqPreambleVec, rxVec_noEQ] = ...
        channelEstimationAndEqualization(rxGrid, rxDataGrid, txGrid, dataMask, ...
                                       pilotMask, numPreamble, preambleRef, ...
                                       rxPreambleAvg, pilotBoost, pilotVal, noiseVar, fs);
    
    % EVM Measurement
    rmsEVM_vec_noEQ(snr_idx) = evmObj(refGrid_data_only, rxVec_noEQ);
    rmsEVM_vec_pilot(snr_idx) = evmObj(refGrid_data_only, eqPilotVec);
    rmsEVM_vec_preamble(snr_idx) = evmObj(refGrid_data_only, eqPreambleVec);
    
    % BER Calculation
    rxMsg = qamdemod(eqPilotVec, ModOrder, 'UnitAveragePower', true);
    rxBits = int2bit(rxMsg, k);
    ber_results = berObj(txBits, rxBits);
    ber_vec_pilot(snr_idx) = ber_results(1);
    reset(berObj);
    
    %% -------- LFM Processing (Matched Filtering) --------
    % Use the original LFM pulse for matched filtering
    h_lfm = conj(fliplr(s_lfm_complex.'));
    z_lfm = conv(rxSig_Combined, h_lfm);
    
    % Find peak
    [~, idx_max] = max(abs(z_lfm));
    est_delay_samples = idx_max - length(s_lfm_complex);
    est_tau = est_delay_samples / fs;
    R_est = c * est_tau / 2;
    range_error_vec(snr_idx) = abs(R_est - R_true);
    
    % Store last iteration for plotting
    if snr_idx == length(SNR_dB_vec)
        rxSig_final = rxSig_Combined;
        rxGrid_final = rxGrid;
        H_est_final = H_est_pilot_grid;
        eqPilotVec_final = eqPilotVec;
        eqPreambleVec_final = eqPreambleVec;
        rxVec_noEQ_final = rxVec_noEQ;
        z_lfm_final = z_lfm;
        idx_max_final = idx_max;
        R_est_final = R_est;
    end
    
    if mod(snr_idx, 5) == 0
        fprintf('SNR = %d dB: EVM = %.2f%%, BER = %.2e, Range Error = %.2f m\n', ...
                current_SNR_dB, rmsEVM_vec_pilot(snr_idx), ber_vec_pilot(snr_idx), range_error_vec(snr_idx));
    end
end

%% ---------------- Calculate PAPR ----------------
PAPR_tx_ofdm = calculatePAPR(txSig_OFDM);
PAPR_tx_lfm = calculatePAPR(s_lfm_complex);
PAPR_tx_combined = calculatePAPR(txSig_Combined);
PAPR_rx = calculatePAPR(rxSig_final);

fprintf('\nPAPR Results:\n');
fprintf('  OFDM alone: %.2f dB\n', PAPR_tx_ofdm);
fprintf('  LFM alone: %.2f dB\n', PAPR_tx_lfm);
fprintf('  Combined TX: %.2f dB\n', PAPR_tx_combined);
fprintf('  Combined RX: %.2f dB\n', PAPR_rx);

%% ---------------- Plotting ----------------

% Figure 1: Time Domain Signals
figure(1);
subplot(3,1,1);
plot(real(txSig_OFDM(1:3*samples_per_symbol)));
title('Transmitted OFDM Signal (First 3 Symbols)');
xlabel('Samples'); ylabel('Amplitude'); grid on;

subplot(3,1,2);
plot(real(s_lfm_complex_padded(1:3*samples_per_symbol)));
title('Transmitted LFM Chirp Signal (First 3 OFDM Symbol Durations)');
xlabel('Samples'); ylabel('Amplitude'); grid on;

subplot(3,1,3);
plot(real(txSig_Combined(1:3*samples_per_symbol)));
title('Combined Signal: OFDM + LFM (First 3 Symbols)');
xlabel('Samples'); ylabel('Amplitude'); grid on;

%% ---------------- REPLACED: PSD Analysis using Welch Method ----------------
% Figure 2: Power Spectral Density (Welch) - Individual and Combined
figure(2);
set(gcf, 'Position', [100 100 1200 900]);

% PSD parameters for consistent analysis
window_length = 1024;
overlap = 512;
nfft_psd = 4096;

% Calculate PSD for all signals
[psd_ofdm, f_ofdm] = pwelch(txSig_OFDM, window_length, overlap, nfft_psd, fs, 'centered');
[psd_lfm, f_lfm] = pwelch(s_lfm_complex_padded, window_length, overlap, nfft_psd, fs, 'centered');
[psd_combined, f_combined] = pwelch(txSig_Combined, window_length, overlap, nfft_psd, fs, 'centered');

% Plot 1: OFDM PSD
subplot(3,1,1);
plot(f_ofdm/1e6, 10*log10(psd_ofdm), 'b', 'LineWidth', 1.5);
xlabel('Frequency (MHz)'); ylabel('PSD (dB/Hz)');
title('Power Spectral Density: OFDM Signal Only (Welch Method)');
grid on;
xlim([-fs/2e6*0.6, fs/2e6*0.6]); % Focus on central frequencies
% Add OFDM bandwidth annotation
yl = ylim;
hold on;
ofdm_bw = (numSC * SCS) / 1e6;
plot([-ofdm_bw/2, ofdm_bw/2], [yl(1)+0.9*diff(yl), yl(1)+0.9*diff(yl)], 'r--', 'LineWidth', 2);
text(0, yl(1)+0.85*diff(yl), sprintf('OFDM BW = %.1f MHz', ofdm_bw), ...
     'HorizontalAlignment', 'center', 'BackgroundColor', 'white');
hold off;

% Plot 2: LFM PSD
subplot(3,1,2);
plot(f_lfm/1e6, 10*log10(psd_lfm), 'r', 'LineWidth', 1.5);
xlabel('Frequency (MHz)'); ylabel('PSD (dB/Hz)');
title('Power Spectral Density: LFM Signal Only (Welch Method)');
grid on;
xlim([-B_lfm*1.5/1e6, B_lfm*1.5/1e6]); % Focus on LFM bandwidth
% Add LFM bandwidth annotation
yl = ylim;
hold on;
plot([-B_lfm/2e6, B_lfm/2e6], [yl(1)+0.9*diff(yl), yl(1)+0.9*diff(yl)], 'k--', 'LineWidth', 2);
text(0, yl(1)+0.85*diff(yl), sprintf('LFM Bandwidth = %.1f MHz', B_lfm/1e6), ...
     'HorizontalAlignment', 'center', 'BackgroundColor', 'white');
hold off;

% Plot 3: Combined PSD
subplot(3,1,3);
plot(f_combined/1e6, 10*log10(psd_combined), 'k', 'LineWidth', 1.5);
hold on;
plot(f_ofdm/1e6, 10*log10(psd_ofdm), 'b--', 'LineWidth', 1, 'DisplayName', 'OFDM PSD');
plot(f_lfm/1e6, 10*log10(psd_lfm), 'r--', 'LineWidth', 1, 'DisplayName', 'LFM PSD');
hold off;
xlabel('Frequency (MHz)'); ylabel('PSD (dB/Hz)');
title('Power Spectral Density: Combined Signal (OFDM + LFM) - Welch Method');
legend('Combined', 'OFDM', 'LFM', 'Location', 'best');
grid on;
xlim([-fs/2e6*0.6, fs/2e6*0.6]);

% NEW Figure: Detailed LFM Analysis
figure(11);
set(gcf, 'Position', [100 100 1200 800]);

% Plot 1: LFM Time Domain
subplot(2,3,1);
plot(t_lfm*1e6, real(s_lfm_complex));
xlabel('Time (µs)'); ylabel('Amplitude');
title('LFM Chirp - Time Domain');
grid on;

% Plot 2: LFM Instantaneous Frequency
subplot(2,3,2);
f_inst = f0_lfm + K_lfm * t_lfm;
plot(t_lfm*1e6, f_inst/1e6);
xlabel('Time (µs)'); ylabel('Frequency (MHz)');
title('LFM Instantaneous Frequency');
grid on;
% Add bandwidth annotation
yl = ylim;
hold on;
plot([0, T_lfm*1e6], [f0_lfm/1e6, (f0_lfm+B_lfm)/1e6], 'r--', 'LineWidth', 1);
text(T_lfm*1e6/2, mean(yl), sprintf('Δf = %.1f MHz', B_lfm/1e6), ...
     'HorizontalAlignment', 'center', 'BackgroundColor', 'white');
hold off;

% Plot 3: LFM Frequency Spectrum (FFT for comparison)
subplot(2,3,3);
plot(f_lfm_axis/1e6, abs(S_lfm));
xlabel('Frequency (MHz)'); ylabel('Magnitude');
title('LFM Frequency Spectrum (FFT)');
grid on;
xlim([-B_lfm*1.2/1e6, B_lfm*1.2/1e6]);
% Add bandwidth annotation
yl = ylim;
hold on;
plot([-B_lfm/2e6, B_lfm/2e6], [yl(1)+0.9*diff(yl), yl(1)+0.9*diff(yl)], 'r--', 'LineWidth', 2);
text(0, yl(1)+0.85*diff(yl), sprintf('Bandwidth = %.1f MHz', B_lfm/1e6), ...
     'HorizontalAlignment', 'center', 'BackgroundColor', 'white');
hold off;

% Plot 4: LFM Phase
subplot(2,3,4);
phase_unwrapped = unwrap(angle(s_lfm_complex));
plot(t_lfm*1e6, phase_unwrapped);
xlabel('Time (µs)'); ylabel('Phase (rad)');
title('LFM Unwrapped Phase');
grid on;

% Plot 5: LFM Spectrogram
subplot(2,3,5);
spectrogram(s_lfm_complex, 256, 200, 1024, fs, 'yaxis');
title('LFM Spectrogram');
colorbar;

% Plot 6: LFM PSD using Welch method (high resolution)
subplot(2,3,6);
[psd_lfm_detailed, f_lfm_detailed] = pwelch(s_lfm_complex, 2048, 1024, 4096, fs, 'centered');
plot(f_lfm_detailed/1e6, 10*log10(psd_lfm_detailed), 'r', 'LineWidth', 1.5);
xlabel('Frequency (MHz)'); ylabel('PSD (dB/Hz)');
title('LFM Power Spectral Density (Welch - High Res)');
grid on;
xlim([-B_lfm*1.2/1e6, B_lfm*1.2/1e6]);
% Add bandwidth annotation
yl = ylim;
hold on;
plot([-B_lfm/2e6, B_lfm/2e6], [yl(1)+0.9*diff(yl), yl(1)+0.9*diff(yl)], 'k--', 'LineWidth', 2);
text(0, yl(1)+0.85*diff(yl), sprintf('80 MHz Bandwidth', B_lfm/1e6), ...
     'HorizontalAlignment', 'center', 'BackgroundColor', 'white');
hold off;

% Continue with the rest of your figures (3-10)...
% Figure 3: OFDM Performance - EVM vs SNR
figure(3);
semilogy(SNR_dB_vec, rmsEVM_vec_pilot, 'b-o', 'LineWidth', 1.5, 'DisplayName', 'Pilot MMSE');
hold on;
semilogy(SNR_dB_vec, rmsEVM_vec_preamble, 'r-s', 'LineWidth', 1.5, 'DisplayName', 'Preamble MMSE');
semilogy(SNR_dB_vec, rmsEVM_vec_noEQ, 'g-^', 'LineWidth', 1.5, 'DisplayName', 'No EQ');
hold off;
xlabel('SNR (dB)'); ylabel('RMS EVM (%)');
title('OFDM Performance: EVM vs SNR (with LFM interference)');
legend('show'); grid on;

% Figure 4: OFDM Performance - BER vs SNR
figure(4);
semilogy(SNR_dB_vec, ber_vec_pilot, 'b-o', 'LineWidth', 1.5);
xlabel('SNR (dB)'); ylabel('Bit Error Rate');
title('OFDM Performance: BER vs SNR (with LFM interference)');
grid on;

% Figure 5: LFM Performance - Range Error vs SNR
figure(5);
semilogy(SNR_dB_vec, range_error_vec, 'r-o', 'LineWidth', 1.5);
xlabel('SNR (dB)'); ylabel('Range Error (m)');
title('LFM Radar Performance: Range Error vs SNR (with OFDM interference)');
grid on;

% Figure 6: Constellation Diagrams (at highest SNR)
figure(6);
nplot = min(6000, numel(rxVec_noEQ_final));
sel = randperm(numel(rxVec_noEQ_final), nplot);

subplot(1,3,1);
plot(real(rxVec_noEQ_final(sel)), imag(rxVec_noEQ_final(sel)), '.');
title(sprintf('No EQ\nEVM: %.2f%%', rmsEVM_vec_noEQ(end)));
axis square; grid on;

subplot(1,3,2);
plot(real(eqPilotVec_final(sel)), imag(eqPilotVec_final(sel)), '.');
title(sprintf('Pilot MMSE\nEVM: %.2f%%', rmsEVM_vec_pilot(end)));
axis square; grid on;

subplot(1,3,3);
plot(real(eqPreambleVec_final(sel)), imag(eqPreambleVec_final(sel)), '.');
title(sprintf('Preamble MMSE\nEVM: %.2f%%', rmsEVM_vec_preamble(end)));
axis square; grid on;

% Figure 7: LFM Matched Filter Output
figure(7);
t_z = (0:length(z_lfm_final)-1) / fs;
plot(t_z*1e6, abs(z_lfm_final));
hold on;
plot(t_z(idx_max_final)*1e6, abs(z_lfm_final(idx_max_final)), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
hold off;
xlabel('Time (μs)'); ylabel('|Matched Filter Output|');
title(sprintf('LFM Matched Filter Output (Estimated Range: %.2f m, True: %.2f m)', R_est_final, R_true));
grid on;
legend('Matched Filter Output', 'Peak', 'Location', 'best');

% Figure 8: PAPR Comparison
figure(8);
bar_data = [PAPR_tx_ofdm, PAPR_tx_lfm, PAPR_tx_combined, PAPR_rx];
bar(bar_data);
set(gca, 'XTickLabel', {'OFDM TX', 'LFM TX', 'Combined TX', 'Combined RX'});
ylabel('PAPR (dB)');
title('PAPR Comparison');
grid on;
for i = 1:4
    text(i, bar_data(i)+0.1, sprintf('%.2f dB', bar_data(i)), 'HorizontalAlignment', 'center');
end

% Figure 9: Spectrogram Comparison
figure(9);
subplot(2,2,1);
spectrogram(txSig_OFDM, 256, 200, 1024, fs, 'yaxis');
title('Spectrogram: OFDM Signal');
colorbar;

subplot(2,2,2);
spectrogram(s_lfm_complex, 256, 200, 1024, fs, 'yaxis');
title('Spectrogram: LFM Signal');
colorbar;

subplot(2,2,3);
spectrogram(txSig_Combined, 256, 200, 1024, fs, 'yaxis');
title('Spectrogram: Combined Signal');
colorbar;

subplot(2,2,4);
spectrogram(rxSig_final, 256, 200, 1024, fs, 'yaxis');
title('Spectrogram: Received Combined Signal');
colorbar;

fprintf('\n=== Simulation Complete! ===\n');
fprintf('Generated 11 figures showing combined OFDM + LFM performance\n');
fprintf('Figure 2: PSD Analysis using Welch Method (Replaced FFT)\n');
fprintf('Figure 11: Detailed LFM Analysis\n');

%% ---------------- Support Functions ----------------
function papr = calculatePAPR(signal)
    power = abs(signal).^2;
    peak_power = max(power);
    avg_power = mean(power);
    papr = 10*log10(peak_power / avg_power);
end

function [H_est_pilot_grid, eqPilotVec, eqPreambleVec, rxVec_noEQ] = ...
    channelEstimationAndEqualization(rxGrid, rxDataGrid, txGrid, dataMask, ...
                                   pilotMask, numPreamble, preambleRef, ...
                                   rxPreambleAvg, pilotBoost, pilotVal, noiseVar, fs)
    
    % [Your existing function code remains the same...]
    numSC = size(rxGrid, 1);
    numSym = size(rxGrid, 2);
    NSymbols_data = size(rxDataGrid, 2);
    
    % Pilot-based LS & interpolation
    [pRows_rel, pCols_rel] = find(pilotMask);
    pCols_full = pCols_rel + numPreamble;
    Yp = rxGrid(sub2ind(size(rxGrid), pRows_rel, pCols_full));
    Xp = pilotBoost * pilotVal * ones(size(Yp));
    Hp = Yp ./ (Xp + eps);
    
    xPil = double(pCols_full(:));
    yPil = double(pRows_rel(:));
    pvals_re = real(Hp(:));
    pvals_im = imag(Hp(:));
    
    [symGrid, scGrid] = meshgrid(1:numSym, 1:numSC);
    F_re = scatteredInterpolant(xPil, yPil, pvals_re, 'linear', 'nearest');
    F_im = scatteredInterpolant(xPil, yPil, pvals_im, 'linear', 'nearest');
    
    H_est_pilot_grid = F_re(symGrid, scGrid) + 1j * F_im(symGrid, scGrid);
    H_est_pilot_grid = fillmissing(H_est_pilot_grid, 'nearest');
    
    % Smoothing
    k_smooth = [1 4 6 4 1]; k_smooth = k_smooth/sum(k_smooth);
    H_est_pilot_grid = conv2(k_smooth, k_smooth, H_est_pilot_grid, 'same');
    
    % Signal power for MMSE
    txDataGrid = txGrid(:, numPreamble+1:end);
    signalPower = mean(abs(txDataGrid(dataMask)).^2, 'all');
    
    % MMSE equalizer (pilot-based)
    H_use = H_est_pilot_grid;
    mmseWeight = conj(H_use) ./ (abs(H_use).^2 + noiseVar / (signalPower + eps));
    eqGrid_pilotMMSE = mmseWeight .* rxGrid;
    eqData_pilotMMSE = eqGrid_pilotMMSE(:, numPreamble+1:end);
    
    % CPE correction (pilot-based)
    eqData_pilotMMSE_CPE = eqData_pilotMMSE;
    for symRel = 1:NSymbols_data
        if any(pilotMask(:, symRel))
            rows = find(pilotMask(:, symRel));
            colFull = symRel + numPreamble;
            receivedPilotsEq = eqGrid_pilotMMSE(rows, colFull);
            ph = angle(mean(receivedPilotsEq ./ (pilotBoost * pilotVal + eps)));
            eqData_pilotMMSE_CPE(:, symRel) = eqData_pilotMMSE_CPE(:, symRel) * exp(-1j * ph);
        end
    end
    
    % Preamble-based MMSE
    H_preamble = rxPreambleAvg ./ (mean(preambleRef,2) + eps);
    H_preamble_grid = repmat(H_preamble, 1, numSym);
    mmseWeight_preamble = conj(H_preamble_grid) ./ (abs(H_preamble_grid).^2 + noiseVar / (signalPower + eps));
    eqGrid_preambleMMSE = mmseWeight_preamble .* rxGrid;
    eqData_preambleMMSE = eqGrid_preambleMMSE(:, numPreamble+1:end);
    
    % CPE correction (preamble-based)
    eqData_preambleMMSE_CPE = eqData_preambleMMSE;
    for symRel = 1:NSymbols_data
        if any(pilotMask(:, symRel))
            rows = find(pilotMask(:, symRel));
            colFull = symRel + numPreamble;
            receivedPilotsEq = eqGrid_preambleMMSE(rows, colFull);
            ph = angle(mean(receivedPilotsEq ./ (pilotBoost * pilotVal + eps)));
            eqData_preambleMMSE_CPE(:, symRel) = eqData_preambleMMSE_CPE(:, symRel) * exp(-1j * ph);
        end
    end
    
    % Extract data vectors
    rxVec_noEQ = rxDataGrid(dataMask);
    eqPilotVec = eqData_pilotMMSE_CPE(dataMask);
    eqPreambleVec = eqData_preambleMMSE_CPE(dataMask);
end