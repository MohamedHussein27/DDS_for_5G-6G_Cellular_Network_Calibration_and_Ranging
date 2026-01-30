%% OFDM + LFM simultaneous transmission with AWGN only
clear; clc; close all;

%% ---------------- Parameters ----------------
NFFT = 4096;
SCS = 120e3;                 % Subcarrier spacing [Hz]
NumSCinRB = 12;
NusedFFT = 264 * NumSCinRB;  % active subcarriers (3168)
NGuardSC = NFFT - NusedFFT;
NSymbols_data = 100;         % number of data symbols
numPreamble = 3;             % number of identical preamble OFDM symbols
totalSymbols = numPreamble + NSymbols_data; % OFDM modulator NumSymbols
ModOrder = 256;              % 256-QAM

% AWGN / SNR
SNR_dB_vec = 0:2:60;   % Sweep from 0 to 60 dB in 2 dB steps

% Pilot / estimation settings
pilotSpacingSC = 8;
pilotSpacingSym = 2;
pilotBoostdB = 6;
pilotVal = 1+0j;

%% ---------------- Derived & CP ----------------
fs = NFFT * SCS;  % Sampling frequency
Ts = 1/fs;
cpLen = 128;      % choose CP length (can be adjusted)
fprintf('CP length = %d samples\n', cpLen);

%% ---------------- Frequency Plan ----------------
ofdm_bw = NusedFFT * SCS;
ofdm_edge_freq = ofdm_bw / 2;

% Guard band and LFM placement
lfm_guard = 10e6; % 10 MHz
lfm_bw = 20e6;    % 20 MHz
lfm_start_freq = ofdm_edge_freq + lfm_guard;        
lfm_stop_freq = lfm_start_freq + lfm_bw;           
lfm_center_freq = (lfm_start_freq + lfm_stop_freq) / 2;

fprintf('OFDM Band: +/- %.2f MHz\n', ofdm_edge_freq/1e6);
fprintf('LFM Band: %.2f MHz to %.2f MHz (BW=%.2f MHz)\n', ...
        lfm_start_freq/1e6, lfm_stop_freq/1e6, lfm_bw/1e6);

% RX Low-pass filter to pass OFDM and reject LFM
lpf_order = 128;
lpf_cutoff = (ofdm_edge_freq + lfm_start_freq)/2;
b_lpf = fir1(lpf_order, lpf_cutoff/(fs/2), 'low');
a_lpf = 1;

%% ---------------- OFDM objects ----------------
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
fprintf('OFDM grid: %d subcarriers x %d symbols (preamble=%d, data=%d)\n', ...
        numSC, numSym, numPreamble, NSymbols_data);

%% ---------------- Preamble ----------------
preambleRef = ones(numSC, numPreamble);

%% ---------------- Data & Pilots ----------------
msg = randi([0 ModOrder-1], numSC, NSymbols_data);
dataQAM = qammod(msg, ModOrder, 'UnitAveragePower', true);

pilotMask = false(numSC, NSymbols_data);
pilotMask(1:pilotSpacingSC:end, 1:pilotSpacingSym:end) = true;

txGrid = zeros(numSC, numSym);
txGrid(:, 1:numPreamble) = preambleRef;
txGrid(:, numPreamble+1:end) = dataQAM;

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

%% ---------------- OFDM Modulate ----------------
txSig_ofdm = ofdmmod(txGrid);
ofdm_duration_s = length(txSig_ofdm)/fs;
fprintf('Total OFDM signal length: %d samples (%.2f us)\n', length(txSig_ofdm), ofdm_duration_s*1e6);

%% ---------------- LFM Chirp & Combine ----------------
chirp_obj = dsp.Chirp('SampleRate', fs, ...
    'TargetTime', ofdm_duration_s, ...
    'InitialFrequency', lfm_start_freq, ...
    'TargetFrequency', lfm_stop_freq, ...
    'SamplesPerFrame', length(txSig_ofdm));

txSig_lfm = chirp_obj();
% Power normalization
P_ofdm = var(txSig_ofdm);
P_lfm_unscaled = var(txSig_lfm);
txSig_lfm_scaled = txSig_lfm * sqrt(P_ofdm / P_lfm_unscaled) * 0.1; % scale LFM smaller

% Combined signal (OFDM complex + LFM real)
tx_combined = txSig_ofdm + txSig_lfm_scaled;
fprintf('OFDM Power: %.2f | LFM Power: %.2f\n', P_ofdm, var(txSig_lfm_scaled));

%% ---------------- Plot Combined TX Spectrum ----------------
figure;
subplot(2,1,1);
plot(real(tx_combined(1:NFFT+cpLen)));
title('Combined TX Signal (Time Domain, 1st Symbol)');
xlabel('Samples'); ylabel('Amplitude'); grid on;

subplot(2,1,2);
pwelch(tx_combined, [], [], [], fs, 'centered', 'power');
title('Combined TX Signal (PSD)');

%% ---------------- SNR Sweep (AWGN only) ----------------
rmsEVM_vec_noEQ = zeros(size(SNR_dB_vec));
rmsEVM_vec_pilot = zeros(size(SNR_dB_vec));
rmsEVM_vec_preamble = zeros(size(SNR_dB_vec));
ber_vec_pilot = zeros(size(SNR_dB_vec));
evmObj = comm.EVM('AveragingDimensions', [1 2]);
berObj = comm.ErrorRate;

for snr_idx = 1:length(SNR_dB_vec)
    current_SNR_dB = SNR_dB_vec(snr_idx);
    
    % Add AWGN
    rxSig_noisy = awgn(tx_combined, current_SNR_dB, 'measured');
    
    % LPF for OFDM extraction
    rxSig_filtered_ofdm = filtfilt(b_lpf, a_lpf, rxSig_noisy);
    
    % OFDM demodulate
    rxGrid = ofdmdemod(rxSig_filtered_ofdm);
    
    % Separate preamble
    rxPreamble = rxGrid(:, 1:numPreamble);
    rxPreambleAvg = mean(rxPreamble, 2);
    rxDataGrid = rxGrid(:, numPreamble+1:end);
    
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
    H_est_pilot_grid = F_re(symGrid, scGrid) + 1j*F_im(symGrid, scGrid);
    
    % MMSE equalizer
    signalPower = mean(abs(txDataGrid(dataMask)).^2, 'all');
    noiseVar = signalPower / (10^(current_SNR_dB/10));
    mmseWeight = conj(H_est_pilot_grid) ./ (abs(H_est_pilot_grid).^2 + noiseVar / (signalPower+eps));
    eqGrid_pilotMMSE = mmseWeight .* rxGrid;
    eqData_pilotMMSE = eqGrid_pilotMMSE(:, numPreamble+1:end);
    
    % CPE correction (pilot-based)
    eqData_pilotMMSE_CPE = eqData_pilotMMSE;
    for symRel = 1:NSymbols_data
        if any(pilotMask(:, symRel))
            rows = find(pilotMask(:, symRel));
            colFull = symRel + numPreamble;
            receivedPilotsEq = eqGrid_pilotMMSE(rows, colFull);
            ph = angle(mean(receivedPilotsEq ./ (pilotBoost*pilotVal + eps)));
            eqData_pilotMMSE_CPE(:, symRel) = eqData_pilotMMSE_CPE(:, symRel) * exp(-1j*ph);
        end
    end
    
    % EVM & BER
    rxVec_noEQ = rxDataGrid(dataMask);
    eqPilotVec = eqData_pilotMMSE_CPE(dataMask);
    rmsEVM_vec_noEQ(snr_idx) = evmObj(refGrid_data_only, rxVec_noEQ);
    rmsEVM_vec_pilot(snr_idx) = evmObj(refGrid_data_only, eqPilotVec);
    
    rxMsg = qamdemod(eqPilotVec, ModOrder, 'UnitAveragePower', true);
    rxBits = int2bit(rxMsg, k);
    ber_results = berObj(txBits, rxBits);
    ber_vec_pilot(snr_idx) = ber_results(1);
    reset(berObj);
    
    fprintf('SNR: %2.0f dB  | Pilot EVM: %.3f %% | BER: %.2e\n', ...
            current_SNR_dB, rmsEVM_vec_pilot(snr_idx), ber_vec_pilot(snr_idx));
end
%% ---------------- NEW: RADAR SIMULATION ----------------
% (This runs once, using the noise level from the *last* SNR iteration)
fprintf('\n--- Running Radar Simulation ---\n');
c = 3e8;               % speed of light
R_true = 150000;         % target distance (m) -> increased for this Fs
alpha = 0.1;           % reflection coefficient
tau = 2*R_true/c;      % round-trip delay
delay_samples = round(tau*fs);

% Create radar return signal
rx_radar = zeros(delay_samples + length(tx_combined), 1);
rx_radar(delay_samples + (1:length(tx_combined))) = alpha*tx_combined;

% Add noise (using the last noise_power from the loop)
signal_power = mean(abs(tx_combined).^2);
SNR_linear = 10^(0/10);
noise_power = signal_power / SNR_linear;
noise_radar = sqrt(noise_power/2)*(randn(size(rx_radar)) + 1j*randn(size(rx_radar)));
rx_radar_noisy = rx_radar + noise_radar;

% Matched filter (using the full combined waveform as the template)
h_mf = conj(flipud(tx_combined));
mf_out = conv(rx_radar_noisy, h_mf, 'full');

[~, idx_max] = max(abs(mf_out));
est_delay_samples = idx_max - length(tx_combined);
R_est = c * est_delay_samples / 2 / fs;

fprintf('Target R_true = %.2f m, Estimated R_est = %.2f m\n', R_true, R_est);
% --- END NEW ---


%% ---------------- Plot EVM & BER ----------------
figure;
semilogy(SNR_dB_vec, rmsEVM_vec_noEQ, 'x-', 'LineWidth', 1.5, 'DisplayName','No EQ'); hold on;
semilogy(SNR_dB_vec, rmsEVM_vec_pilot, 'o-', 'LineWidth',1.5,'DisplayName','Pilot MMSE+CPE');
xlabel('SNR (dB)'); ylabel('RMS EVM (%)'); title('EVM vs SNR'); grid on; legend show;

figure;
semilogy(SNR_dB_vec, ber_vec_pilot, 's-', 'LineWidth',1.5); xlabel('SNR (dB)'); ylabel('BER'); title('BER vs SNR'); grid on;
