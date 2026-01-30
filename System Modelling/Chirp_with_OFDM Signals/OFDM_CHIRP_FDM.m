%% OFDM + LFM simultaneous transmission with AWGN only (preamble-only)
clear; clc; close all;

%% ---------------- Parameters ----------------
NFFT = 4096;
SCS = 120e3;                 
NumSCinRB = 12;
NusedFFT = 264 * NumSCinRB;  
NGuardSC = NFFT - NusedFFT;
NSymbols_data = 100;         
numPreamble = 3;             
totalSymbols = numPreamble + NSymbols_data; 
ModOrder = 256;              
SNR_dB_vec = 0:2:40;  

%% ---------------- Derived & CP ----------------
fs = NFFT * SCS;  
Ts = 1/fs;
cpLen = 128;      
fprintf('CP length = %d samples\n', cpLen);

%% ---------------- Frequency Plan ----------------
ofdm_bw = NusedFFT * SCS;
ofdm_edge_freq = ofdm_bw / 2;

lfm_guard = 10e6; 
lfm_bw = 80e6;    
lfm_start_freq = ofdm_edge_freq + lfm_guard;        
lfm_stop_freq = lfm_start_freq + lfm_bw;           

fprintf('OFDM Band: +/- %.2f MHz\n', ofdm_edge_freq/1e6);
fprintf('LFM Band: %.2f MHz to %.2f MHz (BW=%.2f MHz)\n', ...
        lfm_start_freq/1e6, lfm_stop_freq/1e6, lfm_bw/1e6);

lpf_order = 128;
lpf_cutoff = (ofdm_edge_freq + lfm_start_freq)/2;
b_lpf = fir1(lpf_order, lpf_cutoff/(fs/2), 'low');
a_lpf = 1;

hp_order = 128;
hp_cutoff = (ofdm_edge_freq + lfm_start_freq)/2;
b_hp = fir1(hp_order, hp_cutoff/(fs/2), 'high');
a_hp = 1;

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

%% ---------------- Preamble ----------------
preambleRef = ones(numSC, numPreamble);

%% ---------------- Data ----------------
msg = randi([0 ModOrder-1], numSC, NSymbols_data);
dataQAM = qammod(msg, ModOrder, 'UnitAveragePower', true);

txGrid = zeros(numSC, numSym);
txGrid(:, 1:numPreamble) = preambleRef;
txGrid(:, numPreamble+1:end) = dataQAM;

%% ---------------- OFDM Modulate ----------------
txSig_ofdm = ofdmmod(txGrid);
ofdm_duration_s = length(txSig_ofdm)/fs;

%% ---------------- LFM Chirp Generation (Baseband) ----------------
t_vec = (0:length(txSig_ofdm)-1)'/fs;

% Baseband LFM: from -BW/2 to +BW/2 around carrier
lfm_baseband = exp(1j*2*pi*((-lfm_bw/2)*t_vec + (lfm_bw/(2*ofdm_duration_s))*t_vec.^2));

% Upconvert to LFM center frequency
lfm_center_freq = (lfm_start_freq + lfm_stop_freq)/2;
txSig_lfm = lfm_baseband .* exp(1j*2*pi*lfm_center_freq*t_vec);

% Power normalization: LFM = OFDM
P_ofdm = mean(abs(txSig_ofdm).^2);
P_lfm_unscaled = mean(abs(txSig_lfm).^2);
txSig_lfm_scaled = txSig_lfm * sqrt(P_ofdm / P_lfm_unscaled);

% Combined signal
tx_combined = txSig_ofdm + txSig_lfm_scaled;
fprintf('OFDM Power: %.2f | LFM Power: %.2f\n', P_ofdm, mean(abs(txSig_lfm_scaled).^2));

%% ---------------- Transmit Signal Plots ----------------
figure;
subplot(2,1,1);
plot(real(tx_combined(1:min(length(tx_combined), NFFT+cpLen+1000))));
title('Combined TX Signal (Time Domain, 1st Symbol)');
xlabel('Samples'); ylabel('Amplitude'); grid on;

subplot(2,1,2);
[pxx,freq] = pwelch(tx_combined, 4096, 2048, 4096, fs, 'centered');
plot(freq/1e6, 10*log10(pxx)); grid on;
xlabel('Frequency (MHz)'); ylabel('Power/Frequency (dB/Hz)');
title('Combined TX Signal PSD');

%% ---------------- SNR Sweep (AWGN only) ----------------
rmsEVM_preamble = zeros(size(SNR_dB_vec));
evmObj = comm.EVM('AveragingDimensions',[1 2]);

for snr_idx = 1:length(SNR_dB_vec)
    SNRdB = SNR_dB_vec(snr_idx);
    rxSig = awgn(tx_combined, SNRdB, 'measured');
    rxSig_filtered = filtfilt(b_lpf, a_lpf, rxSig);
    rxGrid = ofdmdemod(rxSig_filtered);
    rxPreamble = rxGrid(:,1:numPreamble);
    H_est = mean(rxPreamble ./ preambleRef, 2);
    eqPreamble = rxPreamble ./ H_est;
    for symIdx = 1:numPreamble
        eqPreamble(:,symIdx) = eqPreamble(:,symIdx) * exp(-1j*angle(mean(eqPreamble(:,symIdx))));
    end
    rmsEVM_preamble(snr_idx) = evmObj(preambleRef, eqPreamble);
    fprintf('SNR: %2.0f dB | Preamble EVM: %.3f %%\n', SNRdB, rmsEVM_preamble(snr_idx));
end

figure;
semilogy(SNR_dB_vec, rmsEVM_preamble,'s-','LineWidth',1.5);
xlabel('SNR (dB)'); ylabel('RMS EVM (%)');
title('Preamble-only MMSE + CPE EVM vs SNR'); grid on;

%% ---------------- Corrected RADAR SIMULATION ----------------
c = 3e8;
R_true = 302;           % True range
alpha = 0.5;            % Target reflection factor
tau_true = 2*R_true/c;  % Round-trip delay
delay_samples = round(tau_true*fs);

% Received radar signal with delay
rx_radar = zeros(delay_samples + length(txSig_lfm_scaled), 1);
rx_radar(delay_samples + (1:length(txSig_lfm_scaled))) = alpha*txSig_lfm_scaled;

% Add AWGN
signal_power = mean(abs(txSig_lfm_scaled).^2);
SNR_dB_fixed = 30;
SNR_linear = 10^(SNR_dB_fixed/10);
noise_power = signal_power / SNR_linear;
noise_radar = sqrt(noise_power/2)*(randn(size(rx_radar)) + 1j*randn(size(rx_radar)));
rx_radar_noisy = rx_radar + noise_radar;

%% ---------------- High-pass Filter to Extract LFM ----------------
rx_radar_lfm = filtfilt(b_hp, a_hp, rx_radar_noisy);

%% ---------------- Low-pass Filter to Extract OFDM ----------------
rxSig_ofdm_only = filtfilt(b_lpf, a_lpf, rx_radar_noisy);

%% ---------------- Correlation-based LFM Detection ----------------
% Use cross-correlation instead of manual matched filter
[corr_out, lags] = xcorr(rx_radar_lfm, txSig_lfm_scaled);

[~, idx_peak_corr] = max(abs(corr_out));
tau_est_corr = lags(idx_peak_corr) / fs;   % Estimated delay
R_est_corr = tau_est_corr * c / 2;         % Estimated range

error_corr = R_est_corr - R_true;

fprintf('--- LFM Radar (Correlation) ---\n');
fprintf('True R = %.2f m | Estimated R = %.2f m | Error = %.2f m | Radar resolution = %.2f m\n', ...
    R_true, R_est_corr, error_corr, c/(2*lfm_bw));

%% ---------------- Plot Correlation Output ----------------
figure;
corr_time_us = lags/fs*1e6;
zoom_window = 200; 
zoom_range = max(1, idx_peak_corr-zoom_window):min(length(corr_out), idx_peak_corr+zoom_window);
plot(corr_time_us(zoom_range), abs(corr_out(zoom_range)),'LineWidth',1.5);
xlabel('Time (\mus)'); ylabel('Amplitude');
title('LFM Correlation Output (Zoomed Around Peak)');
grid on;
hold on;
plot(corr_time_us(idx_peak_corr), abs(corr_out(idx_peak_corr)),'ro','MarkerSize',8,'LineWidth',1.5);
text(corr_time_us(idx_peak_corr), abs(corr_out(idx_peak_corr))*1.05, ...
    sprintf('t=%.4fus', tau_est_corr*1e6));

%% ---------------- Frequency Domain of HPF Output (LFM) ----------------
figure;
[pxx_hpf,f_hpf] = pwelch(rx_radar_lfm, 4096, 2048, 4096, fs, 'centered');
plot(f_hpf/1e6, 10*log10(pxx_hpf),'LineWidth',1.5);
xlabel('Frequency (MHz)'); ylabel('Power/Frequency (dB/Hz)');
title('Spectrum of Received Radar Signal after High-Pass Filter (LFM)');
grid on;
xlim([lfm_start_freq/1e6-10, lfm_stop_freq/1e6+10]);

%% ---------------- Frequency Domain of LPF Output (OFDM) ----------------
figure;
[pxx_lpf,f_lpf] = pwelch(rxSig_ofdm_only, 4096, 2048, 4096, fs, 'centered');
plot(f_lpf/1e6, 10*log10(pxx_lpf),'LineWidth',1.5);
xlabel('Frequency (MHz)'); ylabel('Power/Frequency (dB/Hz)');
title('Spectrum of Received Radar Signal after Low-Pass Filter (OFDM)');
grid on;
xlim([0 ofdm_edge_freq/1e6+10]);
