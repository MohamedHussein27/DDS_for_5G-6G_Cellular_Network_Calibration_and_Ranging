%% OFDM + Chirp + Multipath + Noise (Robust with Pilot-Based Equalization)
clear; close all; clc;

%% -------------------------
% System Parameters
c = 3e8;                 % speed of light (m/s)
Fs = 491.5e6;            % sampling freq (Hz)
Ts = 1/Fs;

% OFDM
Nfft = 4096;
Ncp  = 128;
Nsym = 100;
M    = 256;

% --- MODIFIED OFDM PARAMETERS ---
subcarrier_spacing = 120e3;  
num_subcarriers_per_RB = 12; 
num_RB = 264;
numDataSubc = num_RB * num_subcarriers_per_RB;
% --- END MODIFIED PARAMETERS ---

% Pilot parameters
pilot_spacing = 10;      

% Chirp / radar
chirp_bw = 80e6;
f0_chirp = 0;             % Start frequency 0 Hz

% Noise / channel
SNR_combined_dB = 30;    
channel_gain = 1.0;

% Main target
target_distance = 302;   % meters
roundtrip_time = 2*target_distance/c;

%% Derived indices
half = floor(numDataSubc/2);
center = Nfft/2 + 1;
subc_idx = (center-half):(center-half+numDataSubc-1);
subc_idx = round(subc_idx);

% Pilot indices
pilot_idx_data = 1:pilot_spacing:numDataSubc;
pilot_idx_global = subc_idx(pilot_idx_data);
num_pilots = length(pilot_idx_data);

% Pilot sequence
pilot_seq = exp(1j*pi*(0:num_pilots-1).^2/num_pilots).';

numFrames = Nsym + 1;
frame_len = Nfft + Ncp;
L_total = frame_len * numFrames;
sample_delay = round(roundtrip_time*Fs);
if sample_delay >= L_total
    error('Echo delay too long for capture length');
end
fprintf('Capture length = %d samples, sample_delay = %d (target = %.2f m)\n', L_total, sample_delay, target_distance);

%% -------------------------
% 1) OFDM transmit with embedded pilots
bitsPerSym = log2(M);
numQAM = numDataSubc * Nsym;
tx_bits = randi([0 1], bitsPerSym * numQAM, 1);
tx_idx = bi2de(reshape(tx_bits,bitsPerSym,[]).','left-msb');
tx_qam = qammod(tx_idx, M, 'UnitAveragePower', true);
tx_qam = reshape(tx_qam, numDataSubc, Nsym);

FD = zeros(Nfft, numFrames);

% Preamble symbol
pilot_vec_full = exp(1j*2*pi*(0:numDataSubc-1)/numDataSubc).';
FD(subc_idx,1) = pilot_vec_full;

% Data symbols with embedded pilots
for k = 1:Nsym
    symbol_data = tx_qam(:,k);
    symbol_data(pilot_idx_data) = pilot_seq;
    FD(subc_idx,k+1) = symbol_data;
end

TD_frames = zeros(frame_len, numFrames);
for k = 1:numFrames
    xt = ifft(ifftshift(FD(:,k)));
    xt_cp = [xt(end-Ncp+1:end); xt];
    TD_frames(:,k) = xt_cp;
end
tx_ofdm = TD_frames(:);

%% -------------------------
% 2) Chirp signal
t = (0:L_total-1).' / Fs;
T_sim = t(end); 
% Linear chirp from 0 to chirp_bw
tx_chirp = chirp(t, f0_chirp, T_sim, chirp_bw, 'linear');

% Plot generated chirp (real part)
figure; plot(t, real(tx_chirp)); title('Generated Chirp Signal (Real Part)'); xlabel('Time (s)'); ylabel('Amplitude'); grid on;

% Scale chirp relative to OFDM
target_chirp_rel_db = 10;
P_ofdm = mean(abs(tx_ofdm).^2);
P_chirp = mean(abs(tx_chirp).^2);
tx_chirp = tx_chirp * sqrt(P_ofdm*10^(target_chirp_rel_db/10)/(P_chirp+eps));

%% -------------------------
% Chirp zoom & PSD
Nfft_plot = 4096;         
zoom_time_us = 90e-6;      
Nzoom = min(length(t), round(zoom_time_us * Fs));

figure;
subplot(3,1,1);
plot(t(1:Nzoom)*1e6, real(tx_chirp(1:Nzoom)), 'b');
xlabel('Time (\mus)'); ylabel('Amplitude');
title(sprintf('Chirp Time-Domain Zoom (first %.1f \\mus)', zoom_time_us*1e6));
grid on;

% Instantaneous frequency
analytic_chirp = hilbert(tx_chirp);
inst_phase = unwrap(angle(analytic_chirp));
inst_freq = (diff(inst_phase)/(2*pi)) * Fs; 
subplot(3,1,2);
plot(t(1:end-1)*1e6, inst_freq/1e6, 'r');
xlabel('Time (\mus)'); ylabel('Instantaneous Frequency (MHz)');
title('Instantaneous Frequency of Chirp (via Hilbert)');
grid on;
xlim([0 zoom_time_us*1e6]);

% PSD using Welch
pwelch_win = 2048;
pwelch_ovlp = round(0.5*pwelch_win);
pwelch_nfft = 4096;
[psd_chirp, f_psd] = pwelch(tx_chirp, pwelch_win, pwelch_ovlp, pwelch_nfft, Fs, 'centered');
subplot(3,1,3);
plot(f_psd/1e6, 10*log10(psd_chirp + eps), 'k', 'LineWidth', 1.2);
xlabel('Frequency (MHz)'); ylabel('PSD (dB/Hz)');
title('Chirp PSD (Welch method)');
grid on;
xlim([-Fs/2/1e6 Fs/2/1e6]);

%% -------------------------
% 3) Composite TX
tx_direct = tx_ofdm + tx_chirp;

%% -------------------------
% 4) Channel: Target Echo + AWGN
echo = zeros(L_total,1);
echo_atten = 0.1; 
len_main = L_total - sample_delay;
echo(sample_delay+1:end) = echo_atten * tx_chirp(1:len_main);

rx_before_noise = channel_gain * tx_direct + echo;
P_rx = mean(abs(rx_before_noise).^2);
SNR_combined_lin = 10^(SNR_combined_dB/10);
noise_pow_rx = P_rx / SNR_combined_lin;
noise_rx = sqrt(noise_pow_rx/2)*(randn(L_total,1)+1j*randn(L_total,1));
rx = rx_before_noise + noise_rx;

figure; plot(t, real(rx)); title('Received Signal (Direct + Echo + AWGN)'); xlabel('Time (s)'); ylabel('Amplitude'); grid on;

%% -------------------------
% 5) OFDM/Chirp Separation
rx_for_ofdm = rx - channel_gain*tx_chirp;
rx_residual = rx - channel_gain*tx_direct;

%% -------------------------
% 6) ROBUST OFDM RX (Pilot-Based Equalization)
rx_frames = reshape(rx_for_ofdm, frame_len, numFrames);
RX_FD = zeros(Nfft, numFrames);
for k = 1:numFrames
    rx_nocp = rx_frames(Ncp+1:end,k);
    RX_FD(:,k) = fftshift(fft(rx_nocp,Nfft));
end

rx_uneq = RX_FD(subc_idx, 2:end);
rx_uneq_vec = rx_uneq(:);

% Preamble EVM
Y_preamble_rx = RX_FD(subc_idx, 1);
X_preamble_tx = pilot_vec_full;
H_est_initial_fs = Y_preamble_rx ./ (X_preamble_tx + eps);
H_imperfect = mean(H_est_initial_fs);
X_eq_preamble_imperfect = Y_preamble_rx ./ (H_imperfect + eps);
err_preamble = X_eq_preamble_imperfect - X_preamble_tx;
evm_preamble = sqrt(mean(abs(err_preamble).^2)/mean(abs(X_preamble_tx).^2));
fprintf('Preamble EVM = %.4f (%.2f%%)\n', evm_preamble, evm_preamble*100);

% Enhanced equalization
Ypilot = RX_FD(:,1);
H_est_initial = ones(Nfft,1);
H_est_initial(subc_idx) = Ypilot(subc_idx)./(pilot_vec_full+eps);

rx_eq = zeros(numDataSubc, Nsym);
H_est_per_symbol = zeros(Nfft, Nsym);

for k = 1:Nsym
    Yk = RX_FD(:,k+1);
    Y_pilots = Yk(pilot_idx_global);
    H_pilots = Y_pilots ./ (pilot_seq + eps);

    if k == 1
        H_est_symbol = H_est_initial;
        H_est_symbol(pilot_idx_global) = H_pilots;
        H_est_interp = interp1(pilot_idx_global, H_pilots, subc_idx, 'linear', 'extrap');
        H_est_symbol(subc_idx) = H_est_interp.';
    else
        H_prev = H_est_per_symbol(:,k-1);
        H_est_symbol = H_prev;
        H_est_symbol(pilot_idx_global) = 0.7*H_prev(pilot_idx_global) + 0.3*H_pilots;
        H_est_interp = interp1(pilot_idx_global, H_pilots, subc_idx, 'linear', 'extrap');
        alpha = 0.8;
        H_est_symbol(subc_idx) = alpha*H_prev(subc_idx) + (1-alpha)*H_est_interp.';
    end

    H_est_per_symbol(:,k) = H_est_symbol;

    Y_data = Yk(subc_idx);
    X_est = Y_data ./ (H_est_symbol(subc_idx) + eps);
    X_est(pilot_idx_data) = pilot_seq;
    rx_eq(:,k) = X_est;
end

rx_eq_vec = rx_eq(:);
tx_sym_vec = tx_qam(:);

data_positions = true(numDataSubc,1);
data_positions(pilot_idx_data) = false;
rx_data_only = [];
tx_data_only = [];
for k=1:Nsym
    rx_data_only = [rx_data_only; rx_eq(data_positions,k)];
    tx_data_only = [tx_data_only; tx_qam(data_positions,k)];
end

%% -------------------------
% 7) Enhanced EVM
evm_per_symbol = zeros(Nsym,1);
for k=1:Nsym
    err_sym = rx_eq(data_positions,k) - tx_qam(data_positions,k);
    evm_per_symbol(k) = sqrt(mean(abs(err_sym).^2)/mean(abs(tx_qam(data_positions,k)).^2));
end
fprintf('EVM per symbol: [%s]%%\n', sprintf('%.2f ', evm_per_symbol*100));

%% -------------------------
% 8) Chirp echo detection (Matched Filter) vs time
[xc_res, lags_res] = xcorr(rx_residual, tx_chirp);
time_axis = lags_res / Fs;  % seconds

figure;
plot(time_axis, abs(xc_res), 'b', 'LineWidth', 1.2);
xlabel('Time (s)');
ylabel('Correlation Magnitude');
title('Matched Filter Output of Chirp (Residual) vs Time');
grid on;

% Peak detection
[~, idx_global] = max(abs(xc_res));
lag_global = lags_res(idx_global);

win = 50;
idx_window = find(lags_res >= sample_delay-win & lags_res <= sample_delay+win);
[~, rel_idx] = max(abs(xc_res(idx_window)));
idx_near = idx_window(rel_idx);
lag_near = lags_res(idx_near);


dist_global = lag_global*Ts*c/2;
dist_near   = lag_near*Ts*c/2;
fprintf('True distance = %.2f m (sample_delay=%d)\n', target_distance, sample_delay);
fprintf('Estimated distance (global peak) = %.4f m (lag=%d)\n', dist_global, lag_global);
fprintf('Estimated distance (near expected) = %.4f m (lag=%d)\n', dist_near, lag_near);

%% -------------------------
% 9) Enhanced plots (channel evolution, EVM, constellation)
figure;
subplot(2,1,1);
plot(subc_idx, 20*log10(abs(H_est_initial(subc_idx))), 'b-', 'LineWidth', 2);
hold on;
for k = 1:Nsym
    plot(subc_idx, 20*log10(abs(H_est_per_symbol(subc_idx,k))), 'Color', [0.5 0.5 0.5 0.3]);
end
plot(pilot_idx_global, 20*log10(abs(H_est_per_symbol(pilot_idx_global,1))), 'ro', 'MarkerSize', 4);
xlabel('Subcarrier Index'); ylabel('|H(f)| (dB)');
title('Channel Estimate Evolution (Pilot positions marked)');
legend('Initial Estimate','Per-symbol estimates','Pilot positions');
grid on;

subplot(2,1,2);
plot(1:Nsym, evm_per_symbol*100, 'bo-', 'LineWidth',2);
xlabel('Symbol Index'); ylabel('EVM (%)');
title('EVM vs Symbol Index with Enhanced Equalization');
grid on;

figure;
subplot(1,4,1);
scatter(real(tx_sym_vec(1:min(500,length(tx_sym_vec)))), imag(tx_sym_vec(1:min(500,length(tx_sym_vec)))), 8, 'filled');
title('1. Transmitted (256-QAM)'); xlabel('I'); ylabel('Q'); grid on; axis equal; xlim([-1.5 1.5]); ylim([-1.5 1.5]);
subplot(1,4,2);
scatter(real(rx_uneq_vec(1:min(500,length(rx_uneq_vec)))), imag(rx_uneq_vec(1:min(500,length(rx_uneq_vec)))), 8, 'filled');
title('2. Received (Before EQ)'); xlabel('I'); ylabel('Q'); grid on; axis equal; xlim([-1.5 1.5]); ylim([-1.5 1.5]);
subplot(1,4,3);
scatter(real(rx_eq_vec(1:min(500,length(rx_eq_vec)))), imag(rx_eq_vec(1:min(500,length(rx_eq_vec)))), 8, 'filled');
title('3. Received (After EQ)'); xlabel('I'); ylabel('Q'); grid on; axis equal; xlim([-1.5 1.5]); ylim([-1.5 1.5]);
subplot(1,4,4);
scatter(real(rx_data_only(1:min(500,length(rx_data_only)))), imag(rx_data_only(1:min(500,length(rx_data_only)))), 8, 'filled');
title('4. Data Only (Pilots Removed)'); xlabel('I'); ylabel('Q'); grid on; axis equal; xlim([-1.5 1.5]); ylim([-1.5 1.5]);

%% -------------------------
% Original transmitted chirp
figure;
subplot(2,1,1);
plot(t, abs(tx_chirp), 'b'); xlabel('Time (s)'); ylabel('Magnitude'); title('Original Transmitted Chirp Signal'); grid on;
subplot(2,1,2);
plot(t, unwrap(angle(tx_chirp)), 'r'); xlabel('Time (s)'); ylabel('Phase (rad)'); title('Unwrapped Phase of Original Transmitted Chirp'); grid on;

%% -------------------------
% 10) Summary
fprintf('\n=== System Summary ===\n');
fprintf('Subcarrier spacing: %.2f Hz\n', subcarrier_spacing);
fprintf('Number of resource blocks: %d\n', num_RB);
fprintf('Pilot spacing: every %d subcarriers\n', pilot_spacing);
fprintf('Number of pilots per symbol: %d\n', num_pilots);
fprintf('Total data subcarriers: %d\n', numDataSubc);
fprintf('Data rate improvement with pilot-based tracking: Enabled\n');
fprintf('Chirp bandwidth: 0 to %.2f MHz\n', chirp_bw/1e6);
