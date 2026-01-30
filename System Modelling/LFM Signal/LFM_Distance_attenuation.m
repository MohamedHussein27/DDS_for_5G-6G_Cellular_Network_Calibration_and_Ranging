% LFM radar simulation with proper range-dependent attenuation (R^4 law) - COMPLETE PLOTS
clear; close all; clc;

%% ----------------- Parameters -----------------
c = 3e8;                 % speed of light (m/s)
A = 1;                   % amplitude of transmitted chirp
f0 = 0;                  % start frequency (Hz) (baseband-ish example)
B = 200e6;                % bandwidth (Hz)
T = 50e-6;               % pulse duration (seconds)
K = B / T;               % chirp rate (Hz/s) -> K = B/T
theta = 0;               % initial phase (rad)
Fs = 5e8;                % sampling frequency (Hz) (must be >> f0+B)
R_true = 300;            % true target range in meters

% Radar system parameters
Pt = 100000;               % Transmitted power (W)
Gt = 50;                 % Transmitter antenna gain (dBi)
Gr = 50;                 % Receiver antenna gain (dBi)
lambda = 0.03;           % Wavelength at 10 GHz (m)
sigma_rcs = 1;           % Target Radar Cross Section (m²)
L = 2;                   % System losses (linear factor > 1)

% Reference SNR at 1 meter (for calibration)
SNR_ref_dB = 35;         % SNR at reference range of 1 meter (dB)

fprintf('Radar System Parameters:\n');
fprintf('Transmitted Power: %.0f W\n', Pt);
fprintf('Antenna Gains: %.0f dBi\n', Gt);
fprintf('Wavelength: %.3f m (%.1f GHz)\n', lambda, c/lambda/1e9);
fprintf('Target RCS: %.1f m²\n', sigma_rcs);
fprintf('System Losses: %.1f (%.1f dB)\n', L, 10*log10(L));

%% ----------------- Time vectors -----------------
t = 0:1/Fs:(T - 1/Fs);   % time vector for one pulse (length N)
N = length(t);

%% ----------------- Transmitted signal -----------------
s = A*cos(2*pi*(f0*t + 0.5*K*t.^2) + theta);    % real chirp signal
s_complex = A*exp(1j*(2*pi*(f0.*t + 0.5*K.*t.^2) + theta));

%% ----------------- Frequency analysis -----------------
f = (-N/2:N/2-1)*(Fs/N);
S = fftshift(fft(s_complex));

%% ----------------- Calculate Range-Dependent Attenuation -----------------
% Convert gains from dB to linear
Gt_lin = 10^(Gt/10);
Gr_lin = 10^(Gr/10);

% Calculate received power using Radar Range Equation
Pr = (Pt * Gt_lin * Gr_lin * lambda^2 * sigma_rcs) / ((4*pi)^3 * R_true^4 * L);

% Calculate attenuation factor (alpha) from radar equation
% Received power = alpha^2 * Transmitted power
alpha_range = sqrt(Pr / Pt);

% Convert reference SNR to actual SNR at target range
SNR_ref_lin = 10^(SNR_ref_dB/10);
SNR_actual_lin = SNR_ref_lin * (1/R_true^4);  % SNR scales with 1/R^4
SNR_actual_dB = 10*log10(SNR_actual_lin);

fprintf('\n=== Range-Dependent Parameters ===\n');
fprintf('Target Range: %.0f m\n', R_true);
fprintf('Received Power: %.2e W\n', Pr);
fprintf('Power Attenuation Factor (alpha): %.6f\n', alpha_range);
fprintf('Reference SNR at 1m: %.1f dB\n', SNR_ref_dB);
fprintf('Actual SNR at %.0fm: %.1f dB\n', R_true, SNR_actual_dB);
fprintf('Attenuation due to range: %.1f dB\n', SNR_ref_dB - SNR_actual_dB);

%% ----------------- Create received echo (with proper attenuation) ---------------
tau = 2*R_true / c;                     % round-trip delay (s)
delay_samples = round(tau * Fs);        % integer sample delay
fprintf('\nTrue round-trip delay = %.3e s -> delay samples = %d\n', tau, delay_samples);

% Received vector length
rec_len = delay_samples + N + 100;
received = zeros(1, rec_len);
received_complex = zeros(1, rec_len);

% Place the delayed, range-attenuated transmitted pulse
start_idx = delay_samples + 1;
received(start_idx : start_idx + N - 1) = alpha_range * s;
received_complex(start_idx : start_idx + N - 1) = alpha_range * s_complex;

%% ----------------- Add AWGN with Range-Dependent SNR -----------------
% Compute echo signal power (real)
Psig_actual = mean((alpha_range * s).^2);

% Compute noise sigma for the ACTUAL range-dependent SNR
noise_variance_actual = Psig_actual / SNR_actual_lin;
sigma_actual = sqrt(noise_variance_actual);

% Add white Gaussian noise
noise_real = sigma_actual * randn(size(received));
received_noisy = received + noise_real;

noise_complex = sigma_actual/sqrt(2) * (randn(size(received)) + 1j*randn(size(received)));
received_complex_noisy = received_complex + noise_complex;

% Compute ACTUAL SNR inside echo window
echo_window = received_noisy(start_idx:start_idx+N-1);
noise_only_window = received_noisy(1:delay_samples);  % Noise before echo
actual_SNR_measured_dB = 10*log10(mean(echo_window.^2) / mean(noise_only_window.^2));

fprintf('\n=== SNR Verification ===\n');
fprintf('Theoretical SNR at %.0f m: %.1f dB\n', R_true, SNR_actual_dB);
fprintf('Measured SNR in echo window: %.1f dB\n', actual_SNR_measured_dB);

%% ----------------- Matched filter (correlation) --------------------------
h = fliplr(s);                  % matched filter impulse response
z = conv(received_noisy, h);    % length = rec_len + N - 1
t_z = (0:length(z)-1) / Fs;

% Find peak location (global max)
[~, idx_max] = max(abs(z));
est_delay_samples = idx_max - N;
est_tau = est_delay_samples / Fs;
R_est = c * est_tau / 2;

fprintf('\n=== Range Estimation ===\n');
fprintf('Estimated delay samples = %d\n', est_delay_samples);
fprintf('Estimated time delay = %.3e s\n', est_tau);
fprintf('Estimated range = %.3f m\n', R_est);
fprintf('Range error = %.3f m\n', R_est - R_true);

%% ----------------- PSD Analysis -----------------
% Calculate PSD using pwelch method for better spectral estimation
[pxx_trans, f_trans] = pwelch(s_complex, hanning(1024), 512, 1024, Fs, 'centered', 'psd');
[pxx_rec, f_rec] = pwelch(received_complex_noisy, hanning(1024), 512, 1024, Fs, 'centered', 'psd');
[pxx_noise, f_noise] = pwelch(noise_complex, hanning(1024), 512, 1024, Fs, 'centered', 'psd');

% Convert to dB scale for better visualization
pxx_trans_db = 10*log10(pxx_trans);
pxx_rec_db = 10*log10(pxx_rec);
pxx_noise_db = 10*log10(pxx_noise);

%% ----------------- PLOTS: Time and Frequency Domain --------------------------
figure('Name','LFM Signal Analysis - Time & Frequency','Position',[100 100 1200 900]);

% Plot 1: Transmitted signal in time domain
subplot(3,2,1);
plot(t*1e6, s);
xlabel('Time (µs)'); ylabel('Amplitude'); 
title('Transmitted LFM Chirp - Time Domain');
grid on;

% Plot 2: Transmitted signal in frequency domain
subplot(3,2,2);
plot(f/1e3, abs(S));
xlabel('Frequency (kHz)'); ylabel('Magnitude');
title('Transmitted LFM Chirp - Frequency Domain');
xlim([(f0-B/2)/1e3 (f0+B*1.5)/1e3]); % Zoom around chirp frequencies
grid on;

% Plot 3: Instantaneous frequency of transmitted signal
subplot(3,2,3);
% Theoretical instantaneous frequency: f_inst = f0 + K*t
f_inst = f0 + K*t;
plot(t*1e6, f_inst/1e3);
xlabel('Time (µs)'); ylabel('Frequency (kHz)');
title('Instantaneous Frequency of LFM Chirp');
grid on;

% Plot 4: Received signal in time domain - CORRECTED TIME AXIS
subplot(3,2,4);
% CORRECTED: Create proper time axis that accounts for the delay
t_received = (0:length(received_noisy)-1) / Fs * 1e6;  % Time in microseconds
plot(t_received, received_noisy);
xlabel('Time (µs)'); ylabel('Amplitude'); 
title(sprintf('Received Signal (Range: %.0f m, SNR: %.1f dB)', R_true, SNR_actual_dB));
xlim([0 (delay_samples + N + 50)/Fs*1e6]);
grid on;
% Add vertical line to show where the echo starts
hold on;
xline(delay_samples/Fs*1e6, '--r', 'LineWidth', 2, 'DisplayName', 'Echo start');
xline((delay_samples + N)/Fs*1e6, '--r', 'LineWidth', 2, 'DisplayName', 'Echo end');
legend('show', 'Location', 'best');
hold off;

% Plot 5: Received signal in frequency domain
subplot(3,2,5);
% FFT of received signal (complex version for better frequency analysis)
S_received = fftshift(fft(received_complex_noisy));
f_received = (-length(S_received)/2:length(S_received)/2-1)*(Fs/length(S_received));
plot(f_received/1e3, abs(S_received));
xlabel('Frequency (kHz)'); ylabel('Magnitude');
title('Received Signal - Frequency Domain');
xlim([(f0-B/2)/1e3 (f0+B*1.5)/1e3]);
grid on;

% Plot 6: Matched filter output
subplot(3,2,6);
plot(t_z*1e6, abs(z));
xlabel('Time (µs)'); ylabel('|z(t)|'); 
title('Matched Filter Output |z(t)| (correlation)');
grid on;
% highlight peak
hold on;
plot(est_delay_samples / Fs * 1e6, abs(z(idx_max)), 'ro','MarkerSize',8,'LineWidth',1.5);
text(est_delay_samples / Fs * 1e6, abs(z(idx_max))*0.9, sprintf(' Est peak -> R=%.2f m', R_est));

%% ----------------- NEW: Time Domain Comparison Plot -----------------
figure('Name','Time Domain Signal Comparison','Position',[400 200 1000 600]);

% Create proper time axes for comparison
t_transmitted = t * 1e6;  % Transmitted signal time (µs)
t_received_full = (0:length(received_noisy)-1) / Fs * 1e6;  % Received signal time (µs)

subplot(2,1,1);
plot(t_transmitted, s, 'b', 'LineWidth', 1.5);
xlabel('Time (µs)'); ylabel('Amplitude');
title('Transmitted LFM Signal');
grid on;
xlim([0 max(t_transmitted)]);

subplot(2,1,2);
plot(t_received_full, received_noisy, 'r', 'LineWidth', 1);
xlabel('Time (µs)'); ylabel('Amplitude');
title(sprintf('Received Signal (Range: %.0f m, SNR: %.1f dB)', R_true, SNR_actual_dB));
grid on;
% Mark the echo region
hold on;
echo_start = delay_samples/Fs*1e6;
echo_end = (delay_samples + N)/Fs*1e6;
xline(echo_start, '--g', 'LineWidth', 2, 'DisplayName', 'Echo start');
xline(echo_end, '--g', 'LineWidth', 2, 'DisplayName', 'Echo end');
% Shade the echo region
xFill = [echo_start, echo_end, echo_end, echo_start];
yFill = [min(ylim), min(ylim), max(ylim), max(ylim)];
fill(xFill, yFill, 'g', 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'DisplayName', 'Echo region');
legend('show', 'Location', 'best');
hold off;

%% ----------------- NEW: PSD Analysis Plots -----------------
figure('Name','Power Spectral Density Analysis','Position',[300 100 1200 800]);

% Plot 1: PSD comparison (linear scale)
subplot(2,2,1);
plot(f_trans/1e6, pxx_trans, 'b', 'LineWidth', 1.5, 'DisplayName', 'Transmitted Chirp');
hold on;
plot(f_rec/1e6, pxx_rec, 'r', 'LineWidth', 1, 'DisplayName', 'Received Signal');
plot(f_noise/1e6, pxx_noise, 'g', 'LineWidth', 0.5, 'DisplayName', 'Noise Only');
xlabel('Frequency (MHz)'); ylabel('Power Spectral Density (W/Hz)');
title('PSD Comparison - Linear Scale');
legend;
grid on;
xlim([-B/1e6*1.2, B/1e6*1.2]);

% Plot 2: PSD comparison (dB scale)
subplot(2,2,2);
plot(f_trans/1e6, pxx_trans_db, 'b', 'LineWidth', 1.5, 'DisplayName', 'Transmitted Chirp');
hold on;
plot(f_rec/1e6, pxx_rec_db, 'r', 'LineWidth', 1, 'DisplayName', 'Received Signal');
plot(f_noise/1e6, pxx_noise_db, 'g', 'LineWidth', 0.5, 'DisplayName', 'Noise Only');
xlabel('Frequency (MHz)'); ylabel('Power Spectral Density (dB/Hz)');
title('PSD Comparison - dB Scale');
legend;
grid on;
xlim([-B/1e6*1.2, B/1e6*1.2]);

% Plot 3: Zoomed PSD around chirp bandwidth
subplot(2,2,3);
plot(f_trans/1e6, pxx_trans_db, 'b', 'LineWidth', 2);
xlabel('Frequency (MHz)'); ylabel('PSD (dB/Hz)');
title('Transmitted LFM Chirp PSD (Zoomed)');
grid on;
xlim([-B/1e6, B/1e6]);
% Add bandwidth annotation
yl = ylim;
hold on;
plot([-B/2e6, B/2e6], [yl(1)+0.9*diff(yl), yl(1)+0.9*diff(yl)], 'r--', 'LineWidth', 2);
text(0, yl(1)+0.85*diff(yl), sprintf('Bandwidth = %.1f MHz', B/1e6), ...
     'HorizontalAlignment', 'center', 'BackgroundColor', 'white');

% Plot 4: PSD of noise component
subplot(2,2,4);
plot(f_noise/1e6, pxx_noise_db, 'k', 'LineWidth', 1);
xlabel('Frequency (MHz)'); ylabel('PSD (dB/Hz)');
title('Noise PSD (White Gaussian Noise)');
grid on;
xlim([-Fs/2e6*0.1, Fs/2e6*0.1]); % Zoom to see flat characteristic
% Calculate and display noise floor
noise_floor_dB = mean(pxx_noise_db);
text(0, noise_floor_dB+3, sprintf('Noise Floor = %.1f dB/Hz', noise_floor_dB), ...
     'HorizontalAlignment', 'center', 'BackgroundColor', 'white');

%% ----------------- Additional Detailed Frequency Analysis -----------------
figure('Name','Detailed Frequency Analysis','Position',[200 200 1000 700]);

% Spectrogram of transmitted signal
subplot(2,2,1);
window = 256;
noverlap = 200;
nfft = 1024;
spectrogram(s_complex, window, noverlap, nfft, Fs, 'yaxis');
title('Spectrogram of Transmitted LFM Chirp');
colorbar;

% Frequency content comparison
subplot(2,2,2);
hold on;
plot(f/1e3, abs(S)/max(abs(S)), 'b', 'LineWidth', 1.5, 'DisplayName', 'Transmitted');
plot(f_received/1e3, abs(S_received)/max(abs(S_received)), 'r', 'LineWidth', 1, 'DisplayName', 'Received');
xlabel('Frequency (kHz)'); ylabel('Normalized Magnitude');
title('Normalized Frequency Spectrum Comparison');
xlim([(f0-B/2)/1e3 (f0+B*1.5)/1e3]);
legend;
grid on;

% Phase progression
subplot(2,2,3);
phase_unwrapped = unwrap(angle(s_complex));
plot(t*1e6, phase_unwrapped);
xlabel('Time (µs)'); ylabel('Phase (rad)');
title('Unwrapped Phase of LFM Chirp');
grid on;

% NEW: PSD of matched filter output
subplot(2,2,4);
% Calculate PSD of matched filter output
[pxx_mf, f_mf] = pwelch(z, hanning(1024), 512, 1024, Fs, 'centered', 'psd');
plot(f_mf/1e6, 10*log10(pxx_mf));
xlabel('Frequency (MHz)'); ylabel('PSD (dB/Hz)');
title('PSD of Matched Filter Output');
grid on;

%% ----------------- Matched Filter Peak Zoom (original plot) -----------------
figure('Name','Matched Filter Peak Zoom','Position',[200 200 500 350]);
win = max(1, idx_max - 200) : min(length(z), idx_max + 200);
plot((win - N)/Fs*1e6, abs(z(win)));
xlabel('Delay (µs) relative to pulse)'); ylabel('|z|'); 
title('Matched filter peak (zoom)');
grid on;
xline(est_delay_samples / Fs * 1e6, '--r', sprintf('Est delay = %.2f \\mus', est_tau*1e6));

%% ----------------- NEW: Range-Dependent Analysis Plots -----------------
figure('Name','Range-Dependent Radar Performance','Position',[100 50 1400 900]);

% Multiple Range Analysis
ranges = [100, 300, 1000, 3000, 10000];  % Different target ranges
SNR_values = zeros(size(ranges));
alpha_values = zeros(size(ranges));

for i = 1:length(ranges)
    R = ranges(i);
    Pr_i = (Pt * Gt_lin * Gr_lin * lambda^2 * sigma_rcs) / ((4*pi)^3 * R^4 * L);
    alpha_i = sqrt(Pr_i / Pt);
    SNR_i_dB = SNR_ref_dB + 40*log10(1/R);  % 40 because 20*log10(alpha) and alpha ~ 1/R^2
    
    SNR_values(i) = SNR_i_dB;
    alpha_values(i) = alpha_i;
end

% Plot 1: Range vs Attenuation (R^4 law)
subplot(2,3,1);
semilogx(ranges, 20*log10(alpha_values), 'b-o', 'LineWidth', 2);
xlabel('Range (m)'); ylabel('Attenuation (dB)');
title('Range vs Signal Attenuation (R^4 Law)');
grid on;

% Plot 2: Range vs SNR
subplot(2,3,2);
semilogx(ranges, SNR_values, 'r-o', 'LineWidth', 2);
xlabel('Range (m)'); ylabel('SNR (dB)');
title('Range vs SNR');
grid on;

% Plot 3: Detection threshold analysis
subplot(2,3,3);
ranges_detection = logspace(2, 5, 100);  % 100 m to 100 km
SNR_detection = SNR_ref_dB + 40*log10(1./ranges_detection);

semilogx(ranges_detection, SNR_detection, 'b', 'LineWidth', 2);
hold on;
% Typical detection thresholds
xline(10^(SNR_ref_dB/40), '--r', 'LineWidth', 1, 'DisplayName', '0 dB SNR range');
yline(0, '--k', 'LineWidth', 1, 'DisplayName', 'Detection Threshold');
xlabel('Range (m)'); ylabel('SNR (dB)');
title('Detection Range Analysis');
legend; grid on;

% Plot 4: Signal power vs range
subplot(2,3,4);
Pr_detection = (Pt * Gt_lin * Gr_lin * lambda^2 * sigma_rcs) ./ ((4*pi)^3 * ranges_detection.^4 * L);
semilogx(ranges_detection, 10*log10(Pr_detection*1e3), 'r', 'LineWidth', 2);
xlabel('Range (m)'); ylabel('Received Power (dBm)');
title('Received Power vs Range');
grid on;

% Plot 5: Radar performance metrics
subplot(2,3,5);
range_resolution = c/(2*B);
PRI_max = 1/1000;  % Maximum PRI for 1 kHz PRF
R_unambiguous = c * PRI_max / 2;
SNR_min_dB = 0;  % Minimum detectable SNR
R_max = 10^((SNR_ref_dB - SNR_min_dB)/40);

bar_data = [range_resolution, R_unambiguous, R_max];
bar_categories = {'Range Resolution', 'Unambiguous Range', 'Max Detection Range'};
bar(bar_data);
set(gca, 'XTickLabel', bar_categories);
ylabel('Range (m)');
title('Radar Performance Metrics');
grid on;

% Add values on bars
for i = 1:3
    text(i, bar_data(i)*1.05, sprintf('%.0f m', bar_data(i)), ...
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
end

% Plot 6: Performance at different ranges
subplot(2,3,6);
test_ranges = [100, 300, 1000, 3000];
detection_probability = [0.99, 0.95, 0.7, 0.3];  % Example values
plot(test_ranges, detection_probability, 'g-o', 'LineWidth', 2);
xlabel('Range (m)'); ylabel('Detection Probability');
title('Detection Performance vs Range');
grid on;

%% ----------------- Ranging estimation explanation & performance ----------
% Range resolution (theoretical) = c/(2*B)
range_resolution = c / (2*B);
fprintf('\nTheoretical range resolution (c / (2B)) = %.3f m\n', range_resolution);

% Ranging estimation: here we compute simple error for this single run
range_error = R_est - R_true;
fprintf('Range estimate error (R_est - R_true) = %.3f m\n', range_error);

%% ----------------- Verification: Show Actual Delay -----------------
fprintf('\n=== Delay Verification ===\n');
fprintf('Theoretical time delay for %.0f m target: %.6f µs\n', R_true, tau*1e6);
fprintf('Actual echo starts at: %.6f µs\n', delay_samples/Fs*1e6);
fprintf('Time difference: %.6f µs\n', abs(tau*1e6 - delay_samples/Fs*1e6));

%% ----------------- Verification of R^4 Law -----------------
fprintf('\n=== Verification of R^4 Law ===\n');
fprintf('At 100m:  SNR = %.1f dB\n', SNR_ref_dB + 40*log10(1/100));
fprintf('At 200m:  SNR = %.1f dB (should be %.1f dB less than 100m)\n', ...
        SNR_ref_dB + 40*log10(1/200), 40*log10(200/100));
fprintf('At 400m:  SNR = %.1f dB (should be %.1f dB less than 100m)\n', ...
        SNR_ref_dB + 40*log10(1/400), 40*log10(400/100));
fprintf('At 800m:  SNR = %.1f dB (should be %.1f dB less than 100m)\n', ...
        SNR_ref_dB + 40*log10(1/800), 40*log10(800/100));

fprintf('\nEvery doubling of range reduces SNR by %.1f dB\n', 40*log10(2));

% Maximum detection range calculation
fprintf('\nMaximum Detection Range: %.0f m (SNR > %.0f dB)\n', R_max, SNR_min_dB);
fprintf('Unambiguous Range: %.0f m (PRF = %.0f Hz)\n', R_unambiguous, 1/PRI_max);