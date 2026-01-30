% LFM radar simulation with matched filtering (AWGN only) - WITH PSD
clear; close all; clc;

%% ----------------- Parameters -----------------
c = 3e8;                 % speed of light (m/s)
A = 1;                   % amplitude of transmitted chirp
f0 = 0;                % start frequency (Hz) (baseband-ish example)
B = 80e6;                 % bandwidth (Hz)
T = 100e-6;              % pulse duration (seconds)
K = B / T;               % chirp rate (Hz/s) -> K = B/T
theta = 0;               % initial phase (rad)
Fs = 5e8;                % sampling frequency (Hz) (must be >> f0+B)
SNR_dB = 10;             % desired SNR at the receiver (dB) for the echo
R_true = 300;            % true target range in meters
alpha = 0.5;             % target amplitude attenuation (reflection coef)

%% ----------------- Time vectors -----------------
t = 0:1/Fs:(T - 1/Fs);   % time vector for one pulse (length N)
N = length(t);

%% ----------------- Transmitted signal -----------------
% Note: phase = 2*pi*(f0*t + K/2 * t.^2) + theta
s = A*cos(2*pi*(f0*t + 0.5*K*t.^2) + theta);    % real chirp signal
%s = chirp(t, f0, T, f0+B, 'linear');           %built-in 
% Create complex analytic version for frequency analysis
s_complex = A*exp(1j*(2*pi*(f0.*t + 0.5*K.*t.^2) + theta));

%% ----------------- Frequency analysis -----------------
% Frequency vector for FFT
f = (-N/2:N/2-1)*(Fs/N);

% FFT of transmitted signal
S = fftshift(fft(s_complex));

%% ----------------- Create received echo (delayed + attenuated) ---------------
tau = 2*R_true / c;                     % round-trip delay (s)
delay_samples = round(tau * Fs);        % integer sample delay
fprintf('True round-trip delay = %.3e s -> delay samples = %d\n', tau, delay_samples);

% Received vector length: we will create a received record long enough to contain delay + pulse
rec_len = delay_samples + N + 100;      % pad a little after pulse
received = zeros(1, rec_len);
received_complex = zeros(1, rec_len);

% place the delayed, attenuated transmitted pulse in received
start_idx = delay_samples + 1;          % MATLAB indexing
received(start_idx : start_idx + N - 1) = alpha * s;
received_complex(start_idx : start_idx + N - 1) = alpha * s_complex;

%% ----------------- Add AWGN Manually (Correct Method) -----------------
% Compute echo signal power (real)
Psig = mean((alpha*s).^2);

% Compute noise sigma for the CHOSEN SNR
SNR_lin = 10^(SNR_dB/10);
noise_variance = Psig / SNR_lin;
sigma = sqrt(noise_variance);

% Add white Gaussian noise (real and complex)
noise_real = 100 * sigma * randn(size(received));
received_noisy = received + noise_real;

noise_complex = sigma/sqrt(2) * (randn(size(received)) + 1j*randn(size(received)));
received_complex_noisy = received_complex + noise_complex;

% Compute ACTUAL SNR inside echo window
echo_window = received_noisy(start_idx:start_idx+N-1);
noise_window = echo_window - alpha*s;

actual_SNR_dB = 10*log10( Psig / mean(noise_window.^2) );
fprintf('Requested SNR = %.1f dB, Actual SNR = %.1f dB\n', SNR_dB, actual_SNR_dB);


%% ----------------- Matched filter (correlation) --------------------------
% matched filter for real signal s is h = time-reversed s (no conjugate needed for real)
h = fliplr(s);                  % matched filter impulse response

% convolution = matched filtering
z = conv(received_noisy, h);    % length = rec_len + N - 1

% time axis for z
t_z = (0:length(z)-1) / Fs;

% Find peak location (global max)
[~, idx_max] = max(abs(z));
% convert to equivalent delay in samples: idx_max corresponds to (delay_samples + N)
est_delay_samples = idx_max - N;    % as discussed in analysis
est_tau = est_delay_samples / Fs;
R_est = c * est_tau / 2;

fprintf('Estimated delay samples = %d -> est_tau = %.3e s -> R_est = %.3f m\n',...
        est_delay_samples, est_tau, R_est);

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
title('Received Signal (delayed & noisy) - Time Domain');
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
title('Received Signal (Showing Actual Delay)');
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
xlabel('Delay (µs) relative to pulse)'); ylabel('|z|'); title('Matched filter peak (zoom)');
grid on;
xline(est_delay_samples / Fs * 1e6, '--r', sprintf('Est delay = %.2f \\mus', est_tau*1e6));

%% ----------------- Ranging estimation explanation & performance ----------
% Range resolution (theoretical) = c/(2*B)
range_resolution = c / (2*B);
fprintf('Theoretical range resolution (c / (2B)) = %.3f m\n', range_resolution);

% Ranging estimation: here we compute simple error for this single run
range_error = R_est - R_true;
fprintf('Range estimate error (R_est - R_true) = %.3f m\n', range_error);

%% ----------------- Verification: Show Actual Delay -----------------
fprintf('\n=== Delay Verification ===\n');
fprintf('Theoretical time delay for %.0f m target: %.6f µs\n', R_true, tau*1e6);
fprintf('Actual echo starts at: %.6f µs\n', delay_samples/Fs*1e6);
fprintf('Time difference: %.6f µs\n', abs(tau*1e6 - delay_samples/Fs*1e6));