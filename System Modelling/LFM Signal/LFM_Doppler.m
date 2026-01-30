% LFM radar simulation with matched filtering (AWGN + multipath + Doppler)
% Author: ChatGPT (example)
clear; close all; clc;

%% ----------------- Parameters -----------------
c = 3e8;                 % speed of light (m/s)
A = 1;                   % amplitude of transmitted chirp
f0 = 0;                % start frequency (Hz) (baseband-ish example)
B = 80e6;                 % bandwidth (Hz)
T = 50e-6;              % pulse duration (seconds)
K = B / T;               % chirp rate (Hz/s) -> K = B/T
theta = 0;               % initial phase (rad)
Fs = 5e8;                % sampling frequency (Hz) (must be >> f0+B)
SNR_dB = 0;             % desired SNR at the receiver (dB) for the echo
R_true = 300;            % true target range in meters
alpha = 0.3;             % target amplitude attenuation (reflection coef)

%% ----------------- Multipath parameters (editable) -----------------
% num_mp is number of multipath components (additional to direct path).
num_mp = 5;

% Define extra-range (meters) for each multipath relative to the direct path range.
% Example: [30, 55] means one reflection travels an extra 30 m round-trip distance
% compared to the direct path (not total path). Adjust as needed.
mp_extra_range = [30, 25, 10, 20, 25];        % extra meters (positive values)
mp_alpha       = [0.4, 0.2, 0.1, 0.01, 0.03]; % multipath amplitudes (linear, <1 typically)
mp_phase       = [0, pi/6, pi/4, pi/6, pi/2]; % multipath phase offsets in radians

% Validate size
if length(mp_extra_range) < num_mp
    error('mp_extra_range length must be >= num_mp');
end
if length(mp_alpha) < num_mp
    error('mp_alpha length must be >= num_mp');
end
if length(mp_phase) < num_mp
    mp_phase = zeros(1, num_mp);
end

%% ----------------- Time vectors -----------------
t = 0:1/Fs:(T - 1/Fs);   % time vector for one pulse (length N)
N = length(t);

%% ----------------- Transmitted signal -----------------
% Use MATLAB chirp for the real waveform (same instantaneous phase used elsewhere)
s = chirp(t, f0, T, f0 + B, 'linear') * A;

% Create complex analytic version for frequency analysis and complex processing
s_complex = A * exp(1j*(2*pi*(f0.*t + 0.5*K.*t.^2) + theta));

%% ----------------- Pre-compute delays (main + multipath) -----------------
% main (direct) path
tau_main = 2 * R_true / c;                  % round-trip delay (s)
delay_samples_main = round(tau_main * Fs);

% multipath extra delays (convert extra ranges to extra round-trip delays)
% Here mp_extra_range is the extra round-trip distance difference in meters
mp_tau_extra = 2 * mp_extra_range / c;     % seconds (round-trip extra)
mp_delay_samples_extra = round(mp_tau_extra * Fs);

% compute the required record length to hold the latest multipath
all_delays_samples = delay_samples_main + [0, mp_delay_samples_extra];
max_delay_samples = max(all_delays_samples);

% allocate received buffer with safety padding
rec_len = max_delay_samples + N + 200;  % extra padding after last echo
received = zeros(1, rec_len);
received_complex = zeros(1, rec_len);

fprintf('Main delay (s) = %.3e  -> main delay samples = %d\n', tau_main, delay_samples_main);
for k = 1:num_mp
    fprintf('MP %d extra range = %.1f m -> extra delay (samples) = %d\n', ...
        k, mp_extra_range(k), mp_delay_samples_extra(k));
end
fprintf('rec_len = %d samples (Fs = %.3g Hz, duration %.3g s)\n', rec_len, Fs, rec_len/Fs);

%% ----------------- DOPPLER PARAMETERS -----------------
v_target = 30;                   % radial velocity (m/s), positive = approaching
f_c = f0 + B/2;                  % carrier frequency (central frequency of chirp)
f_d = 5 * v_target * f_c / c;    % Doppler shift (Hz)
fprintf('Doppler shift = %.2f Hz\n', f_d);

%% ----------------- Insert direct-path (delayed + attenuated) with Doppler -----------------
start_idx = delay_samples_main + 1;          % MATLAB indexing
end_idx   = start_idx + N - 1;

if end_idx <= rec_len
    % Apply Doppler to main path
    doppler_phase_main = 2*pi*f_d*t;
    s_main_complex = alpha * s_complex .* exp(1j*doppler_phase_main);
    s_main_real = real(s_main_complex);
    
    received(start_idx : end_idx) = received(start_idx : end_idx) + s_main_real;
    received_complex(start_idx : end_idx) = received_complex(start_idx : end_idx) + s_main_complex;
else
    error('Direct path exceeds buffer length. Increase rec_len or Fs.');
end

%% ----------------- Insert multipath echoes with Doppler -----------------
for k = 1:num_mp
    mp_start = start_idx + mp_delay_samples_extra(k);
    mp_end   = mp_start + N - 1;
    if mp_end <= rec_len
        % Apply Doppler to multipath (same Doppler shift for all paths from same target)
        doppler_phase_mp = 2*pi*f_d*t;
        s_mp_complex = mp_alpha(k) * exp(1j * mp_phase(k)) .* s_complex .* exp(1j*doppler_phase_mp);
        s_mp_real = real(s_mp_complex);
        
        % add real multipath copy (real-domain)
        received(mp_start:mp_end) = received(mp_start:mp_end) + s_mp_real;
        % add complex multipath (complex baseband)
        received_complex(mp_start:mp_end) = received_complex(mp_start:mp_end) + s_mp_complex;
    else
        warning('Multipath %d exceeds buffer length and was not added.', k);
    end
end

fprintf('Added direct path and %d multipath components into received buffer.\n', num_mp);

%% ----------------- Add AWGN using built-in function -----------------
% Using MATLAB's built-in awgn function (measured power)
received_noisy = awgn(received, SNR_dB, 'measured');
received_complex_noisy = awgn(received_complex, SNR_dB, 'measured');

% For logging: estimate actual noise
actual_noise = received_noisy - received;
Psig_recv = mean((alpha * s).^2);      % average power of attenuated chirp (real)
actual_sigma = std(actual_noise);
actual_SNR_dB = 10*log10(Psig_recv / (actual_sigma^2));

fprintf('Requested SNR = %.1f dB, Actual (measured) SNR = %.2f dB\n', SNR_dB, actual_SNR_dB);
fprintf('Received echo nominal power (approx) = %.3e ; noise sigma (measured) = %.3e\n',...
        Psig_recv, actual_sigma);

%% ----------------- Matched filter (correlation) --------------------------
h = fliplr(s);                  % matched filter impulse response (real)
z = conv(received_noisy, h);    % matched filter output (length = rec_len + N - 1)
t_z = (0:length(z)-1) / Fs;     % time axis for z

% Find peaks: we expect multiple peaks now due to multipath
[pk_vals, pk_locs] = findpeaks(abs(z), 'SortStr', 'descend', 'NPeaks', num_mp+3, 'MinPeakHeight', max(abs(z))*0.02);

% Convert first (largest) peak to range estimate
if ~isempty(pk_locs)
    idx_max = pk_locs(1);
    est_delay_samples = idx_max - N;    % account for convolution shift
    est_tau = est_delay_samples / Fs;
    R_est = c * est_tau / 2;
    fprintf('Estimated (strongest) delay samples = %d -> est_tau = %.3e s -> R_est = %.3f m\n',...
        est_delay_samples, est_tau, R_est);
else
    error('No peaks found in matched filter output.');
end

% For debugging: print all detected peaks and their estimated ranges
fprintf('Detected peaks (value, sample index, est range m):\n');
for i=1:length(pk_locs)
    dsamp = pk_locs(i) - N;
    rr = c * (dsamp / Fs) / 2;
    fprintf('  Peak %d: val=%.3e, idx=%d, range=%.3f m\n', i, pk_vals(i), pk_locs(i), rr);
end

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
f = (-N/2:N/2-1)*(Fs/N);
S = fftshift(fft(s_complex));
plot(f/1e3, abs(S));
xlabel('Frequency (kHz)'); ylabel('Magnitude');
title('Transmitted LFM Chirp - Frequency Domain');
xlim([(f0-B/2)/1e3 (f0+B*1.5)/1e3]);
grid on;

% Plot 3: Instantaneous frequency of transmitted signal
subplot(3,2,3);
f_inst = f0 + K*t;
plot(t*1e6, f_inst/1e3);
xlabel('Time (µs)'); ylabel('Frequency (kHz)');
title('Instantaneous Frequency of LFM Chirp');
grid on;

% Plot 4: Received signal (real) in time domain
subplot(3,2,4);
plot((0:length(received_noisy)-1)/Fs*1e6, received_noisy);
xlabel('Time (µs)'); ylabel('Amplitude');
title('Received Signal (delayed, multipath & noisy) - Time Domain');
xlim([0 (max_delay_samples + N + 50)/Fs*1e6]);
grid on;

% Plot 5: Received signal (complex) in frequency domain
subplot(3,2,5);
S_received = fftshift(fft(received_complex_noisy));
f_received = (-length(S_received)/2:length(S_received)/2-1)*(Fs/length(S_received));
plot(f_received/1e3, abs(S_received));
xlabel('Frequency (kHz)'); ylabel('Magnitude');
title('Received Signal - Frequency Domain');
xlim([(f0-B/2)/1e3 (f0+B*1.5)/1e3]);
grid on;

% Plot 6: Matched filter output (magnitude) with detected peaks
subplot(3,2,6);
plot(t_z*1e6, abs(z));
xlabel('Time (µs)'); ylabel('|z(t)|');
title('Matched Filter Output |z(t)| (correlation)');
grid on;
hold on;
% mark detected peaks
for i=1:length(pk_locs)
    plot((pk_locs(i))/Fs*1e6, pk_vals(i), 'ro', 'MarkerSize',8, 'LineWidth',1.2);
    text((pk_locs(i))/Fs*1e6, pk_vals(i)*0.9, sprintf('R=%.1f m', c*(pk_locs(i)-N)/(2*Fs)));
end
% highlight strongest
plot(pk_locs(1)/Fs*1e6, pk_vals(1), 'kp', 'MarkerSize',10, 'MarkerFaceColor','y');

%% ----------------- Detailed Frequency Analysis -----------------
figure('Name','Detailed Frequency Analysis','Position',[200 200 1000 700]);

% Spectrogram of transmitted signal
subplot(2,2,1);
window = 256;
noverlap = 200;
nfft = 1024;
spectrogram(s_complex, window, noverlap, nfft, Fs, 'yaxis');
title('Spectrogram of Transmitted LFM Chirp');
colorbar;

% Frequency content comparison (normalized)
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

% Quadratic phase characteristic (should be parabolic)
subplot(2,2,4);
phase_quadratic = phase_unwrapped - 2*pi*f0*t;
plot(t*1e6, phase_quadratic);
xlabel('Time (µs)'); ylabel('Phase - 2πf₀t (rad)');
title('Quadratic Phase Component (Kπt²)');
grid on;

%% ----------------- Matched Filter Peak Zoom -----------------
figure('Name','Matched Filter Peak Zoom','Position',[200 200 500 350]);
win = max(1, idx_max - 200) : min(length(z), idx_max + 200);
plot((win - N)/Fs*1e6, abs(z(win)));
xlabel('Delay (µs) relative to pulse'); ylabel('|z|'); title('Matched filter peak (zoom)');
grid on;
xline(est_delay_samples / Fs * 1e6, '--r', sprintf('Est delay = %.2f \\mus', est_tau*1e6));

%% ----------------- Ranging estimation explanation & performance ----------
range_resolution = c / (2*B);
fprintf('Theoretical range resolution (c / (2B)) = %.3f m\n', range_resolution);

% Range error (strongest detected)
range_error = R_est - R_true;
fprintf('Range estimate error (R_est - R_true) = %.3f m\n', range_error);

%% ----------------- Additional Analysis: Doppler Effect -----------------
figure('Name','Doppler Effect Analysis','Position',[300 300 1000 600]);

% Plot 1: Doppler phase shift over time
subplot(2,3,1);
doppler_phase = 2*pi*f_d*t;
plot(t*1e6, doppler_phase);
xlabel('Time (µs)'); ylabel('Doppler Phase (rad)');
title('Doppler Phase Shift vs Time');
grid on;

% Plot 2: Frequency shift due to Doppler
subplot(2,3,2);
f_original = f0 + K*t;
f_doppler = f_original + f_d;
plot(t*1e6, f_original/1e3, 'b-', 'LineWidth', 2, 'DisplayName', 'Original');
hold on;
plot(t*1e6, f_doppler/1e3, 'r--', 'LineWidth', 2, 'DisplayName', 'With Doppler');
xlabel('Time (µs)'); ylabel('Frequency (kHz)');
title('Instantaneous Frequency: Original vs Doppler');
legend;
grid on;

% Plot 3: Multipath range distribution
subplot(2,3,3);
multipath_ranges = R_true + mp_extra_range;
stem([R_true, multipath_ranges], [alpha, mp_alpha], 'filled', 'LineWidth', 2);
xlabel('Range (m)'); ylabel('Attenuation');
title('Multipath Range Distribution');
grid on;

% Plot 4: Phase shifts of multipaths
subplot(2,3,4);
stem(1:num_mp, mp_phase, 'filled', 'LineWidth', 2);
xlabel('Multipath Index'); ylabel('Phase Shift (rad)');
title('Multipath Phase Shifts');
grid on;
xticks(1:num_mp);

% Plot 5: Matched filter output around main peak (complex)
subplot(2,3,5);
win_range = max(1, idx_max-100):min(length(z), idx_max+100);
plot((win_range - N)/Fs*1e6, real(z(win_range)), 'b-', 'LineWidth', 1.5, 'DisplayName', 'Real');
hold on;
plot((win_range - N)/Fs*1e6, imag(z(win_range)), 'r-', 'LineWidth', 1.5, 'DisplayName', 'Imag');
plot((win_range - N)/Fs*1e6, abs(z(win_range)), 'k-', 'LineWidth', 2, 'DisplayName', 'Magnitude');
xlabel('Delay (µs)'); ylabel('Amplitude');
title('Complex MF Output around Main Peak');
legend;
grid on;

% Plot 6: SNR vs Detection Performance
subplot(2,3,6);
peak_magnitudes = pk_vals / max(pk_vals);
stem(1:length(peak_magnitudes), peak_magnitudes, 'filled', 'LineWidth', 2);
xlabel('Peak Index'); ylabel('Normalized Magnitude');
title('Normalized Peak Magnitudes');
grid on;

fprintf('\nSimulation completed successfully!\n');
fprintf('Summary: %d multipaths + Doppler (%.1f Hz) + AWGN (%.1f dB SNR)\n', ...
    num_mp, f_d, SNR_dB);