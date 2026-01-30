% LFM radar simulation with multiple pulses and PRF - WITH PSD
clear; close all; clc;

%% ----------------- Parameters -----------------
c = 3e8;                 % speed of light (m/s)
A = 1;                   % amplitude of transmitted chirp
f0 = 0;                  % start frequency (Hz) (baseband-ish example)
B = 80e6;                % bandwidth (Hz)
T = 50e-6;               % pulse duration (seconds)
K = B / T;               % chirp rate (Hz/s) -> K = B/T
theta = 0;               % initial phase (rad)
Fs = 5e8;                % sampling frequency (Hz) (must be >> f0+B)
SNR_dB = 10;             % desired SNR at the receiver (dB) for the echo
R_true = 300;            % true target range in meters
alpha = 0.5;             % target amplitude attenuation (reflection coef)

%% ----------------- Radar Pulse Parameters -----------------
PRF = 1000;              % Pulse Repetition Frequency (Hz)
PRI = 1/PRF;             % Pulse Repetition Interval (seconds)
num_pulses = 3;          % Number of pulses to transmit
fprintf('PRF = %d Hz, PRI = %.3f ms\n', PRF, PRI*1e3);

%% ----------------- Single Pulse Generation -----------------
t_pulse = 0:1/Fs:(T - 1/Fs);   % time vector for one pulse
N_pulse = length(t_pulse);

% Create single LFM pulse
s_pulse = A*cos(2*pi*(f0*t_pulse + 0.5*K*t_pulse.^2) + theta);
s_pulse_complex = A*exp(1j*(2*pi*(f0.*t_pulse + 0.5*K.*t_pulse.^2) + theta));

%% ----------------- Transmit Multiple Pulses -----------------
% Calculate samples per PRI
samples_per_PRI = round(PRI * Fs);
fprintf('Samples per pulse: %d\n', N_pulse);
fprintf('Samples per PRI: %d\n', samples_per_PRI);

% Create transmitted signal with multiple pulses
tx_signal = zeros(1, num_pulses * samples_per_PRI);
tx_signal_complex = zeros(1, num_pulses * samples_per_PRI);

for pulse_idx = 1:num_pulses
    start_sample = (pulse_idx-1) * samples_per_PRI + 1;
    end_sample = start_sample + N_pulse - 1;
    tx_signal(start_sample:end_sample) = s_pulse;
    tx_signal_complex(start_sample:end_sample) = s_pulse_complex;
end

fprintf('Transmitted signal: %d pulses, total length: %d samples (%.3f ms)\n', ...
        num_pulses, length(tx_signal), length(tx_signal)/Fs*1e3);

%% ----------------- Create Received Echoes -----------------
tau = 2*R_true / c;                     % round-trip delay (s)
delay_samples = round(tau * Fs);        % integer sample delay
fprintf('True round-trip delay = %.3e s -> delay samples = %d\n', tau, delay_samples);

% Received vector length: long enough to contain all echoes
rec_len = length(tx_signal) + delay_samples + N_pulse;
received = zeros(1, rec_len);
received_complex = zeros(1, rec_len);

% Place delayed, attenuated pulses in received signal
for pulse_idx = 1:num_pulses
    tx_start = (pulse_idx-1) * samples_per_PRI + 1;
    rx_start = tx_start + delay_samples;
    rx_end = rx_start + N_pulse - 1;
    
    if rx_end <= rec_len
        received(rx_start:rx_end) = received(rx_start:rx_end) + alpha * s_pulse;
        received_complex(rx_start:rx_end) = received_complex(rx_start:rx_end) + alpha * s_pulse_complex;
    end
end

%% ----------------- Add AWGN -----------------
% Compute signal power for SNR calculation
Psig = mean((alpha*s_pulse).^2);
SNR_lin = 10^(SNR_dB/10);
noise_variance = Psig / SNR_lin;
sigma = sqrt(noise_variance);

% Add white Gaussian noise
noise_real = sigma * randn(size(received));
received_noisy = received + noise_real;

noise_complex = sigma/sqrt(2) * (randn(size(received)) + 1j*randn(size(received)));
received_complex_noisy = received_complex + noise_complex;

%% ----------------- Matched Filter Processing -----------------
% Matched filter for single pulse
h = fliplr(s_pulse);                  % matched filter impulse response

% Process each pulse separately using pulse-Doppler processing
pulse_starts = 1:samples_per_PRI:(num_pulses-1)*samples_per_PRI + 1;
range_profiles = zeros(num_pulses, samples_per_PRI + N_pulse - 1);

figure('Name','Pulse-Doppler Processing','Position',[100 100 1200 800]);

for pulse_idx = 1:num_pulses
    % Extract segment for this pulse (include time for maximum possible echo)
    start_segment = pulse_starts(pulse_idx);
    end_segment = min(start_segment + samples_per_PRI + N_pulse - 2, length(received_noisy));
    segment = received_noisy(start_segment:end_segment);
    
    % Matched filtering
    z = conv(segment, h);
    range_profiles(pulse_idx, 1:length(z)) = abs(z);
    
    % Time axis for this pulse
    t_z = (0:length(z)-1) / Fs * 1e6;
    
    % Find peak and estimate range
    [peak_val, idx_max] = max(abs(z));
    est_delay_samples = idx_max - N_pulse;
    est_tau = est_delay_samples / Fs;
    R_est = c * est_tau / 2;
    
    % Plot individual pulse processing
    subplot(2, 2, pulse_idx);
    plot(t_z, abs(z), 'b', 'LineWidth', 1);
    hold on;
    plot(t_z(idx_max), peak_val, 'ro', 'MarkerSize', 8, 'LineWidth', 2);
    xlabel('Time (µs)'); ylabel('|Matched Filter Output|');
    title(sprintf('Pulse %d: Estimated Range = %.1f m', pulse_idx, R_est));
    grid on;
    
    fprintf('Pulse %d: Estimated range = %.2f m (error = %.2f m)\n', ...
            pulse_idx, R_est, abs(R_est - R_true));
end

% Plot all range profiles together
subplot(2,2,4);
for pulse_idx = 1:num_pulses
    t_profile = (0:size(range_profiles,2)-1) / Fs * 1e6;
    plot(t_profile, range_profiles(pulse_idx,:), 'LineWidth', 1.5, ...
         'DisplayName', sprintf('Pulse %d', pulse_idx));
    hold on;
end
xlabel('Time (µs)'); ylabel('|Matched Filter Output|');
title('All Range Profiles');
legend('show'); grid on;

%% ----------------- PSD Analysis -----------------
% Calculate PSD using pwelch method
[pxx_trans, f_trans] = pwelch(tx_signal_complex, hanning(1024), 512, 1024, Fs, 'centered', 'psd');
[pxx_rec, f_rec] = pwelch(received_complex_noisy, hanning(1024), 512, 1024, Fs, 'centered', 'psd');
[pxx_noise, f_noise] = pwelch(noise_complex, hanning(1024), 512, 1024, Fs, 'centered', 'psd');

% Convert to dB scale
pxx_trans_db = 10*log10(pxx_trans);
pxx_rec_db = 10*log10(pxx_rec);
pxx_noise_db = 10*log10(pxx_noise);

%% ----------------- PLOTS: Time Domain Signals -----------------
figure('Name','Multiple Pulse LFM Radar','Position',[100 50 1400 900]);

% Plot 1: Transmitted signal (first few pulses)
subplot(3,2,1);
t_tx = (0:length(tx_signal)-1) / Fs * 1e3;  % Time in milliseconds
plot(t_tx, tx_signal);
xlabel('Time (ms)'); ylabel('Amplitude'); 
title('Transmitted Signal (Multiple LFM Pulses)');
grid on;
xlim([0 min(5, max(t_tx))]); % Show first 5 ms
% Mark pulse locations
hold on;
for pulse_idx = 1:num_pulses
    pulse_time = (pulse_idx-1) * PRI * 1e3;
    xline(pulse_time, '--r', 'Alpha', 0.5);
    xline(pulse_time + T*1e3, '--r', 'Alpha', 0.5);
end
hold off;

% Plot 2: Single pulse zoom
subplot(3,2,2);
plot(t_pulse*1e6, s_pulse);
xlabel('Time (µs)'); ylabel('Amplitude'); 
title('Single LFM Pulse (Zoom)');
grid on;

% Plot 3: Instantaneous frequency of single pulse
subplot(3,2,3);
f_inst = f0 + K*t_pulse;
plot(t_pulse*1e6, f_inst/1e6);
xlabel('Time (µs)'); ylabel('Frequency (MHz)');
title('Instantaneous Frequency of LFM Chirp');
grid on;

% Plot 4: Received signal with echoes
subplot(3,2,4);
t_rx = (0:length(received_noisy)-1) / Fs * 1e3;  % Time in milliseconds
plot(t_rx, received_noisy);
xlabel('Time (ms)'); ylabel('Amplitude'); 
title('Received Signal (with Multiple Echoes)');
grid on;
xlim([0 min(5, max(t_rx))]);
% Mark echo locations
hold on;
for pulse_idx = 1:num_pulses
    echo_start = ((pulse_idx-1) * PRI + tau) * 1e3;
    echo_end = echo_start + T*1e3;
    xline(echo_start, '--g', 'Alpha', 0.7, 'LineWidth', 1.5);
    xline(echo_end, '--g', 'Alpha', 0.7, 'LineWidth', 1.5);
end
legend('Received', 'Echo start/end', 'Location', 'best');
hold off;

% Plot 5: Frequency domain - transmitted
subplot(3,2,5);
S_tx = fftshift(fft(tx_signal_complex));
f_tx = (-length(S_tx)/2:length(S_tx)/2-1)*(Fs/length(S_tx));
plot(f_tx/1e6, abs(S_tx));
xlabel('Frequency (MHz)'); ylabel('Magnitude');
title('Transmitted Signal Spectrum');
xlim([-B/1e6*1.2, B/1e6*1.2]);
grid on;

% Plot 6: Frequency domain - received
subplot(3,2,6);
S_rx = fftshift(fft(received_complex_noisy));
f_rx = (-length(S_rx)/2:length(S_rx)/2-1)*(Fs/length(S_rx));
plot(f_rx/1e6, abs(S_rx));
xlabel('Frequency (MHz)'); ylabel('Magnitude');
title('Received Signal Spectrum');
xlim([-B/1e6*1.2, B/1e6*1.2]);
grid on;

%% ----------------- PSD Analysis Plots -----------------
figure('Name','Power Spectral Density Analysis','Position',[300 100 1200 800]);

% Plot 1: PSD comparison (dB scale)
subplot(2,2,1);
plot(f_trans/1e6, pxx_trans_db, 'b', 'LineWidth', 1.5, 'DisplayName', 'Transmitted');
hold on;
plot(f_rec/1e6, pxx_rec_db, 'r', 'LineWidth', 1, 'DisplayName', 'Received');
plot(f_noise/1e6, pxx_noise_db, 'g', 'LineWidth', 0.5, 'DisplayName', 'Noise');
xlabel('Frequency (MHz)'); ylabel('PSD (dB/Hz)');
title('PSD Comparison - Multiple Pulses');
legend;
grid on;
xlim([-B/1e6*1.2, B/1e6*1.2]);

% Plot 2: Spectrogram of transmitted signal
subplot(2,2,2);
spectrogram(tx_signal_complex, 256, 200, 1024, Fs, 'yaxis');
title('Spectrogram: Transmitted Multiple Pulses');
colorbar;

% Plot 3: Spectrogram of received signal
subplot(2,2,3);
spectrogram(received_complex_noisy, 256, 200, 1024, Fs, 'yaxis');
title('Spectrogram: Received Signal with Echoes');
colorbar;

% Plot 4: Range-Doppler map (simplified)
subplot(2,2,4);
imagesc((0:size(range_profiles,2)-1)/Fs*1e6, 1:num_pulses, range_profiles);
xlabel('Range (µs)'); ylabel('Pulse Number');
title('Range Profiles vs Pulse Number');
colorbar;
xlim([0 tau*2e6]); % Focus on range around target

%% ----------------- Performance Analysis -----------------
fprintf('\n=== Performance Summary ===\n');
fprintf('Pulse Duration: %.1f µs\n', T*1e6);
fprintf('PRI: %.1f µs\n', PRI*1e6);
fprintf('Duty Cycle: %.1f%%\n', (T/PRI)*100);
fprintf('Unambiguous Range: %.1f km\n', (c * PRI / 2) / 1000);
fprintf('Range Resolution: %.2f m\n', c/(2*B));

% Calculate average range estimate
final_ranges = zeros(1, num_pulses);
for pulse_idx = 1:num_pulses
    z = range_profiles(pulse_idx, :);
    [~, idx_max] = max(z);
    est_delay_samples = idx_max - N_pulse;
    final_ranges(pulse_idx) = c * (est_delay_samples / Fs) / 2;
end

fprintf('Average estimated range: %.2f m\n', mean(final_ranges));
fprintf('Range estimation std: %.2f m\n', std(final_ranges));