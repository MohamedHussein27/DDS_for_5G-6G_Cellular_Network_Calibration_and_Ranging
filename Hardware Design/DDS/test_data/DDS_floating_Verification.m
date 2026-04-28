% Ideal Mathematical Chirp vs RTL Comparison (Corrected for FTW_start)
clear; clc;

% --- Hardware Parameters ---
Fs = 491.52e6;             % System clock
N_cycles = 4096;           % Number of simulated cycles
T_dur = N_cycles / Fs;     % Total duration in seconds
FTW_step = 0;         % Hardware FTW step (Must match RTL testbench)
FTW_start = 174762666;             % Hardware FTW start (Must match RTL testbench)

% Calculate the exact physical bandwidth reached by the hardware
max_ftw = FTW_start + (FTW_step * N_cycles);
B = ((max_ftw-FTW_start) * Fs) / (2^32); 

% --- 1. Load RTL Output ---
try
    rtl_data = load('rtl_output.txt');
catch
    error('Could not find rtl_output.txt. Run the Verilog testbench first.');
end

% --- 2. Align Pipeline Latency ---
latency = 1; 

% Extract the valid window of RTL data to match the ideal 0-indexed vectors
% We grab N_cycles of data starting after the latency zeros.
rtl_aligned = rtl_data(1+latency : N_cycles+latency)';

% --- 3. Generate Ideal Discrete-Time Model ---
n = 0:(N_cycles-1); 
t = n / Fs;

% Discrete Phase Accumulation: 
% Replicates the exact cycle-by-cycle addition in the Verilog PHASE_ACC.
% Now includes the initial FTW_start offset!
% We use mod(..., 2^32) to simulate the 32-bit register overflow.
ideal_discrete_phase_word = mod((FTW_start .* n) + (FTW_step * (n .* (n - 1)) / 2), 2^32);

% Convert 32-bit tuning word phase to radians
ideal_phase_rad = 2 * pi * (ideal_discrete_phase_word / 2^32);

% Ideal Amplitude (Scaled to match 8-bit signed hardware, max = 127)
max_amp = 127;
ideal_out = max_amp * sin(ideal_phase_rad);

% --- 4. Plotting Time Domain ---
figure('Name', 'DDS Output: RTL vs Ideal');

% Subplot 1: The Waveforms
subplot(2,1,1);
plot(t * 1e6, ideal_out, 'b', 'LineWidth', 1.5); hold on;
plot(t * 1e6, rtl_aligned, 'r--', 'LineWidth', 1.5);
title('LFM Chirp: RTL vs Discrete Mathematical Model');
xlabel('Time (\mus)');
ylabel('Amplitude');
legend('Ideal Discrete Model', 'Verilog RTL (Aligned)');
grid on;

% Subplot 2: The Error (Quantization & Phase Noise)
subplot(2,1,2);
noise_signal = ideal_out - rtl_aligned;
plot(t * 1e6, noise_signal, 'k');
title('Hardware Error (Quantization Noise & Phase Truncation)');
xlabel('Time (\mus)');
ylabel('Amplitude Error (LSBs)');
grid on;

% --- 5. Signal Quality Metrics ---
% Calculate the Signal-to-Noise Ratio introduced by the hardware limits
signal_power = rms(ideal_out)^2;
noise_power = rms(noise_signal)^2;
sqnr_db = 10 * log10(signal_power / noise_power);

fprintf('Calculated Physical Bandwidth: %.2f MHz\n', B / 1e6);
fprintf('Measured Hardware SQNR: %.2f dB\n', sqnr_db);

% --- 6. Frequency Spectrum (FFT) ---
% Compute the Fast Fourier Transform (FFT) of the aligned RTL signal
Y = fft(rtl_aligned);

% Calculate the two-sided spectrum P2, then compute the single-sided spectrum P1
P2 = abs(Y / N_cycles);
P1 = P2(1:N_cycles/2+1);
P1(2:end-1) = 2 * P1(2:end-1); % Double the amplitude to account for positive/negative frequencies

% Convert to Decibels Full Scale (dBFS), relative to the maximum 8-bit amplitude (127)
% eps is added to prevent log10(0) errors
P1_dBFS = 20 * log10(P1 / max_amp + eps); 

% Define the frequency domain axis (0 to Nyquist frequency)
f = Fs * (0:(N_cycles/2)) / N_cycles;

% Plotting the Frequency Domain
figure('Name', 'RTL Output Frequency Spectrum');
plot(f / 1e6, P1_dBFS, 'b', 'LineWidth', 1.2); hold on;
title('Single-Sided Amplitude Spectrum of RTL Chirp');
xlabel('Frequency (MHz)');
ylabel('Magnitude (dBFS)');
ylim([-100 0]); % Lock the Y-axis to standard dBFS range
grid on;

% Add a vertical line to visually verify the chirp's bandwidth matches the math
xline(B / 1e6, 'r--', 'Calculated Max Bandwidth', 'LabelHorizontalAlignment', 'left', 'LineWidth', 1.5);