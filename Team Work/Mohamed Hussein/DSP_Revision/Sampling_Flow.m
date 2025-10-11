%% DSP Sampling, DTFT, and DFT Visualization
% Author: Mohamed Ahmed (Grad Project)
% Description: This script shows the full chain:
% Continuous signal → Sampling → Impulse train → DTFT → DFT → IDFT

clear; clc; close all;

%% 1️⃣ Generate a "continuous-time" signal
% In practice, we simulate continuous-time using a very fine sampling rate
Fs_cont = 10000;             % very high sampling frequency (simulating continuous)
t_cont = 0:1/Fs_cont:0.01;   % 10 ms duration
x_cont = sin(2*pi*500*t_cont) + 0.5*cos(2*pi*1200*t_cont); % example signal

figure;
plot(t_cont*1000, x_cont, 'b', 'LineWidth', 1.5); hold on;
title('Original Continuous-Time Signal');
xlabel('Time (ms)'); ylabel('Amplitude');
grid on;

%% 2️⃣ Sampling Process
Fs = 4000;                   % actual sampling frequency (must satisfy Nyquist > 2*fmax)
Ts = 1/Fs;
n = 0:Ts:0.01;
x_samp = sin(2*pi*500*n) + 0.5*cos(2*pi*1200*n);

% Create impulse train (for visualization only)
impulse_train = zeros(size(t_cont));
for i = 1:length(n)
    [~, idx] = min(abs(t_cont - n(i)));
    impulse_train(idx) = 1;
end

% Multiply the continuous signal by the impulse train
x_sampled_cont = x_cont .* impulse_train;

% Plot continuous, impulses, and sampled points
figure;
plot(t_cont*1000, x_cont, 'b', 'LineWidth', 1.2); hold on;
stem(n*1000, x_samp, 'r', 'filled');
stem(t_cont*1000, impulse_train*max(x_cont), 'k:');
legend('x(t) continuous','Sampled values','Impulse train');
title('Sampling in Time Domain');
xlabel('Time (ms)'); ylabel('Amplitude');
grid on;

%% 3️⃣ Frequency-domain view (conceptual)
% The multiplication in time by impulse train causes replication in frequency.
% We'll compute FFT of continuous and sampled (for illustration).

Nfft = 8192;
f_cont = linspace(-Fs_cont/2, Fs_cont/2, Nfft);
X_cont = fftshift(abs(fft(x_cont, Nfft)));

f_samp = linspace(-1*(Fs/2), 1*(Fs/2), Nfft);
X_samp = fftshift(abs(fft(x_samp, Nfft)));

figure;
subplot(2,1,1);
plot(f_cont, X_cont / max(X_cont), 'b');
title('Continuous-Time Spectrum');
xlabel('Frequency (Hz)'); ylabel('Normalized |X(f)|'); grid on;

subplot(2,1,2);
plot(f_samp, X_samp / max(X_samp), 'r');
title('Sampled Signal Spectrum (Periodic every Fs)');
xlabel('Frequency (Hz)'); ylabel('Normalized |X_s(f)|'); grid on;

%% 4️⃣ DTFT vs DFT
% The DTFT is continuous and periodic in frequency.
% We can only approximate it numerically (cannot compute truly continuous DTFT).
% That’s why the DTFT is not directly applicable on computers.

% Approximate DTFT by zero-padding FFT of the discrete signal
N_dtft = 2048;
omega = linspace(-pi, pi, N_dtft);
X_dtft = fftshift(abs(fft(x_samp, N_dtft)));

figure;
plot(omega, X_dtft / max(X_dtft), 'LineWidth', 1.5);
title('DTFT (approximated via zero-padded FFT)');
xlabel('\omega (rad/sample)'); ylabel('|X(e^{j\omega})|');
grid on;

% DFT (finite length)
N = 64;
X_dft = fft(x_samp, N);
k = 0:N-1;
omega_k = 2*pi*k/N;

figure;
stem(omega_k, abs(X_dft)/max(abs(X_dft)), 'filled');
title('DFT Samples (Discrete Frequency Points)');
xlabel('\omega (rad/sample)'); ylabel('|X[k]|');
grid on;

%% 5️⃣ IDFT Reconstruction
x_idft = ifft(X_dft, N);

figure;
subplot(2,1,1);
stem(0:N-1, real(x_idft), 'filled');
title('Reconstructed Signal from IDFT (Periodic in Time)');
xlabel('n'); ylabel('Amplitude'); grid on;

subplot(2,1,2);
stem(0:length(x_samp)-1, x_samp(1:length(x_samp)), 'r');
title('Original Sampled Signal');
xlabel('n'); ylabel('Amplitude'); grid on;

%% 6️⃣ Notes (Important)
% - The DTFT is *continuous* in frequency and *periodic* with 2π.
%   But it cannot be exactly computed — we approximate it using FFT.
%
% - The DFT assumes the time-domain signal is *periodic* with period N.
%   Its frequency spectrum is *discrete* (N points) and *periodic* with period N.
%
% - The IDFT reconstructs one period of the time-domain signal.
