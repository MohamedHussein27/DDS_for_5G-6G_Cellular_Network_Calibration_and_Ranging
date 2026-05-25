% Combined pipeline: DDS -> FFT -> MUX -> IFFT
% Only final IFFT output vectors are saved to .txt files.
clear; close all; clc;

%% ----------------- DDS: Transmitted chirp (hardware-equivalent) -----------------
c = 3e8;                 % speed of light (m/s)
f0 = 0;                 % start frequency (Hz)
B = 200e6;              % bandwidth (Hz)
T = 8.333e-6;           % pulse duration (s)
K = B / T;              % chirp rate (Hz/s)
theta = 0;              % initial phase (rad)
Fs = 491.52e6;          % sampling frequency (Hz)

N = round(Fs * T);      % number of samples, here N = 4096
n = 1:N;                % discrete clock cycles

% Phase accumulation matching RTL accumulator (n*(n-1)/2 instead of t^2)
phase_rad = 2 * pi * ( (f0/Fs).*n + (B/(Fs*N)) .* (n .* (n-1) / 2) );

A = 127 / 256;          % amplitude ~0.496 (matches hardware Quarter-Wave LUT)
s = A * sin(phase_rad + theta);   % DDS output vector

%% ----------------- FFT: Built-in floating-point FFT -----------------
fft_out = fft(s);               % complex FFT result
fft_complex = fft_out;          % keep complex directly for MUX

%% ----------------- MUX: Combine OFDM and radar/FFT channels -----------------
% Load OFDM symbols (must exist as text files)
ofdm_re = load('ofdm_data_float_re.txt');
ofdm_im = load('ofdm_data_float_im.txt');
ofdm_complex = ofdm_re + 1j * ofdm_im;

N_OFDM   = 2048;
N_ZEROS  = 381;
N_RADAR  = 1667;                % = N - (N_OFDM + N_ZEROS)

% Build 4096-point MUX frame
mux_frame = zeros(N, 1);

% 1) OFDM occupies first 2048 bins
mux_frame(1 : N_OFDM) = ofdm_complex(1 : N_OFDM);

% 2) Chirp (radar) occupies bins 2430 to 4096, scaled by 1/128 (hardware shift)
radar_start_idx = N_OFDM + N_ZEROS + 1;   % 2048+381+1 = 2430
mux_frame(radar_start_idx : end) = fft_complex(1 : N_RADAR) / 128.0;

%% ----------------- IFFT: Final time-domain transform -----------------
ifft_out = ifft(mux_frame);

% MATLAB's ifft divides by N; hardware typically does not, so multiply back.
ifft_unscaled = ifft_out * N;

ifft_re = real(ifft_unscaled);
ifft_im = imag(ifft_unscaled);

%% ----------------- Export final IFFT output vectors only -----------------
writematrix(ifft_re, 'matlab_out_re_float.txt');
writematrix(ifft_im, 'matlab_out_im_float.txt');

disp('Pipeline complete. Final IFFT outputs saved as:');
disp('  matlab_out_re_float.txt');
disp('  matlab_out_im_float.txt');