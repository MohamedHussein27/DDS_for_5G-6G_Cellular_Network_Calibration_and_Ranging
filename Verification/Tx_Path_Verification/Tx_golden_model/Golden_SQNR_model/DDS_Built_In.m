% LFM radar simulation with matched filtering (AWGN only) - WITH PSD
clear; close all; clc;

%% ----------------- Parameters -----------------
c = 3e8;                 % speed of light (m/s)
f0 = 0;                % start frequency (Hz) (baseband-ish example)
B = 200e6;                 % bandwidth (Hz)
T = 8.333e-6;              % pulse duration (seconds)
K = B / T;               % chirp rate (Hz/s) -> K = B/T
theta = 0;               % initial phase (rad)
Fs = 491.52e6;                % sampling frequency (Hz) (must be >> f0+B)
SNR_dB = 10;             % desired SNR at the receiver (dB) for the echo
R_true = 300;            % true target range in meters
alpha = 0.5;             % target amplitude attenuation (reflection coef)

%% ----------------- Discrete Hardware Time vectors -----------------
N = round(Fs * T);       % Equals exactly 4096
n = 1:N;                 % Discrete clock cycles (matches RTL i = 1 to N_cycles)

%% ----------------- Transmitted signal (Hardware Equivalent) -------
% Instead of t^2, hardware accumulates as n*(n-1)/2
% We divide by 2*pi to normalize the frequency tuning words exactly like the RTL
phase_rad = 2 * pi * ( (f0/Fs).*n + (B/(Fs*N)) .* (n .* (n-1) / 2) );

A = 127 / 256;  % Roughly 0.496
% Use sin() to match the RTL's Quarter-Wave LUT!
s = A * sin(phase_rad + theta); 

% Save to a text file so Python can read it
writematrix(s.', 'matlab_float_out.txt');