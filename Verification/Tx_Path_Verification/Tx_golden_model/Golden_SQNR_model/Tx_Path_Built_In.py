import numpy as np
import os

# =============================================================================
# DDS: Transmitted chirp (hardware-equivalent)
# =============================================================================
c = 3e8                 # speed of light (m/s)
f0 = 0.0                # start frequency (Hz)
B = 200e6               # bandwidth (Hz)
T = 8.333e-6            # pulse duration (s)
K = B / T               # chirp rate (Hz/s)
theta = 0.0             # initial phase (rad)
Fs = 491.52e6           # sampling frequency (Hz)

N = int(round(Fs * T))  # number of samples (should be 4096)
n = np.arange(1, N + 1) # discrete clock cycles (1-indexed like RTL)

# Phase accumulation: n*(n-1)/2 exactly as in hardware
phase_rad = 2 * np.pi * ((f0 / Fs) * n +
                         (B / (Fs * N)) * (n * (n - 1) / 2))

A = 127 / 256           # amplitude ~0.496
s = A * np.sin(phase_rad + theta)  # DDS output vector

# =============================================================================
# FFT: Built-in floating-point FFT
# =============================================================================
fft_out = np.fft.fft(s)          # complex FFT result
fft_complex = fft_out            # keep complex directly for MUX

# =============================================================================
# MUX: Combine OFDM and radar/FFT channels
# =============================================================================
# Load OFDM symbols (must exist as text files)
ofdm_re = np.loadtxt('ofdm_data_float_re.txt')
ofdm_im = np.loadtxt('ofdm_data_float_im.txt')
ofdm_complex = ofdm_re + 1j * ofdm_im

N_OFDM   = 2048
N_ZEROS  = 381
N_RADAR  = 1667          # = N - (N_OFDM + N_ZEROS)

# Build 4096-point MUX frame
mux_frame = np.zeros(N, dtype=complex)

# 1) OFDM occupies first 2048 bins
mux_frame[:N_OFDM] = ofdm_complex[:N_OFDM]

# 2) Chirp (radar) occupies bins 2430 to 4096 (0-indexed: 2429..end)
#    scaled by 1/128 to match hardware arithmetic shift
radar_start_idx = N_OFDM + N_ZEROS          # 2429 in 0-indexed
mux_frame[radar_start_idx:] = fft_complex[:N_RADAR] / 128.0

# =============================================================================
# IFFT: Final time-domain transform
# =============================================================================
ifft_out = np.fft.ifft(mux_frame)

# MATLAB's ifft divides by N; hardware typically does not, so multiply back
ifft_unscaled = ifft_out * N

ifft_re = np.real(ifft_unscaled)
ifft_im = np.imag(ifft_unscaled)

# =============================================================================
# Export final IFFT output vectors only
# =============================================================================
np.savetxt('py_gold_out_re_float.txt', ifft_re, fmt='%.8e')
np.savetxt('py_gold_out_im_float.txt', ifft_im, fmt='%.8e')

print("Pipeline complete. Final IFFT outputs saved as:")
print("  py_gold_out_re_float.txt")
print("  py_gold_out_im_float.txt")
