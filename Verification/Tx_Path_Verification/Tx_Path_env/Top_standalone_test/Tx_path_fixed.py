"""
system_builtin.py
=================
Python golden model — equivalent of SYSTEM_builtin.m

Full TX/Channel/RX pipeline for the Split-Band Radar-Communication system:
  - TX  : DDS Chirp + 256-QAM OFDM → FFT → MUX → IFFT
  - CH  : Circular shift delay (ideal channel)
  - RX  : FFT → DeMUX → Correlation → IFFT (range profile)
  - Metrics : SQNR, power balance

Dependencies (same folder):
    mytypes.py
    dds_core.py
    fft_fixed.py   (contains radix22_dif_fft_fixed + radix22_dif_ifft_fixed)
    mux.py
"""

import numpy as np
from mytypes   import quantize_fixed, reinterpretcast
from dds_core  import dds_core
from fft_fixed import radix22_dif_fft_fixed
from ifft_fixed import radix22_dif_ifft_fixed
from mux       import mux
import matplotlib.pyplot as plt

# ===========================================================================
# 0. WORD-LENGTH PARAMETERS
# ===========================================================================
WL_reduction = 0
scale        = 64           # Amplitude scale factor used in chirp reinterpret

WL_FFT  = 16 - WL_reduction;   FL_FFT  =  8 - WL_reduction
WL_FFT1 = 16 - WL_reduction;   FL_FFT1 =  8 - WL_reduction
WL_FFT2 = 16 - WL_reduction;   FL_FFT2 = 14 - WL_reduction
WL_IFFT2= 16 - WL_reduction;   FL_IFFT2=  5 - WL_reduction


# ===========================================================================
# 1. DESIGN PARAMETERS
# ===========================================================================
Fs      = 491.52e6        # System clock / sampling frequency (Hz)
T_dur   = 8.333e-6        # Signal duration (s)
Nacc    = 10              # Phase accumulator bits
LUT_bits= 10              # LUT address bits
Ns      = int(round(T_dur * Fs))  # Number of samples
t       = np.arange(Ns) / Fs     # Time vector

f0 = 0.0
B  = 200e6                # Chirp bandwidth

# Instantaneous frequency law
f_inst = f0 + B * (t / T_dur)

# Nyquist check
f_inst[f_inst >= Fs / 2] = Fs / 2 - 1e3

# Tuning word
M_dds = np.round(f_inst * (2 ** Nacc) / Fs).astype(np.uint64)


# ===========================================================================
# 2. FREQUENCY CONFIGURATION
# ===========================================================================
fs = 491.52e6             # Sampling frequency
N  = 4096                 # FFT size
df = fs / N               # Subcarrier spacing (~120 kHz)
M_qam = 256               # QAM order

# Frequency axis (centred)
f_axis = (np.arange(N) - N // 2) * df    # -245.76 MHz … +245.76 MHz

# Band limits
idx_chirp_start = np.searchsorted(f_axis, -200e6)
idx_center      = N // 2                  # index of 0 Hz  (0-based)
idx_ofdm_end    = np.searchsorted(f_axis,  200e6, side='right') - 1

# Allocation index vectors (0-based, like Python)
idx_chirp_active = np.arange(idx_chirp_start, idx_center)
idx_ofdm_active  = np.arange(idx_center, idx_ofdm_end + 1)


# ===========================================================================
# 3. CHIRP INDICES (standard / unshifted FFT order)
# ===========================================================================
# Negative frequencies in unshifted FFT live at indices N/2 … N-1
idx_start_chirp_std = N + int(round(-200e6 / df))  # 2429 (0-based)
idx_stop_chirp_std  = N - 1                         # 4095 (0-based)
idx_chirp_standard  = np.arange(idx_start_chirp_std, idx_stop_chirp_std + 1)


# ===========================================================================
# 4. TX PATH
# ===========================================================================
print("--- TX PATH ---")

# --- 4.1 DDS Chirp (time domain) ---
chirp_time = dds_core(M_dds, Nacc, LUT_bits, 'fixed').reshape(-1, 1)

# --- 4.2 FFT of chirp ---
chirp_freq_full = radix22_dif_fft_fixed(chirp_time.flatten(), WL_FFT1, FL_FFT1)

# Reinterpret cast: numerictype(1, WL_FFT, FL_FFT + log2(scale))
log2_scale = int(round(np.log2(scale)))
chirp_freq_full = reinterpretcast(
    chirp_freq_full,
    src_WL=WL_FFT1, src_FL=FL_FFT1,
    dst_WL=WL_FFT,  dst_FL=FL_FFT + log2_scale
)
chirp_freq_full = chirp_freq_full.flatten()

# fftshift for masking / visualisation
chirp_freq_shifted = np.fft.fftshift(chirp_freq_full)

# --- 4.3 OFDM subcarrier indices ---
f_start_ofdm = 10e6
f_stop_ofdm  = f_start_ofdm + 180e6    # 190 MHz

idx_start_ofdm = int(round(f_start_ofdm / df))     # 0-based
idx_stop_ofdm  = int(round(f_stop_ofdm  / df))     # 0-based
idx_ofdm_active_calc = np.arange(idx_start_ofdm, idx_stop_ofdm + 1)

# --- 4.4 256-QAM OFDM symbols ---
num_ofdm_subcarriers = len(idx_ofdm_active_calc)
data_integers = np.random.randint(0, M_qam, size=num_ofdm_subcarriers)

# Manual 256-QAM constellation (no toolbox required – mirrors MATLAB fallback)
I_grid, Q_grid = np.meshgrid(np.arange(-15, 16, 2), np.arange(-15, 16, 2))
constellation  = (I_grid + 1j * Q_grid).flatten()
constellation  = constellation / np.sqrt(np.mean(np.abs(constellation) ** 2))
ofdm_symbols   = constellation[data_integers]

# --- 4.5 Fixed-point cast for OFDM stream ---
def fi_cast(x, WL, FL):
    """Convergent rounding + Saturate – matches MATLAB fimath in SYSTEM_builtin."""
    return quantize_fixed(x, WL=WL, FL=FL, signed=True,
                          rounding='convergent', overflow='saturate')

# --- 4.6 Build stream_pos (OFDM, positive frequencies) ---
target_len    = N // 2   # 2048
stream_pos_in = ofdm_symbols.copy()

if len(stream_pos_in) < target_len:
    pad_len       = target_len - len(stream_pos_in)
    # 83 zero-bins guard at DC, remaining zeros at high end
    stream_pos_in = np.concatenate([
        np.zeros(83, dtype=complex),
        stream_pos_in,
        np.zeros(pad_len - 83, dtype=complex)
    ])

# --- 4.7 Build stream_neg (Chirp, negative frequencies) ---
cutoff_freq = 200e6
idx_cutoff  = int(round(cutoff_freq / (fs / N)))   # ~1667 (0-based → 1666)

chirp_freq_half = chirp_freq_full[:N // 2].copy()
chirp_freq_half[idx_cutoff:] = 0               # zero out >200 MHz bins

stream_neg_in = chirp_freq_half[:1667].copy()  # keep only 0..1666 (MATLAB 1..1667)

if len(stream_neg_in) < target_len:
    pad_len       = target_len - len(stream_neg_in)
    stream_neg_in = np.concatenate([
        np.zeros(pad_len, dtype=complex),
        stream_neg_in
    ])

# --- 4.8 MUX ---
X_tx_combined, tx_valid = mux(stream_pos_in, stream_neg_in, N)
X_tx_combined = fi_cast(X_tx_combined, WL_FFT, FL_FFT)

# --- 4.9 IFFT ---
x_tx_time, _, _ = radix22_dif_ifft_fixed(X_tx_combined, WL_FFT, FL_FFT)

print(f"  TX time-domain samples: {len(x_tx_time)}")


# ===========================================================================
# 5. TRUE HARDWARE SPECTRUM (ZERO-PADDING)
# ===========================================================================
print("--- GENERATING HIGH-RESOLUTION HARDWARE SPECTRUM ---")

# 1. Take the time-domain signal and pad it with zeros (4x oversampling)
pad_factor = 4

# Ensure x_tx_time is an array of standard complex floats for the np.fft function
try:
    x_time_float = np.array([complex(x) for x in x_tx_time])
except TypeError:
    x_time_float = np.array([x.get_val() if hasattr(x, 'get_val') else x for x in x_tx_time])

# Pad with zeros at the end
pad_len = (pad_factor - 1) * len(x_time_float)
x_time_padded = np.concatenate([x_time_float, np.zeros(pad_len, dtype=complex)])

# 2. Calculate the high-resolution FFT
N_pad = len(x_time_padded)
X_true_freq = np.fft.fftshift(np.fft.fft(x_time_padded))

# 3. Create a high-resolution frequency axis
f_axis_pad = (np.arange(N_pad) - N_pad // 2) * (fs / N_pad) / 1e6  # in MHz

# 4. Plot the true spectrum
plt.figure("The Hardware Reality: True Transmit Spectrum", figsize=(8, 5))
plt.plot(f_axis_pad, 20 * np.log10(np.abs(X_true_freq) + 1e-12), color='blue', linewidth=0.8)

plt.axvline(0, color='black', linestyle='--', linewidth=2, label='DC Boundary')
plt.title("True Interpolated Spectrum (Showing Out-of-Band Emissions)")
plt.xlabel("Frequency (MHz)")
plt.ylabel("Magnitude (dB)")
plt.xlim([-250, 250])
# Adjust Y-limits to match your MATLAB plot roughly
plt.ylim([0, 90]) 
plt.grid(True, alpha=0.5)
plt.tight_layout()
plt.show()


# ===========================================================================
# 6. EXPOSE OUTPUTS FOR UVM SCOREBOARD
# ===========================================================================
# The variables below are the golden-model outputs your UVM scoreboard
# should compare against DUT results.

golden_outputs = {
    # TX
    "chirp_time"       : chirp_time.flatten(),
    "chirp_freq_full"  : chirp_freq_full,
    "stream_pos_in"    : stream_pos_in,
    "stream_neg_in"    : stream_neg_in,
    "X_tx_combined"    : X_tx_combined,
    "x_tx_time"        : x_tx_time,
}

if __name__ == "__main__":
    print("Golden model run complete.")
    print(f"  chirp_time     shape : {golden_outputs['chirp_time'].shape}")
    print(f"  x_tx_time      shape : {golden_outputs['x_tx_time'].shape}")

