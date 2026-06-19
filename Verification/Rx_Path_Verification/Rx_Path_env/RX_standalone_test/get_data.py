import numpy as np

# Constants
WL = 16
FL = 8
N_FFT = 4096
N_HALF = 2048
N_OFDM_SC = 1500
MAX_VAL = (1 << (WL - 1)) - 1
MIN_VAL = -(1 << (WL - 1))

# Scaling (matching your required 225 peak)
TARGET_MAX_RX_INT = 225
AMP_REF_FLOAT = 120.0

def _f2q(f: float) -> int:
    return max(MIN_VAL, min(MAX_VAL, int(round(f * (1 << FL)))))

def main():
    print("Generating exact Seed 0 arrays...")
    
    # 1. Generate rx_in (OFDM + Chirp)
    rng = np.random.default_rng(0) # Exact Seed 0
    freq = np.zeros(N_FFT, dtype=np.complex128)

    sc_start = (N_HALF - N_OFDM_SC) // 2
    sc_end = sc_start + N_OFDM_SC
    qam_levels = np.array([-7, -5, -3, -1, 1, 3, 5, 7], dtype=np.float64)
    
    freq[sc_start:sc_end] = (1.0 / 7.0) * (
        rng.choice(qam_levels, size=N_OFDM_SC) +
        1j * rng.choice(qam_levels, size=N_OFDM_SC)
    )

    for k in range(N_HALF, N_FFT):
        k_local = k - N_HALF
        phase = np.pi * k_local * k_local / N_HALF
        freq[k] = 1.0 * np.exp(1j * phase)

    td = np.fft.ifft(freq) * N_FFT
    scale = (TARGET_MAX_RX_INT / float(1 << FL)) / np.max(np.abs(td))
    td_scaled = td * scale

    rx_re = np.array([_f2q(float(s.real)) for s in td_scaled])
    rx_im = np.array([_f2q(float(s.imag)) for s in td_scaled])

    # 2. Generate ref_ram (Chirp Reference)
    ref_re = np.zeros(N_HALF, dtype=int)
    ref_im = np.zeros(N_HALF, dtype=int)
    for k in range(N_HALF):
        phase = np.pi * k * k / N_HALF
        c = AMP_REF_FLOAT * np.exp(1j * phase)
        ref_re[k] = _f2q(float(c.real))
        ref_im[k] = _f2q(float(c.imag))

    # 3. Export to text files
    print("Writing files...")
    np.savetxt("matlab_rx_in_re_seed0.txt", rx_re, fmt="%d")
    np.savetxt("matlab_rx_in_im_seed0.txt", rx_im, fmt="%d")
    np.savetxt("matlab_ref_in_re_seed0.txt", ref_re, fmt="%d")
    np.savetxt("matlab_ref_in_im_seed0.txt", ref_im, fmt="%d")

    print("\nDONE! You can now load these 4 files into MATLAB.")

if __name__ == "__main__":
    main()