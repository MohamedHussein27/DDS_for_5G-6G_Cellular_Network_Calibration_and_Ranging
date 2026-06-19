import numpy as np

# ─────────────────────────────────────────────────────────────────────────────
# Constants from mixed_ofdm_radar_optimal_seq
# ─────────────────────────────────────────────────────────────────────────────
WL           = 16
FL           = 8
N_FFT        = 4096
N_HALF       = 2048
N_OFDM_SC    = 1500
MAX_VAL      =  (1 << (WL - 1)) - 1
MIN_VAL      = -(1 << (WL - 1))

# =============================================================================
# THE "GOLDEN RATIO" CONSTANTS
# =============================================================================
TARGET_MAX_RX_INT = 500
AMP_REF_FLOAT     = 10.0

def _f2q(f: float) -> int:
    return max(MIN_VAL, min(MAX_VAL, int(round(f * (1 << FL)))))

def _build_mixed_rx_frame_optimal(seed: int = 0) -> list:
    rng  = np.random.default_rng(seed)
    freq = np.zeros(N_FFT, dtype=np.complex128)

    sc_start   = (N_HALF - N_OFDM_SC) // 2
    sc_end     = sc_start + N_OFDM_SC
    qam_levels = np.array([-7, -5, -3, -1, 1, 3, 5, 7], dtype=np.float64)
    
    freq[sc_start:sc_end] = (1.0 / 7.0) * (
        rng.choice(qam_levels, size=N_OFDM_SC) +
        1j * rng.choice(qam_levels, size=N_OFDM_SC)
    )

    for k in range(N_HALF, N_FFT):
        k_local = k - N_HALF
        phase   = np.pi * k_local * k_local / N_HALF
        freq[k] = 1.0 * np.exp(1j * phase)

    td = np.fft.ifft(freq) * N_FFT
    
    max_td_val = np.max(np.abs(td))
    scale = (TARGET_MAX_RX_INT / float(1 << FL)) / max_td_val
    td_scaled = td * scale

    return [(_f2q(float(s.real)), _f2q(float(s.imag))) for s in td_scaled]

def _build_chirp_reference_optimal(seed: int = 0) -> list:
    ref = []
    for k in range(N_HALF):
        phase = np.pi * k * k / N_HALF
        c     = AMP_REF_FLOAT * np.exp(1j * phase)
        ref.append((_f2q(float(c.real)), _f2q(float(c.imag))))
    return ref

# ─────────────────────────────────────────────────────────────────────────────
# File Generation Execution
# ─────────────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    print("Generating Optimal 'Golden Ratio' Sequence Files...")
    
    # Generate data using the exact same seed (0) as the UVM testbench
    rx_samples  = _build_mixed_rx_frame_optimal(seed=0)
    ref_samples = _build_chirp_reference_optimal(seed=0)
    
    # 1. Write the 4096 RX samples
    with open("rtl_rx_in_re.txt", "w") as fre, open("rtl_rx_in_im.txt", "w") as fim:
        for re, im in rx_samples:
            fre.write(f"{re}\n")
            fim.write(f"{im}\n")
    print(" -> Wrote 4096 samples to rtl_rx_in_re.txt and rtl_rx_in_im.txt")

    # 2. Write the 2048 Reference samples
    with open("rtl_mux_out_re.txt", "w") as fre, open("rtl_mux_out_im.txt", "w") as fim:
        for re, im in ref_samples:
            fre.write(f"{re}\n")
            fim.write(f"{im}\n")
    print(" -> Wrote 2048 samples to rtl_mux_out_re.txt and rtl_mux_out_im.txt")
    
    print("\nSuccess! Files are ready. You can now run your Verilog RX_TOP_TB in Questasim.")