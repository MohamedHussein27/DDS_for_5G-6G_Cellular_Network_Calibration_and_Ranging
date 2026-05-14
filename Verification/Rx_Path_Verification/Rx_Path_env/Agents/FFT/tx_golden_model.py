"""
    Sponsor: Analog Devices, Inc. (ADI)
    Project: DDS for 5G/6G Cellular Network Calibration and Ranging

    Module: tx_top_golden_model.py

    Description:
        Master wrapper for the full TX datapath.
        Chains DDS -> FFT -> MUX (with OFDM hex files) -> IFFT.
        Generates final TX time-domain outputs for the UVM Scoreboard.
"""

import os
import sys
import numpy as np

# --- Import all verified sub-blocks ---
from dds_new_model import generate_dds_golden_model
from fft_fixed import radix2_dif_fft_fixed
from ifft_fixed import radix2_dif_ifft_fixed
from mux_new import mux

N  = 4096
# ===========================================================================
# HELPER FUNCTIONS
# ===========================================================================
def read_hex_file_to_signed(filename, bits=16):
    """Reads a file of hex strings and converts them to signed integers."""
    if not os.path.exists(filename):
        print(f"Error: Could not find {filename}", file=sys.stderr)
        sys.exit(1)
        
    vals = []
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line: continue
            val = int(line, 16)
            if val & (1 << (bits - 1)):
                val -= (1 << bits)
            vals.append(val)
    return np.array(vals)

# ===========================================================================
# MASTER PIPELINE FUNCTION
# ===========================================================================
def run_tx_top_pipeline(
    FTW_start=0, 
    FTW_step=426666, 
    N_cycles=4096, 
    Fs=491.52e6,
    ofdm_re_array=None,  # NEW: Expecting a Python list/numpy array
    ofdm_im_array=None,  # NEW: Expecting a Python list/numpy array
    out_re_file="tx_golden_out_re.txt",
    out_im_file="tx_golden_out_im.txt"
):
    print("==================================================")
    print("  STARTING TX TOP-LEVEL GOLDEN PIPELINE")
    print("==================================================")

    FL_FFT = 8  # Fractional bits for Q8.8


    # ---------------------------------------------------------
    # STAGE 1: DDS
    # ---------------------------------------------------------
    print(">> Stage 1: Running Bit-True DDS...")
    # chirp_time is an array of raw integers (e.g., ±127)
    chirp_time = generate_dds_golden_model(FTW_start, FTW_step, N_cycles)
    print(f"DDS first outputs:", chirp_time[4090:])  

    # =========================================================
    # 1. EXTRACT REAL AND IMAGINARY PARTS
    # =========================================================
    # Safely grab the real and imaginary parts of the chirp
    x_re = np.real(chirp_time).astype(np.int64)
    x_im = np.imag(chirp_time).astype(np.int64)

    # Flatten arrays just in case they were generated as 2D (e.g., column vectors)
    x_re = x_re.flatten()
    x_im = x_im.flatten()

    # =========================================================
    # 2. FORCE EXACT 4096 SAMPLES (VERY IMPORTANT)
    # =========================================================
    N = 4096
    L = min(len(x_re), len(x_im))
    
    if L < N:
        print(f"Padding chirp input: {L} → {N} samples")
        # Pad with zeros up to length N
        x_re = np.pad(x_re, (0, N - L), 'constant')
        x_im = np.pad(x_im, (0, N - L), 'constant')
    elif L > N:
        print(f"Truncating chirp input to {N} samples")
        x_re = x_re[:N]
        x_im = x_im[:N]

    # Recombine into a strictly sized complex array
    chirp_input = x_re + 1j * x_im

    print (f"chirp before fft:", chirp_input[:100])
    # =========================================================
    # 2.5 VERIFY DDS AGAINST RTL (Sanity Check)
    # =========================================================
    print(">> Stage 1.5: Comparing Golden DDS against RTL output...")
    try:
        rtl_dds = np.loadtxt("rtl_dds_out.txt", dtype=np.int64)
        golden_dds = np.real(chirp_input).astype(np.int64)
        
        # Ensure lengths match before comparing (Truncate or pad RTL if necessary)
        if len(rtl_dds) < N:
            rtl_dds = np.pad(rtl_dds, (0, N - len(rtl_dds)), 'constant')
        elif len(rtl_dds) > N:
            rtl_dds = rtl_dds[:N]
            
        # Compare them
        if np.array_equal(golden_dds, rtl_dds):
            print("   [SUCCESS] VERIFICATION PASSED: Python chirp_input perfectly matches rtl_dds_out.txt!")
        else:
            print("   [ERROR] VERIFICATION FAILED: Mismatch detected.")
            # Find exactly where it failed to help with debugging
            mismatches = np.where(golden_dds != rtl_dds)[0]
            first_mismatch = mismatches[0]
            print(f"      -> Total mismatched samples: {len(mismatches)} out of {N}")
            print(f"      -> First mismatch at index {first_mismatch}: Python = {golden_dds[first_mismatch]}, RTL = {rtl_dds[first_mismatch]}")
            
    except Exception as e:
        print(f"   [WARNING] Could not load 'rtl_dds_out.txt' to compare. ({e})")
    # ---------------------------------------------------------
    # STAGE 2: FFT + Reinterpret Cast
    # ---------------------------------------------------------
    print(">> Stage 2: Running FFT...")
    # FFT output is fractional float
    chirp_freq = radix2_dif_fft_fixed(chirp_input, WL=16, FL=8)

    print(f"chirp first values =", chirp_freq[:100])
       
    # Apply the scale=128 cast (log2_scale = 7)
    #chirp_freq_cast = reinterpretcast(chirp_freq, src_WL=16, src_FL=8, dst_WL=16, dst_FL=8 + 7).flatten()

    # FIX: Keep as complex floats. The MUX will safely cast the real/imag components 
    # to integers individually before applying the >> 7 shift!
    #stream_neg_scaled = np.round(chirp_freq_cast * (2**FL_FFT))

    stream_neg_scaled = chirp_freq[:1667]
    # ---------------------------------------------------------
    # STAGE 3: OFDM Memory Fetch (UPDATED)
    # ---------------------------------------------------------
    print(">> Stage 3: Fetching OFDM from Scoreboard Arrays...")
    
    # Check if arrays were passed, else create empty ones to prevent crashing
    if ofdm_re_array is None or ofdm_im_array is None:
        raise ValueError("Golden model requires ofdm_re_array and ofdm_im_array!")

    ofdm_re_int = np.array(ofdm_re_array, dtype=np.int64)
    ofdm_im_int = np.array(ofdm_im_array, dtype=np.int64)

    # Shift values: fill index 1 to the end with values from index 0 to second-to-last
    ofdm_re_int[1:] = ofdm_re_int[:-1]
    ofdm_im_int[1:] = ofdm_im_int[:-1]

    np.savetxt('ofdm_re_from_mem.hex', ofdm_re_int, fmt='%x')    
    np.savetxt('ofdm_im_from_mem.hex', ofdm_im_int, fmt='%x')

    print(f"OFDM real first values =", ofdm_re_int[:100])
    print(f"OFDM imag first values =", ofdm_im_int[:100])
    
    # Combine into complex integers
    stream_pos_int = ofdm_re_int + 1j * ofdm_im_int

    print(f"OFDM combined first values =", stream_pos_int[:100])


    # ---------------------------------------------------------
    # STAGE 4: MUX
    # ---------------------------------------------------------
    print(">> Stage 4: Running 3-State MUX...")
    # The MUX will grab the first 1667 items of stream_neg, shift them >> 7, and pack the frame.
    X_combined_int, _ = mux(stream_pos_int, stream_neg_scaled, N_cycles)

    # Convert back to fractional floats for the IFFT engine

    # ---------------------------------------------------------
    # STAGE 5: IFFT
    # ---------------------------------------------------------
    print(">> Stage 5: Running IFFT...")
    tx_time = radix2_dif_ifft_fixed(X_combined_int, WL=16, FL=5)


    # ---------------------------------------------------------
    # FINAL EXPORT
    # ---------------------------------------------------------
    print(">> Exporting Data for UVM Scoreboard...")
    # Convert final floats back to 16-bit RTL integers for scoreboard comparison
    #tx_out_re = np.round(np.real(tx_time) * (2**FL_FFT)).astype(np.int64)
    #tx_out_im = np.round(np.imag(tx_time) * (2**FL_FFT)).astype(np.int64)

    tx_out_re = tx_time.real
    tx_out_im = tx_time.imag

    np.savetxt(out_re_file, tx_out_re, fmt='%d')
    np.savetxt(out_im_file, tx_out_im, fmt='%d')

    print(f"\n[SUCCESS] Pipeline Complete! Outputs saved to:")
    print(f"   - {out_re_file}")
    print(f"   - {out_im_file}")

    return tx_out_re, tx_out_im

# ===========================================================================
# EXECUTION
# ===========================================================================
if __name__ == "__main__":
    try:
        re_out, im_out = run_tx_top_pipeline()
    except Exception as e:
        print(f"\n[ERROR] Pipeline Failed: {e}")