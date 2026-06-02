"""
    Sponsor: Analog Devices, Inc. (ADI)
    Project: DDS for 5G/6G Cellular Network Calibration and Ranging

    Module: tx_top_golden_model.py

    Description:
        Master wrapper for the full TX datapath.
        Chains DDS -> FFT -> MUX (with OFDM hex files) -> IFFT.
        Generates final TX time-domain outputs for the UVM Scoreboard.
        *Includes comprehensive intermediate stage text file exports for debugging.*
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
    ofdm_re_array=None,  
    ofdm_im_array=None,  
    out_re_file="tx_golden_out_re.txt",
    out_im_file="tx_golden_out_im.txt"
):
    print("==================================================")
    print("  STARTING TX TOP-LEVEL GOLDEN PIPELINE")
    print("==================================================")

    FL_FFT = 8  # Fractional bits for Q8.8

    print(f"[REF] FTW_start={FTW_start}  FTW_step={FTW_step}")
    
    # ---------------------------------------------------------
    # STAGE 1: DDS
    # ---------------------------------------------------------
    print(">> Stage 1: Running Bit-True DDS...")
    chirp_time = generate_dds_golden_model(FTW_start, FTW_step, N_cycles)
    
    # Extract real/imag and force 4096 samples
    x_re = np.real(chirp_time).astype(np.int64).flatten()
    x_im = np.imag(chirp_time).astype(np.int64).flatten()
    
    L = min(len(x_re), len(x_im))
    if L < N:
        x_re = np.pad(x_re, (0, N - L), 'constant')
        x_im = np.pad(x_im, (0, N - L), 'constant')
    elif L > N:
        x_re = x_re[:N]
        x_im = x_im[:N]

    chirp_input = x_re + 1j * x_im
    
    # EXPORT STAGE 1
    np.savetxt("py_fp_stage1_dds_re.txt", x_re, fmt='%d')
    np.savetxt("py_fp_stage1_dds_im.txt", x_im, fmt='%d')

    # =========================================================
    # 1.5 VERIFY DDS AGAINST RTL (Sanity Check)
    # =========================================================
    print(">> Stage 1.5: Comparing Golden DDS against RTL output...")
    try:
        rtl_dds = np.loadtxt("rtl_dds_out.txt", dtype=np.int64)
        golden_dds = x_re
        
        if len(rtl_dds) < N:
            rtl_dds = np.pad(rtl_dds, (0, N - len(rtl_dds)), 'constant')
        elif len(rtl_dds) > N:
            rtl_dds = rtl_dds[:N]
            
        if np.array_equal(golden_dds, rtl_dds):
            print("   [SUCCESS] VERIFICATION PASSED: Python chirp_input matches rtl_dds_out.txt!")
        else:
            print("   [ERROR] VERIFICATION FAILED: Mismatch detected.")
    except Exception as e:
        print(f"   [WARNING] Could not load 'rtl_dds_out.txt' to compare. ({e})")
        
    # ---------------------------------------------------------
    # STAGE 2: FFT 
    # ---------------------------------------------------------
    print(">> Stage 2: Running FFT...")
    chirp_freq = radix2_dif_fft_fixed(chirp_input, WL=16, FL=8)
    
    # EXPORT STAGE 2
    np.savetxt("py_fp_stage2_fft_re.txt", np.real(chirp_freq).astype(np.int64), fmt='%d')
    np.savetxt("py_fp_stage2_fft_im.txt", np.imag(chirp_freq).astype(np.int64), fmt='%d')

    stream_neg_scaled = chirp_freq[:1667]
    
    # ---------------------------------------------------------
    # STAGE 3: OFDM Memory Fetch 
    # ---------------------------------------------------------
    print(">> Stage 3: Fetching OFDM from Scoreboard Arrays...")
    if ofdm_re_array is None or ofdm_im_array is None:
        raise ValueError("Golden model requires ofdm_re_array and ofdm_im_array!")

    ofdm_re_int = np.array(ofdm_re_array, dtype=np.int64)
    ofdm_im_int = np.array(ofdm_im_array, dtype=np.int64)

    ofdm_re_int[1:] = ofdm_re_int[:-1]
    ofdm_im_int[1:] = ofdm_im_int[:-1]
    

    # EXPORT STAGE 3
    np.savetxt('py_fp_stage3_ofdm_re.hex', ofdm_re_int, fmt='%x')    
    np.savetxt('py_fp_stage3_ofdm_im.hex', ofdm_im_int, fmt='%x')
    
    stream_pos_int = ofdm_re_int + 1j * ofdm_im_int

    # ---------------------------------------------------------
    # STAGE 4: MUX
    # ---------------------------------------------------------
    print(">> Stage 4: Running 3-State MUX...")
    X_combined_int, _ = mux(stream_pos_int, stream_neg_scaled, N_cycles)

    # EXPORT STAGE 4
    np.savetxt("py_fp_stage4_mux_re.txt", np.real(X_combined_int).astype(np.int64), fmt='%d')
    np.savetxt("py_fp_stage4_mux_im.txt", np.imag(X_combined_int).astype(np.int64), fmt='%d')

    # ---------------------------------------------------------
    # STAGE 5: IFFT
    # ---------------------------------------------------------
    print(">> Stage 5: Running IFFT...")
    tx_time = radix2_dif_ifft_fixed(X_combined_int, WL=16, FL=5)

    # ---------------------------------------------------------
    # FINAL EXPORT
    # ---------------------------------------------------------
    print(">> Exporting Data for UVM Scoreboard...")
    tx_out_re = np.real(tx_time).astype(np.int64)
    tx_out_im = np.imag(tx_time).astype(np.int64)

    # EXPORT STAGE 5 (FINAL)
    np.savetxt(out_re_file, tx_out_re, fmt='%d')
    np.savetxt(out_im_file, tx_out_im, fmt='%d')

    print(f"\n[SUCCESS] Pipeline Complete! Final Outputs saved to:")
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