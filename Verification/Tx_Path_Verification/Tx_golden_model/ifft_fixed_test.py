"""
    Sponsor: Analog Devices, Inc. (ADI)
    Institution: Faculty of Engineering, Ain Shams University
    Project: DDS for 5G/6G Cellular Network Calibration and Ranging

    Module: ifft_file_test.py

    Description:
        4096-POINT IFFT TEST BENCH (FILE INPUT)
        Reads integer sequences from text files, runs them through the 
        bit-true Python model, and writes the results out to text files.
"""

import os
import sys
import numpy as np

# Import the IFFT wrapper
from ifft_fixed import radix2_dif_ifft_fixed
from fft_fixed  import radix2_dif_fft_fixed

# =========================================================
# CONFIGURATION
# =========================================================
N  = 4096
WL = 16
FL = 5  # Note: I changed this to 5 to match your earlier scripts. 
        # Your MATLAB prompt said FL=8, adjust if needed!

def main():
    # Get the directory where this script is located
    script_dir = os.path.dirname(os.path.abspath(__file__))
    if not script_dir:
        script_dir = os.getcwd()

    # =========================================================
    # 1. READ INPUT FILES
    # =========================================================
    file_re = os.path.join(script_dir, 'rtl_dds_out.txt')
    file_im = os.path.join(script_dir, 'input_imag_test.txt')

    if not os.path.exists(file_re) or not os.path.exists(file_im):
        print(f"Error: Input files not found in {script_dir}", file=sys.stderr)
        sys.exit(1)

    try:
        x_re = np.loadtxt(file_re, dtype=np.int64)
    except Exception as e:
        x_re = np.array([], dtype=np.int64)
        
    try:
        x_im = np.loadtxt(file_im, dtype=np.int64)
    except Exception as e:
        x_im = np.array([], dtype=np.int64)

    # Flatten arrays just in case they were read as 2D
    x_re = x_re.flatten()
    x_im = x_im.flatten()

    # =========================================================
    # 2. FORCE EXACT 4096 SAMPLES (VERY IMPORTANT)
    # =========================================================
    L = min(len(x_re), len(x_im))
    
    if L < N:
        print(f"Padding input: {L} → {N} samples")
        # Pad with zeros up to length N
        x_re = np.pad(x_re, (0, N - L), 'constant')
        x_im = np.pad(x_im, (0, N - L), 'constant')
    else:
        print(f"Truncating input to {N} samples")
        x_re = x_re[:N]
        x_im = x_im[:N]

    x_input = x_re + 1j * x_im
    print(f"x_input:", x_input[:100])

    # =========================================================
    # 3. RUN YOUR BIT-TRUE FFT/IFFT
    # =========================================================
    print(f"Running radix2_dif_ifft_fixed (N={N}, WL={WL}, FL={FL})...")
    X_out = radix2_dif_ifft_fixed(x_input, WL, FL, False)
    print(f"X_out:", X_out[:100])

    # =========================================================
    # 4. WRITE OUTPUT FILES
    # =========================================================
    out_re = np.real(X_out).astype(np.int64)
    out_im = np.imag(X_out).astype(np.int64)

    out_file_re = os.path.join(script_dir, 'out_fft_real.txt')
    out_file_im = os.path.join(script_dir, 'out_fft_imag.txt')

    np.savetxt(out_file_re, out_re, fmt='%d')
    np.savetxt(out_file_im, out_im, fmt='%d')

    # =========================================================
    # DONE
    # =========================================================
    print("\nFFT DONE")
    print(f"Outputs written to:\n  {out_file_re}\n  {out_file_im}")

if __name__ == "__main__":
    main()