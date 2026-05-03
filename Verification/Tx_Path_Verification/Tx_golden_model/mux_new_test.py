"""
    Sponsor: Analog Devices, Inc. (ADI)
    Project: DDS for 5G/6G Cellular Network Calibration and Ranging

    Module: mux_testbench.py

    Description:
        File I/O wrapper for the MUX golden model.
        Reads hardware-level FFT outputs (decimal) and OFDM RAM outputs (hex), 
        slices them to match the RTL memory depths, and passes them through 
        the mux.py function.
"""

import os
import sys
import numpy as np

# Import the core mux function from your mux.py file
from mux_new import mux

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
            
            # Parse hex and apply two's complement
            val = int(line, 16)
            if val & (1 << (bits - 1)):
                val -= (1 << bits)
            vals.append(val)
    return np.array(vals)

def run_mux_testbench(
    fft_re_file="output_real.txt",
    fft_im_file="output_imag.txt",
    ofdm_re_file="ofdm_data_re.hex",
    ofdm_im_file="ofdm_data_im.hex",
    out_re_file="mux_out_real.txt",
    out_im_file="mux_out_imag.txt",
    N=4096
):
    print("--- 1. READING MUX INPUT FILES ---")
    
    # 1A. Read FFT Outputs (Decimal Text Files)
    try:
        fft_re = np.loadtxt(fft_re_file, dtype=np.int64)
        fft_im = np.loadtxt(fft_im_file, dtype=np.int64)
        print(f"  Loaded FFT data: {len(fft_re)} samples")
    except Exception as e:
        print(f"  Error reading FFT files: {e}")
        sys.exit(1)

    # 1B. Read OFDM RAM Outputs (Hexadecimal Files)
    ofdm_re = read_hex_file_to_signed(ofdm_re_file, bits=16)
    ofdm_im = read_hex_file_to_signed(ofdm_im_file, bits=16)
    print(f"  Loaded OFDM RAM data: {len(ofdm_re)} samples")

    # Combine into complex integers (Since the hardware MUX just routes raw wires)
    fft_complex = fft_re + 1j * fft_im
    ofdm_complex = ofdm_re + 1j * ofdm_im

    print("--- 2. FORMATTING STREAMS ---")
    
    # --- STREAM POSITIVE (OFDM) ---
    # The RTL expects up to 2048 OFDM samples. We pass the whole array.
    stream_pos = ofdm_complex

    # --- STREAM NEGATIVE (Radar/Chirp from FFT) ---
    # The RTL radar_ram is strictly 1667 deep. 
    # It grabs the first 1667 bins of the FFT output.
    stream_neg = fft_complex[:1667]

    print("--- 3. RUNNING MUX ---")
    # Feed the arrays into the exact 3-state RTL emulation in mux.py
    X_out, valid_out = mux(stream_pos, stream_neg, N)

    print("--- 4. WRITING MUX OUTPUT FILES ---")
    # Extract the real and imaginary parts back to raw integers
    mux_out_re = np.real(X_out).astype(np.int64)
    mux_out_im = np.imag(X_out).astype(np.int64)

    # Write them out as decimal text files for the IFFT testbench
    np.savetxt(out_re_file, mux_out_re, fmt='%d')
    np.savetxt(out_im_file, mux_out_im, fmt='%d')

    print(f"SUCCESS! MUX outputs written to '{out_re_file}' and '{out_im_file}'.")
    print(f"Total multiplexed frame length: {len(X_out)} samples.")

if __name__ == "__main__":
    # Execute the testbench
    run_mux_testbench()