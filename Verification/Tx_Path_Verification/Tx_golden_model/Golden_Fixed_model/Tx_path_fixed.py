"""
system_builtin_function.py
==========================
Python golden model — perfectly aligned with SYSTEM_builtin.m
Refactored as a callable function.
"""

import numpy as np
import os
from mytypes import quantize_fixed, reinterpretcast
from fft_fixed import radix2_dif_fft_fixed
from ifft_fixed import radix2_dif_ifft_fixed
from mux import mux

def generate_tx_golden_model(FTW_start=0, FTW_step=426666, N_cycles=4096, 
                             ofdm_hex_re_file="ofdm_data_re.hex", 
                             ofdm_hex_im_file="ofdm_data_im.hex", 
                             output_filename_re="matlab_tx_out_re.txt",
                             output_filename_im="matlab_tx_out_im.txt"):
    """
    Generates the golden TX path outputs exactly matching the MATLAB script.
    
    Args:
        FTW_start: Starting Frequency Tuning Word
        FTW_step: Frequency Tuning Word step per clock cycle
        N_cycles: Number of samples / FFT size (e.g., 4096)
        ofdm_hex_re_file: Text file containing real OFDM symbols in hex
        ofdm_hex_im_file: Text file containing imaginary OFDM symbols in hex
        output_filename_re: Name of the text file to save the IFFT Real results
        output_filename_im: Name of the text file to save the IFFT Imag results
    
    Returns:
        ifft_out_real, ifft_out_imag (as raw integers)
    """
    # ===========================================================================
    # 0. WORD-LENGTH PARAMETERS (Matches MATLAB exactly)
    # ===========================================================================
    WL_reduction = 0
    scale        = 128          # CORRECTED: MATLAB uses 128

    WL_FFT  = 16 - WL_reduction;   FL_FFT  =  8 - WL_reduction
    WL_FFT1 = 16 - WL_reduction;   FL_FFT1 =  8 - WL_reduction
    WL_FFT2 = 16 - WL_reduction;   FL_FFT2 = 14 - WL_reduction
    WL_IFFT2= 16 - WL_reduction;   FL_IFFT2=  5 - WL_reduction

    # ===========================================================================
    # 1. DESIGN PARAMETERS
    # ===========================================================================
    Fs       = 491.52e6        # System clock (Hz)
    N        = int(N_cycles)   # 4096
    
    ADDRESS_WIDTH = 16
    MEMORY_WIDTH  = 8

    # ===========================================================================
    # 4. TX PATH
    # ===========================================================================
    print("--- TX PATH GENERATION ---")

    # --- 4.1 DDS Chirp (Exact Bit-True Hardware Logic from MATLAB) ---
    lut_depth = 2 ** (ADDRESS_WIDTH - 2)
    n_lut = np.arange(lut_depth)
    lut_table = np.round((2**(MEMORY_WIDTH-1)-1) * np.sin((np.pi/2) * n_lut / lut_depth)).astype(np.int32)

    # Vectorized Phase Accumulation
    cycles_arr = np.arange(1, N + 1, dtype=np.uint64)
    accumulation_term = (cycles_arr * (cycles_arr - np.uint64(1))) // np.uint64(2)
    phase_acc = (np.uint64(FTW_start) * cycles_arr + np.uint64(FTW_step) * accumulation_term)
    
    # 32-bit Truncation and Quadrant extraction
    phase_32 = (phase_acc & 0xFFFFFFFF).astype(np.uint32)
    phase_16 = phase_32 >> (32 - ADDRESS_WIDTH)
    quadrant = phase_16 >> 14

    mapped_addr = np.zeros_like(phase_16, dtype=np.uint32)
    neg_flag = np.zeros_like(phase_16, dtype=bool)

    # 1st Quad
    q0 = (quadrant == 0)
    mapped_addr[q0] = phase_16[q0]
    neg_flag[q0] = False

    # 2nd Quad
    q1 = (quadrant == 1)
    mapped_addr[q1] = 32768 - 1 - phase_16[q1]
    neg_flag[q1] = False

    # 3rd Quad
    q2 = (quadrant == 2)
    mapped_addr[q2] = phase_16[q2] - 32768
    neg_flag[q2] = True

    # 4th Quad
    q3 = (quadrant == 3)
    mapped_addr[q3] = 65536 - 1 - phase_16[q3]
    neg_flag[q3] = True

    amp = lut_table[mapped_addr]
    chirp_time = np.where(neg_flag, -amp, amp).astype(np.int64)

    # --- 4.2 FFT of chirp ---  
    # Because chirp_time is already integers, it implicitly acts as the fi(..., 16, 8) 
    # fractional representation when passed into our integer-based Python FFT!
    chirp_freq_full = radix2_dif_fft_fixed(chirp_time, WL_FFT1, FL_FFT1)

    # Reinterpret cast (scale=128 -> log2_scale=7)
    log2_scale = int(round(np.log2(scale)))
    chirp_freq_full = reinterpretcast(
        chirp_freq_full,
        src_WL=WL_FFT1, src_FL=FL_FFT1,
        dst_WL=WL_FFT,  dst_FL=FL_FFT + log2_scale
    )
    chirp_freq_full = chirp_freq_full.flatten()

    # --- 4.3 Read OFDM symbols from Hex files ---
    def read_hex_file_to_signed(filename, bits=16):
        """Reads a file of hex strings and converts them to signed integers."""
        if not os.path.exists(filename):
            raise FileNotFoundError(f"Missing OFDM file: {filename}")
            
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

    print(f"Reading OFDM data from {ofdm_hex_re_file} and {ofdm_hex_im_file}...")
    ofdm_real_int = read_hex_file_to_signed(ofdm_hex_re_file, bits=16)
    ofdm_imag_int = read_hex_file_to_signed(ofdm_hex_im_file, bits=16)

    # Convert to complex floats for the MUX
    ofdm_symbols = (ofdm_real_int + 1j * ofdm_imag_int) / (1 << FL_FFT)

    # --- 4.4 Fixed-point cast for OFDM stream ---
    def fi_cast(x, WL, FL):
        """Convergent rounding + Saturate – matches MATLAB fimath."""
        return quantize_fixed(x, WL=WL, FL=FL, signed=True,
                              rounding='convergent', overflow='saturate')

    # --- 4.5 Build stream_pos (OFDM, positive frequencies) ---
    target_len    = N // 2   # 2048
    stream_pos_in = ofdm_symbols.copy()

    if len(stream_pos_in) < target_len:
        pad_len = target_len - len(stream_pos_in)
        # 83 bins of zero-padding at DC (Matches MATLAB)
        stream_pos_in = np.concatenate([
            np.zeros(83, dtype=complex),
            stream_pos_in,
            np.zeros(pad_len - 83, dtype=complex)
        ])

    # --- 4.6 Build stream_neg (Chirp, negative frequencies) ---
    # Fs = 491.52e6, N = 4096. df = Fs/N
    # 200e6 / df = 1666.66 -> round to 1667
    idx_cutoff = 1667 

    chirp_freq_half = chirp_freq_full[:N // 2].copy()
    chirp_freq_half[idx_cutoff:] = 0               # zero out >200 MHz bins

    stream_neg_in = chirp_freq_half[:1667].copy()  # MATLAB: stream_neg_in(1668:end) = []

    # PAD RADAR: Pad at the START (Matches MATLAB)
    if len(stream_neg_in) < target_len:
        pad_len = target_len - len(stream_neg_in)
        stream_neg_in = np.concatenate([
            np.zeros(pad_len, dtype=complex),
            stream_neg_in
        ])

    # --- 4.7 MUX ---
    X_tx_combined, tx_valid = mux(stream_pos_in, stream_neg_in, N)
    X_tx_combined = fi_cast(X_tx_combined, WL_FFT, FL_FFT)

    # --- 4.8 IFFT ---
    x_tx_time = radix2_dif_ifft_fixed(X_tx_combined, WL_FFT, FL_FFT)

    print(f"  TX time-domain samples generated: {len(x_tx_time)}")

    # ===========================================================================
    # 5. EXTRACT & FILE I/O
    # ===========================================================================
    # Extract Real and Imaginary components (Matches MATLAB tx_re_int / tx_im_int)
    # Scale by 2^FL_FFT to bring floats back to exact RTL raw integers
    ifft_out_real = np.round(np.real(x_tx_time) * (2 ** FL_FFT)).astype(np.int64)
    ifft_out_imag = np.round(np.imag(x_tx_time) * (2 ** FL_FFT)).astype(np.int64)

    # Write output to the requested text files
    print(f"--- WRITING OUTPUT TO {output_filename_re} & {output_filename_im} ---")
    
    with open(output_filename_re, 'w') as f_re:
        for r in ifft_out_real:
            f_re.write(f"{r}\n")
            
    with open(output_filename_im, 'w') as f_im:
        for i in ifft_out_imag:
            f_im.write(f"{i}\n")
            
    print("Golden model generation and file write complete.")

    return ifft_out_real, ifft_out_imag

# ===========================================================================
# Example Execution
# ===========================================================================
if __name__ == "__main__":
    
    # Run the function (Make sure your hex files exist in the same directory!)
    out_re, out_im = generate_tx_golden_model(
        FTW_start=0, 
        FTW_step=426666, 
        N_cycles=4096, 
        ofdm_hex_re_file="ofdm_data_re.hex",  
        ofdm_hex_im_file="ofdm_data_im.hex", 
        output_filename_re="matlab_tx_out_re.txt",
        output_filename_im="matlab_tx_out_im.txt"
    )