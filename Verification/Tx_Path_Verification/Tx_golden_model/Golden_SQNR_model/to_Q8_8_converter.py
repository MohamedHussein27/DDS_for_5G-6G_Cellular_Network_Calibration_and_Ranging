import numpy as np
import os

def convert_complex_float_to_q8_8(re_in_file, im_in_file, re_out_file, im_out_file):
    """
    Reads floating-point Real and Imaginary text files, converts them 
    to Q8.8 fixed-point notation, and saves them to new text files.
    """
    if not os.path.exists(re_in_file) or not os.path.exists(im_in_file):
        print("Error: Could not find the input files!")
        return

    print("Reading floating-point FFT data...")
    float_re = np.loadtxt(re_in_file)
    float_im = np.loadtxt(im_in_file)
    
    # 1. Convert to Q8.8 format (Multiply by 2^8 = 256)
    q8_8_re = np.round(float_re * 256.0).astype(np.int64)
    q8_8_im = np.round(float_im * 256.0).astype(np.int64)
    
    # 2. Enforce 16-bit Signed Integer Boundaries [-32768 to 32767]
    # This acts as your hardware saturation logic!
    q8_8_re = np.clip(q8_8_re, -32768, 32767)
    q8_8_im = np.clip(q8_8_im, -32768, 32767)

    # 3. Save the integer vectors to output text files
    np.savetxt(re_out_file, q8_8_re, fmt='%d')
    np.savetxt(im_out_file, q8_8_im, fmt='%d')
    
    print(f"Success! Converted {len(q8_8_re)} samples to Q8.8 integer format.")
    print(f"Saved to:\n -> {re_out_file}\n -> {im_out_file}")

if __name__ == "__main__":
    # Define your filenames here
    RE_IN  = "py_sb_mux_re_float.txt"
    IM_IN  = "py_sb_mux_im_float.txt"
    
    RE_OUT = "py_sb_mux_re_q8_8.txt"
    IM_OUT = "py_sb_mux_im_q8_8.txt"
    
    convert_complex_float_to_q8_8(RE_IN, IM_IN, RE_OUT, IM_OUT)