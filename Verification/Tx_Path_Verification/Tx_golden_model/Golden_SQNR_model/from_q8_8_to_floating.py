import numpy as np
import os

def convert_hex_q8_8_to_float(hex_re_in, hex_im_in, float_re_out, float_im_out):
    """
    Reads Q8.8 formatted Hex files (16-bit Two's Complement), 
    reverses the fixed-point scaling, and saves them as floating-point text files.
    """
    if not os.path.exists(hex_re_in) or not os.path.exists(hex_im_in):
        print("Error: Could not find the input hex files!")
        return

    print("Reading Q8.8 Hex data...")
    
    # 1. Read the raw hex strings from the files
    with open(hex_re_in, 'r') as f:
        hex_re_lines = f.read().split()
    with open(hex_im_in, 'r') as f:
        hex_im_lines = f.read().split()

    # 2. Convert hex strings to base-10 integers
    ints_re = np.array([int(val, 16) for val in hex_re_lines])
    ints_im = np.array([int(val, 16) for val in hex_im_lines])

    # 3. Apply 16-bit Two's Complement
    # If the value is >= 0x8000 (32768), it is a negative number
    ints_re[ints_re >= 32768] -= 65536
    ints_im[ints_im >= 32768] -= 65536

    # 4. Reverse the Q8.8 Scaling (Divide by 2^8 = 256.0)
    floats_re = ints_re / 256.0
    floats_im = ints_im / 256.0

    # 5. Save the floating-point arrays to text files
    np.savetxt(float_re_out, floats_re, fmt='%.6f')
    np.savetxt(float_im_out, floats_im, fmt='%.6f')

    print(f"Success! Reversed {len(floats_re)} Q8.8 hex samples back to floating-point.")
    print(f"Saved to:\n -> {float_re_out}\n -> {float_im_out}")

if __name__ == "__main__":
    # Your input hex files from the RTL / initial generation
    HEX_RE_IN = "ofdm_data_re.hex"
    HEX_IM_IN = "ofdm_data_im.hex"
    
    # The new floating-point text files ready for your MATLAB/Python MUX
    FLOAT_RE_OUT = "ofdm_float_re.txt"
    FLOAT_IM_OUT = "ofdm_float_im.txt"
    
    convert_hex_q8_8_to_float(HEX_RE_IN, HEX_IM_IN, FLOAT_RE_OUT, FLOAT_IM_OUT)