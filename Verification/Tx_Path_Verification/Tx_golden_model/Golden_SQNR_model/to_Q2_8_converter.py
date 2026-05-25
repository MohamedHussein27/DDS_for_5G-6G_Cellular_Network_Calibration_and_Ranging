import numpy as np
import os

def convert_float_to_q2_8(input_file, output_file):
    """
    Reads a text file of floating point numbers, converts them to Q2.8 
    fixed-point integers, and saves them to a new text file.
    """
    if not os.path.exists(input_file):
        print(f"Error: Could not find {input_file}")
        return

    print(f"Reading floating-point data from {input_file}...")
    
    # 1. Load the floating point vector from the text file
    float_data = np.loadtxt(input_file)
    
    # 2. Convert to Q2.8 format
    # Multiply by 2^8 (256) and round to nearest integer
    q2_8_data = np.round(float_data * 256.0).astype(np.int32)
    
    # Optional: Clamp to 16-bit signed range if your system expects 16-bit words
    # q2_8_data = np.clip(q2_8_data, -32768, 32767)

    # 3. Save the integer vector to the output text file
    np.savetxt(output_file, q2_8_data, fmt='%d')
    print(f"Success! Converted {len(q2_8_data)} samples to Q2.8 integer format.")
    print(f"Saved to {output_file}")

if __name__ == "__main__":
    # Define your filenames here
    INPUT_TXT = "py_sb_dds_float.txt"
    OUTPUT_TXT = "py_sb_dds_q2_8.txt"
    
    convert_float_to_q2_8(INPUT_TXT, OUTPUT_TXT)