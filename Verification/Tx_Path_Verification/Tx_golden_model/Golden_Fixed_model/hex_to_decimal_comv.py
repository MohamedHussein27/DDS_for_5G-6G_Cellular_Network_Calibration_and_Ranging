import os

def convert_hex_to_txt(hex_filepath, txt_filepath, bit_width=16, is_signed=True):
    """
    Reads a raw .hex file, converts the signed hex data to decimal, 
    and saves it to a .txt file.
    """
    if not os.path.exists(hex_filepath):
        print(f"Error: File '{hex_filepath}' not found.")
        return

    decimal_values = []
    
    try:
        with open(hex_filepath, 'r') as file:
            for line in file:
                clean_line = line.strip()
                raw_values = clean_line.split()
                
                for val in raw_values:
                    if not val:
                        continue
                        
                    # 1. Convert hex string to integer (Base 16)
                    num = int(val, 16)

                    # 2. Handle signed data (Two's complement)
                    if is_signed:
                        if num >= (1 << (bit_width - 1)):
                            num -= (1 << bit_width)

                    decimal_values.append(num)

        # 3. Write to the output text file
        with open(txt_filepath, 'w') as file:
            for num in decimal_values:
                file.write(f"{num}\n")

        print(f"Success: Converted '{hex_filepath}' to '{txt_filepath}' ({len(decimal_values)} values).")

    except ValueError as e:
        print(f"Error: Found invalid hex data in '{hex_filepath}'. Details: {e}")
    except Exception as e:
        print(f"An error occurred while processing '{hex_filepath}': {e}")


# --- Execution ---
if __name__ == "__main__":
    # IMPORTANT: Adjust this if your RTL output is not 16-bit!
    DATA_BIT_WIDTH = 16 
    TREAT_AS_SIGNED = True

    print("Starting conversion...")

    # Process Real Part
    convert_hex_to_txt(
        hex_filepath='rtl_mux_out_re.hex',
        txt_filepath='rtl_mux_out_re.txt',
        bit_width=DATA_BIT_WIDTH,
        is_signed=TREAT_AS_SIGNED
    )

    # Process Imaginary Part
    convert_hex_to_txt(
        hex_filepath='rtl_mux_out_im.hex',
        txt_filepath='rtl_mux_out_im.txt',
        bit_width=DATA_BIT_WIDTH,
        is_signed=TREAT_AS_SIGNED
    )
    
    print("Done.")