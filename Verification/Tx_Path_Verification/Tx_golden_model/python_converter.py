import re
import os

def convert_coe_to_txt(coe_filepath, txt_filepath, bit_width=16, is_signed=True):
    """
    Reads a .coe file and converts its data to decimal, saving to a .txt file.
    
    Parameters:
    - coe_filepath: Path to the input .coe file.
    - txt_filepath: Path to the output .txt file.
    - bit_width: The bit width of the data (needed for signed two's complement).
    - is_signed: Boolean indicating if the data should be treated as signed.
    """
    if not os.path.exists(coe_filepath):
        print(f"Error: File '{coe_filepath}' not found.")
        return

    try:
        with open(coe_filepath, 'r') as file:
            content = file.read()

        # 1. Extract the radix (base) - usually 2, 10, or 16
        radix_match = re.search(r'memory_initialization_radix\s*=\s*(\d+)\s*;', content, re.IGNORECASE)
        radix = int(radix_match.group(1)) if radix_match else 16 # Default to hex if missing

        # 2. Extract the data vector
        vector_match = re.search(r'memory_initialization_vector\s*=\s*([^;]+);', content, re.IGNORECASE)
        if not vector_match:
            print(f"Error: Could not find 'memory_initialization_vector' in {coe_filepath}")
            return

        vector_data = vector_match.group(1)
        
        # Split the data by commas, spaces, or newlines
        raw_values = re.split(r'[\s,]+', vector_data.strip())

        decimal_values = []
        for val in raw_values:
            if not val: 
                continue
                
            # Convert string to integer based on the radix
            num = int(val, radix)

            # 3. Handle signed data (Two's complement)
            if is_signed:
                # If the sign bit is 1, subtract the full scale to get the negative value
                if num >= (1 << (bit_width - 1)):
                    num -= (1 << bit_width)

            decimal_values.append(num)

        # 4. Write to the output text file
        with open(txt_filepath, 'w') as file:
            for num in decimal_values:
                file.write(f"{num}\n")

        print(f"Success: Converted '{coe_filepath}' to '{txt_filepath}' ({len(decimal_values)} values).")

    except Exception as e:
        print(f"An error occurred while processing '{coe_filepath}': {e}")


# --- Execution ---
if __name__ == "__main__":
    # IMPORTANT: Set your bit width here. 
    # Usually, DAC/ADC data is 12, 14, or 16 bits.
    DATA_BIT_WIDTH = 16 
    TREAT_AS_SIGNED = True

    # Process Real Part
    convert_coe_to_txt(
        coe_filepath='tx_final_out_re.coe',
        txt_filepath='tx_final_out_re.txt',
        bit_width=DATA_BIT_WIDTH,
        is_signed=TREAT_AS_SIGNED
    )

    # Process Imaginary Part
    convert_coe_to_txt(
        coe_filepath='tx_final_out_im.coe',
        txt_filepath='tx_final_out_im.txt',
        bit_width=DATA_BIT_WIDTH,
        is_signed=TREAT_AS_SIGNED
    )