import itertools
import os

def compare_txt_files(file1_path, file2_path):
    """
    Compares two text files line by line and returns the indices of mismatches.
    
    Returns a list of dictionaries containing the index (0-based) and the 
    differing values from each file.
    """
    if not os.path.exists(file1_path) or not os.path.exists(file2_path):
        print("Error: One or both of the files do not exist.")
        return []

    mismatched_data = []

    with open(file1_path, 'r') as f1, open(file2_path, 'r') as f2:
        # zip_longest pairs up lines, filling with None if one file is shorter
        for index, (line1, line2) in enumerate(itertools.zip_longest(f1, f2)):
            
            # Clean up whitespace/newlines. If a file ended, it remains None.
            val1 = line1.strip() if line1 is not None else "<EOF - End of File>"
            val2 = line2.strip() if line2 is not None else "<EOF - End of File>"
            
            # Compare the values
            if val1 != val2:
                mismatched_data.append({
                    'index': index,      # 0-based index
                    'line_num': index + 1, # 1-based line number for easier reading
                    'file1_val': val1,
                    'file2_val': val2
                })

    return mismatched_data

# --- Execution ---
if __name__ == "__main__":
    file_A = 'matlab_mux_im_q8_8.txt'
    file_B = 'py_sb_mux_im_q8_8.txt'
    
    # Create some dummy files for testing if they don't exist
    if not os.path.exists(file_A):
        with open(file_A, 'w') as f: f.write("10\n20\n30\n40\n")
        with open(file_B, 'w') as f: f.write("10\n25\n30\n")

    print(f"Comparing '{file_A}' and '{file_B}'...\n")
    
    mismatches = compare_txt_files(file_A, file_B)

    if not mismatches:
        print("Success: The files match perfectly!")
    else:
        print(f"Found {len(mismatches)} mismatches:")
        print("-" * 50)
        for mismatch in mismatches:
            print(f"Index: {mismatch['index']} (Line {mismatch['line_num']})")
            print(f"  File 1: {mismatch['file1_val']}")
            print(f"  File 2: {mismatch['file2_val']}")
            print("-" * 50)