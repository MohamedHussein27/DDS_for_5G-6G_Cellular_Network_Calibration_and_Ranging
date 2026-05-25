import numpy as np
import os

def calculate_sqnr(ref_vector, dut_vector):
    """
    Calculates the Signal-to-Quantization-Noise Ratio (SQNR) in dB.
    Works seamlessly with both real and complex numpy arrays.
    """
    # Ensure inputs are numpy arrays
    ref = np.asarray(ref_vector)
    dut = np.asarray(dut_vector)

    if ref.shape != dut.shape:
        raise ValueError(f"Shape mismatch: REF is {ref.shape}, DUT is {dut.shape}")

    # Calculate error (Quantization Noise)
    error = dut - ref

    # Calculate power (Sum of squared magnitudes)
    # np.abs() automatically handles the sqrt(re^2 + im^2) for complex numbers
    signal_power = np.sum(np.abs(ref)**2)
    noise_power = np.sum(np.abs(error)**2)

    # Handle perfect bit-true matches to avoid divide-by-zero errors
    if noise_power == 0:
        return float('inf')

    # Compute SQNR in Decibels
    sqnr_db = 10 * np.log10(signal_power / noise_power)
    
    return sqnr_db

# =========================================================================
# Example Usage: Comparing your MATLAB Reference vs RTL Output
# =========================================================================
if __name__ == "__main__":
    
    print("Loading data files...")
    
    # 1. Load the MATLAB Floating-Point Golden Reference
    # (These are the files we just generated in the final MATLAB script)
    try:
        ref_re = np.loadtxt("py_gold_out_re_q8_8.txt")
        ref_im = np.loadtxt("py_gold_out_im_q8_8.txt")
        reference_signal = ref_re + 1j * ref_im
    except Exception as e:
        print(f"Error loading reference files: {e}")
        exit()

    # 2. Load your DUT (RTL) Output
    # (Assuming you dumped your RTL IFFT output to text files during simulation)
    # Replace these filenames with whatever your UVM monitor spits out!
    try:
        dut_re = np.loadtxt("matlab_out_re_q8_8.txt")
        dut_im = np.loadtxt("matlab_out_im_q8_8.txt")
        dut_signal = dut_re + 1j * dut_im
    except Exception as e:
        print(f"Error loading DUT files: {e}")
        exit()

    # 3. Calculate and Report
    print(f"Comparing {len(reference_signal)} samples...")
    sqnr = calculate_sqnr(reference_signal, dut_signal)
    
    if sqnr == float('inf'):
        print("SQNR: INF dB (Perfect bit-true match!)")
    else:
        print(f"SQNR: {sqnr:.2f} dB")
        
        THRESHOLD_DB = 30.0
        if sqnr >= THRESHOLD_DB:
            print("Verdict: [PASS]")           # <-- Removed the checkmark
        else:
            print(f"Verdict: [FAIL] (Below {THRESHOLD_DB} dB threshold)") # <-- Removed the cross