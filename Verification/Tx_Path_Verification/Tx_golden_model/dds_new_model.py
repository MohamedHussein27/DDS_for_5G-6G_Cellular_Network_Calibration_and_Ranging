import numpy as np
import os

def generate_dds_golden_model(FTW_start=0, FTW_step=426666, N_cycles=4096, Fs=491.52e6, mem_file_path="memfile.mem"):
    """
    Bit-True DDS Hardware Emulator.
    Matches the exact truncation, quadrant mapping, and ROM fetching of the Verilog RTL.
    
    Args:
        FTW_start: Starting Frequency Tuning Word
        FTW_step: Frequency Tuning Word step per clock cycle
        N_cycles: Number of clock cycles to simulate
        Fs: System clock frequency (Hz)
        mem_file_path: Path to the Verilog memory file containing the ROM hex data
        
    Returns:
        ideal_out (numpy array): The bit-true time-domain chirp
        B (float): The calculated physical bandwidth reached in Hz
    """
    
    # --- 1. Load ROM Data ---
    if not os.path.exists(mem_file_path):
        raise FileNotFoundError(f"Could not find the ROM file: {mem_file_path}")
        
    # Read the exact hex strings from the Verilog memory file
    with open(mem_file_path, 'r') as f:
        # Ignore empty lines and convert hex string to integer
        rom_data = np.array([int(line.strip(), 16) for line in f if line.strip()], dtype=np.int32)
        
    # --- 2. Calculate Physical Bandwidth ---
    max_ftw = FTW_start + (FTW_step * N_cycles)
    B = ((max_ftw - FTW_start) * Fs) / (2**32)
    
    # --- 3. Generate Ideal Discrete-Time Model (BIT-TRUE DIRECT FETCH) ---
    # Create clock cycle index array: 1 to N_cycles (Shifts the array left by 1)
    n = np.arange(1, N_cycles + 1, dtype=np.uint64)
    
    # Discrete Phase Accumulation (Formula: FTW_start*n + FTW_step*[n(n-1)/2])
    accumulation_term = (n * (n - np.uint64(1))) // np.uint64(2)
    phase_acc = (np.uint64(FTW_start) * n) + (np.uint64(FTW_step) * accumulation_term)
    
    # Apply 32-bit wrapping mask (mod 2^32)
    ideal_discrete_phase_word = (phase_acc & 0xFFFFFFFF).astype(np.uint32)
    
    # ---- EMULATE PHASE_ACC.v ----
    # Truncate to ADDRESS_WIDTH = 16 bits (drop the bottom 16 bits)
    truncated_phase_16b = ideal_discrete_phase_word >> 16
    
    # ---- EMULATE Quadrant Mapper & LUT ----
    # Extract top 2 bits for quadrant, lower 14 bits for the mapped address
    quadrant = truncated_phase_16b >> 14
    addr = truncated_phase_16b & 0x3FFF  # 0x3FFF is 16383 in decimal
    
    # Quadrant Mapping Logic (Emulating first_quad_address.v)
    # 1st quad (0) and 3rd quad (2): mapped_addr = addr
    # 2nd quad (1) and 4th quad (3): mapped_addr = 16383 - addr
    mapped_addr = np.where((quadrant == 1) | (quadrant == 3), 
                           16383 - addr, 
                           addr)
    
    # Sign Mapping (Emulating negative_mux.v)
    # 3rd quad (2) and 4th quad (3) are negative
    neg_flag = np.where((quadrant == 2) | (quadrant == 3), 1, 0)
    
    # ---- EMULATE LUT.v (Direct ROM Read) ----
    # Python is naturally 0-indexed, so we do not need the "+ 1" offset used in MATLAB
    lut_amplitude = rom_data[mapped_addr]
    
    # ---- Apply Negative Mux ----
    ideal_out = np.where(neg_flag == 1, -lut_amplitude, lut_amplitude)
    
    return ideal_out, B

# ===========================================================================
# Example Execution
# ===========================================================================
if __name__ == "__main__":
    try:
        # Run the function
        ideal_chirp, bandwidth = generate_dds_golden_model(
            FTW_start=0, 
            FTW_step=426666, 
            N_cycles=4096, 
            Fs=491.52e6, 
            mem_file_path="memfile.mem"
        )
        
        print(f"Calculated Physical Bandwidth: {bandwidth / 1e6:.2f} MHz")
        print(f"Generated {len(ideal_chirp)} bit-true DDS samples.")
        
        # Sneak peek at the first 10 samples
        print(f"First 10 samples: {ideal_chirp[:12]}")
        
    except FileNotFoundError as e:
        print(e)