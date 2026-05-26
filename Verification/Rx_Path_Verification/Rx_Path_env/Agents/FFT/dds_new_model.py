import numpy as np

def generate_dds_golden_model(FTW_start=0, FTW_step=426666, N_cycles=4096):
    """
    Bit-True DDS Hardware Emulator.
    Translated directly from the updated MATLAB sequential logic.
    Generates the LUT internally and mimics exact combinational address mapping.
    """
    ADDRESS_WIDTH = 16
    MEMORY_WIDTH = 8
    
    # =========================================================
    # 1. INTERNAL LUT GENERATION (Quarter-Wave)
    # =========================================================
    lut_depth = 2**(ADDRESS_WIDTH - 2) 
    n_lut = np.arange(lut_depth)
    
    # Matches MATLAB: round((2^(MEMORY_WIDTH-1)-1) * sin((pi/2) * (0:lut_depth-1) / lut_depth))
    lut_table = np.round((2**(MEMORY_WIDTH - 1) - 1) * np.sin((np.pi / 2) * n_lut / lut_depth)).astype(np.int32)
    
    # =========================================================
    # 2. SEQUENTIAL LOGIC (Phase Accumulation)
    # =========================================================
    # The MATLAB loop models clock cycles i=1 to N.
    # Mathematically, the accumulation at step 'i' is: FTW_start*i + FTW_step*[i(i-1)/2]
    i = np.arange(1, N_cycles + 1, dtype=np.uint64)
    accumulation_term = (i * (i - np.uint64(1))) // np.uint64(2)
    
    phase_acc = (np.uint64(FTW_start) * i) + (np.uint64(FTW_step) * accumulation_term)
    
    # Modulo 2^32 wrap (matches: mod(..., 2^32))
    truncated_phase = (phase_acc & 0xFFFFFFFF).astype(np.uint32)
    
    # =========================================================
    # 3. COMBINATIONAL LOGIC (Quadrant Mapping)
    # =========================================================
    # Matches MATLAB bitshift logic
    phase_16 = truncated_phase >> (32 - ADDRESS_WIDTH)
    quadrant = phase_16 >> 14
    
    mapped_addr = np.zeros_like(phase_16, dtype=np.uint32)
    neg_flag = np.zeros_like(phase_16, dtype=bool)
    
    # Case 0: 1st Quad
    q0 = (quadrant == 0)
    mapped_addr[q0] = phase_16[q0]
    neg_flag[q0] = False
    
    # Case 1: 2nd Quad
    q1 = (quadrant == 1)
    mapped_addr[q1] = 32768 - 1 - phase_16[q1]
    neg_flag[q1] = False
    
    # Case 2: 3rd Quad
    q2 = (quadrant == 2)
    mapped_addr[q2] = phase_16[q2] - 32768
    neg_flag[q2] = True
    
    # Case 3: 4th Quad
    q3 = (quadrant == 3)
    mapped_addr[q3] = 65536 - 1 - phase_16[q3]
    neg_flag[q3] = True
    
    # =========================================================
    # 4. LUT LOOKUP & OUTPUT
    # =========================================================
    # Python is naturally 0-indexed, so we omit the "+ 1" offset used in MATLAB.
    amp = lut_table[mapped_addr]
    
    # Apply negative flag
    chirp_time = np.where(neg_flag, -amp, amp).astype(np.int64)
    
    return chirp_time

# ===========================================================================
# Example Execution
# ===========================================================================
if __name__ == "__main__":
    
    # Run the function (No mem_file needed anymore!)
    ideal_chirp = generate_dds_golden_model(
        FTW_start=0, 
        FTW_step=426666, 
        N_cycles=4096
    )
    
    print(f"Generated {len(ideal_chirp)} bit-true DDS samples.")
    print(f"First 10 samples: {ideal_chirp[:10]}")