import numpy as np

def generate_mock_rom():
    """
    Generates a mathematical quarter-wave ROM on the fly.
    Provides the exact 16,384 values needed for the LUT.
    """
    addresses = np.arange(16384)
    # 14-bit addressing mapped to 16-bit phase for a quarter wave
    sine_float = np.sin(2.0 * np.pi * addresses / 65536.0)
    # Assuming 8-bit output (signed max 127)
    rom_data = np.round(sine_float * 127).astype(np.int32)
    return rom_data

def generate_expected_chirp(ftw_start, ftw_step, cycles, rom_data):
    """
    Bit-true Python reference model for LFM Chirp.
    """
    # ---- 1. EMULATE PHASE_ACC.v ----
    n = np.arange(cycles, dtype=np.int64)
    ideal_phase = (ftw_start * n) + (ftw_step * (n * (n - 1)) // 2)
    ideal_phase = ideal_phase % (2**32) 
    truncated_phase_16b = ideal_phase // 65536 

    # ---- 2. EMULATE Quadrant Mapper & LUT ----
    ideal_out = np.zeros(cycles, dtype=np.int32)

    for k in range(cycles):
        current_phase = truncated_phase_16b[k]
        
        quadrant = current_phase // 16384
        addr = current_phase % 16384
        
        if quadrant == 0:     # 1st quad
            mapped_addr = addr
            neg_flag = 0
        elif quadrant == 1:   # 2nd quad
            mapped_addr = 16383 - addr
            neg_flag = 0
        elif quadrant == 2:   # 3rd quad
            mapped_addr = addr
            neg_flag = 1
        else:                 # 4th quad
            mapped_addr = 16383 - addr
            neg_flag = 1
            
        # ---- 3. EMULATE LUT.v ----
        lut_amplitude = rom_data[mapped_addr]
        
        # ---- 4. EMULATE negative_mux.v ----
        if neg_flag == 1:
            ideal_out[k] = -lut_amplitude
        else:
            ideal_out[k] = lut_amplitude
            
    return ideal_out