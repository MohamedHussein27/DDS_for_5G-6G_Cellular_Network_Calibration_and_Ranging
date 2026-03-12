function dds_out = dds_core(M_start, M_step, Ns, Nacc, LUT_bits, DT_Mode) %#codegen
%DDS_CORE Hardware-accurate DDS model with Lossless Quarter-Wave Symmetry

    % -----------------------------------------------------
    % 1. Load Types
    % -----------------------------------------------------
    T = mytypes(DT_Mode);
    assert(LUT_bits >= 2, 'LUT_bits must be >= 2 for quarter-wave logic');
    assert(LUT_bits <= Nacc, 'LUT_bits must be <= Nacc');
    
    % Cast the HW parameters
    M_start_cast = cast(M_start, 'like', T.M);
    M_step_cast  = cast(M_step, 'like', T.M);
    
    % -----------------------------------------------------
    % 2. Quarter-Wave Sine LUT (N/4 + 1 Rule)
    % -----------------------------------------------------
    % A full wave has 2^LUT_bits points (e.g., 256).
    % A quarter wave has 2^(LUT_bits - 2) points (e.g., 64).
    Q_LUT_bits = LUT_bits - 2; 
    Q_LUT_size = 2^Q_LUT_bits;
    
    % CRITICAL FIX: We store from 0 up to Q_LUT_size (65 points total).
    % This guarantees we store the exact peak (1.0) at pi/2.
    addr_ROM = 0 : Q_LUT_size;
    
    % Generate the quarter wave using the full circle denominator
    sine_LUT = cast(sin(2*pi*addr_ROM / (2^LUT_bits)), 'like', T.sine_LUT);
    
    % -----------------------------------------------------
    % 3. Dual Accumulators
    % -----------------------------------------------------
    phase_acc = uint64(0);
    freq_acc  = uint64(M_start_cast);
    phase     = zeros(1, Ns, 'uint64');
    
    % --- HARDWARE PIPELINE SIMULATION ---
    for n = 1:Ns
        phase(n)  = phase_acc;
        
        phase_acc = phase_acc + freq_acc;
        phase_acc = bitand(phase_acc, uint64(2^Nacc - 1)); 
        
        freq_acc  = freq_acc + uint64(M_step_cast);
    end
    
    % -----------------------------------------------------
    % 4. Phase Truncation
    % -----------------------------------------------------
    addr_shift = Nacc - LUT_bits;
    lut_addr   = bitshift(phase, -addr_shift);
    
    % -----------------------------------------------------
    % 5. Quarter-Wave Decoding Logic
    % -----------------------------------------------------
    dds_out = zeros(1, Ns, 'like', T.dds_out);
    offset_mask = uint64(Q_LUT_size - 1); % e.g., 63
    
    for n = 1:Ns
        current_addr = uint64(lut_addr(n));
        
        % A. Extract top 2 bits for quadrant (0, 1, 2, or 3)
        quadrant = bitshift(current_addr, -Q_LUT_bits);
        
        % B. Extract remaining bits for offset
        offset = bitand(current_addr, offset_mask);
        
        % C. Mirroring logic for odd quadrants (Q2 and Q4)
        if bitand(quadrant, uint64(1)) == uint64(1)
            % Move backward from the peak (e.g., 64 - offset)
            read_addr = uint64(Q_LUT_size) - offset;
        else
            % Move forward from zero
            read_addr = offset;
        end
        
        % D. Fetch absolute value from the Quarter-LUT
        % Because MATLAB is 1-indexed, we add 1 to the read_addr
        val = sine_LUT(double(read_addr) + 1);
        
        % E. Negation logic for bottom half of the circle (Q3 and Q4)
        if bitand(quadrant, uint64(2)) == uint64(2)
            dds_out(n) = -val;
        else
            dds_out(n) = val;
        end
    end
end