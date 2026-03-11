function dds_out = dds_core(M_start, M_step, Ns, Nacc, LUT_bits, DT_Mode) %#codegen
%DDS_CORE Hardware-accurate DDS model with Integrated Frequency Sweep
%
% Inputs:
%   M_start  : Tuning word for the start frequency (f0)
%   M_step   : The amount to increase the tuning word every clock cycle
%   Ns       : Total number of samples to generate
%   Nacc     : Bit-width of the Phase Accumulator
%   LUT_bits : Bit-width of the Sine LUT address
%   DT_Mode  : 'fixed', 'single', or 'double'

    % -----------------------------------------------------
    % 1. Load Types
    % -----------------------------------------------------
    T = mytypes(DT_Mode);
    assert(LUT_bits <= Nacc, 'LUT_bits must be <= Nacc');
    
    % -----------------------------------------------------
    % 2. Cast Inputs
    % -----------------------------------------------------
    % Cast the HW parameters to the defined precision type
    M_start_cast = cast(M_start, 'like', T.M);
    M_step_cast  = cast(M_step, 'like', T.M);
    
    % -----------------------------------------------------
    % 3. Sine LUT (Phase-to-Amplitude Converter)
    % -----------------------------------------------------
    LUT_size = 2^LUT_bits;
    addr_ROM = 0:LUT_size-1;
    sine_LUT = cast(sin(2*pi*addr_ROM/LUT_size), 'like', T.sine_LUT);
    
    % -----------------------------------------------------
    % 4. Dual Accumulators (Frequency and Phase)
    % -----------------------------------------------------
    phase_acc = uint64(0);                   % Accumulates Phase
    freq_acc  = uint64(M_start_cast);        % Accumulates Frequency
    
    phase = zeros(1, Ns, 'uint64');
    
    % --- HARDWARE PIPELINE SIMULATION ---
    for n = 1:Ns
        % Step A: Read current phase for the LUT
        phase(n) = phase_acc;
        
        % Step B: Update Phase Accumulator
        % Add the *current* instantaneous frequency to the phase
        phase_acc = phase_acc + freq_acc;
        phase_acc = bitand(phase_acc, uint64(2^Nacc - 1)); % Hardware wrap-around
        
        % Step C: Update Frequency Accumulator
        % Add the step size to increase the frequency for the next clock cycle
        freq_acc = freq_acc + uint64(M_step_cast);
    end
    
    % -----------------------------------------------------
    % 5. Phase Truncation
    % -----------------------------------------------------
    addr_shift = Nacc - LUT_bits;
    lut_addr = bitshift(phase, -addr_shift);
    lut_addr = cast(lut_addr, 'like', T.lut_addr);
    
    % -----------------------------------------------------
    % 6. ROM Lookup
    % -----------------------------------------------------
    dds_out = zeros(1, Ns, 'like', T.dds_out);
    dds_out(:) = sine_LUT(double(lut_addr) + 1);
end