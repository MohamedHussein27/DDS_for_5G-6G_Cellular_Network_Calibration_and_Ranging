function T = mytypes(dt)
%MYTYPES Defines data types. 
%   Uses coder.inline('always') to ensure types are compiled directly.

    coder.inline('always'); 

    switch dt
        case 'double'
            % --- SIGNAL TYPES ---
            T.M        = double([]);
            T.acc      = double([]);
            T.phase    = double([]);
            T.lut_addr = double([]);
            T.sine_LUT = double([]);
            T.dds_out  = double([]);
            T.f_target_inst = double([]);
            
            % --- STRUCTURAL VARIABLES (Always Double) ---
            T.Ns         = double([]);
            T.LUT_size   = double([]);
            T.addr_ROM   = double([]);
            T.addr_shift = double([]);
            
        case 'single'
            % --- SIGNAL TYPES ---
            T.M        = single([]);
            T.acc      = single([]);
            T.phase    = single([]);
            T.lut_addr = single([]);
            T.sine_LUT = single([]);
            T.dds_out  = single([]);
            T.f_target_inst = single([]);
            
            % --- STRUCTURAL VARIABLES (Always Double) ---
            T.Ns         = double([]);
            T.LUT_size   = double([]);
            T.addr_ROM   = double([]);
            T.addr_shift = double([]);
            
        case 'fixed'
            % --- FIMATH (Behavior) ---
            F = fimath('RoundingMethod', 'Floor', 'OverflowAction', 'Wrap');
            
            % --- SIGNAL TYPES ---
            % M, Acc, Phase: Unsigned 32-bit Integers (Counters)
            T.M        = fi(0, 0, 32, 0, F); 
            T.acc      = fi(0, 0, 32, 0, F);
            T.phase    = fi(0, 0, 32, 0, F);
            
            % LUT Address: Unsigned 16-bit Integer
            T.lut_addr = fi(0, 0, 16, 0, F);
            
            % Sine Output: Signed 16-bit Fraction (Range -1 to 1)
            T.sine_LUT = fi(0, 1, 16, 15, F);
            T.dds_out  = fi(0, 1, 16, 15, F);
            
            % --- STRUCTURAL VARIABLES (Must remain Double) ---
            % These are used for array indexing and loop bounds.
            T.Ns         = double([]);
            T.LUT_size   = double([]);
            T.addr_ROM   = double([]);
            T.addr_shift = double([]);
            
        otherwise
            error('Unknown data type mode');
    end
end