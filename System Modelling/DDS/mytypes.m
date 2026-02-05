function T = mytypes(dt)

    % global var for differnt Q notation sweeping
    global DDS_FRAC_BITS;
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
            if isempty(DDS_FRAC_BITS)
                fracBits = 15;   % Safe default
            else
                fracBits = DDS_FRAC_BITS;
            end
        
            F = fimath('RoundingMethod','Floor', ...
                       'OverflowAction','Wrap');
        
            % Phase path (unchanged)
            T.M     = fi(0, 0, 32, 0, F);
            T.acc   = fi(0, 0, 32, 0, F);
            T.phase = fi(0, 0, 32, 0, F);
        
            % LUT address
            T.lut_addr = fi(0, 0, 16, 0, F);
        
            % Sine amplitude (Q1.fracBits)
            WL = fracBits + 1;
            T.sine_LUT = fi(0, 1, WL, fracBits, F);
            T.dds_out  = fi(0, 1, WL, fracBits, F);


            
        otherwise
            error('Unknown data type mode');
    end
end