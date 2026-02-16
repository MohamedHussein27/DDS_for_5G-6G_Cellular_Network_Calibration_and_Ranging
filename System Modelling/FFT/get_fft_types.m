function T = get_fft_types(mode, WL, FL)
% GET_FFT_TYPES Returns a struct of data types based on the selected mode.
% mode: 'Double', 'Single', or 'Fixed'
% WL: Word Length (for fixed point)
% FL: Fractional Length (for fixed point)

    switch mode
        case 'Double'
            T.Input  = 'double';
            
        case 'Single'
            T.Input  = 'single';
            
        case 'Fixed'
            % Define Fixed-Point properties (Signed, WordLength, FractionLength)
            % Rounding: Nearest, Overflow: Saturate (Hardware friendly)
            T.Input = numerictype(1, WL, FL);
            
        otherwise
            error('Unknown mode');
    end
end