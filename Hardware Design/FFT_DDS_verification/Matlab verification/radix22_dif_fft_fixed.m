function X = radix22_dif_fft_fixed(x, WL, FL)
% Fixed-point Radix-2² DIF FFT with 2 Guard Bits
% Models native 18-bit FPGA Block RAM for maximum 16-bit exit SQNR

    N = length(x);
    stages = log2(N);

    % 1. The Strict 16-bit Exit Door
    T_out = numerictype(1, WL, FL);   

    % 2. The 18-bit Internal Guard-Bit RAM (Free FPGA Hardware!)
    % We add 2 bits of word length, assigned purely to the fraction
    WL_internal = WL + 2;
    FL_internal = FL + 2;
    T_internal = numerictype(1, WL_internal, FL_internal);

    % Use Convergent rounding to kill DC bias
    F = fimath( ...
        'RoundingMethod', 'Convergent', ...
        'OverflowAction', 'Wrap', ...
        'ProductMode', 'FullPrecision', ...
        'SumMode', 'FullPrecision');

    % 3. Initialize X into the 18-bit Guard RAM
    X = fi(x, T_internal, F);

    % Twiddles use the 16-bit optimized precision we found earlier
    T_twiddle = numerictype(1, WL, WL - 2);

    for s = stages:-1:1
        m  = 2^s;
        mh = m/2;

        W_raw = exp(-1j*2*pi*(0:mh-1)/m).';
        W = fi(W_raw, T_twiddle, F);

        for k = 1:m:N
            a = X(k:k+mh-1);
            b = X(k+mh:k+m-1);

            % Math happens in FullPrecision, rounds back to 18-bit Guard RAM
            X(k:k+mh-1)   = a + b;
            X(k+mh:k+m-1) = (a - b) .* W;
        end
    end

    % Bit-reversed output
    X = bitrevorder(X);

    % 4. The Final Squeeze
    % Chop off the 2 guard bits to cleanly exit the block as strict 16-bit
    X = fi(X, T_out, F); 
end













