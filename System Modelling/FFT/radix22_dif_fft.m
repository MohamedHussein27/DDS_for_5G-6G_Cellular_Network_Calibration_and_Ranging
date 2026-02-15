function X = radix22_dif_fft(x)
%#codegen
% RADIX22_DIF_FFT Core Algorithm
% Hardware-friendly Radix-2² DIF FFT implementation.
% Supports: Double, Single, and Fixed-Point (via polymorphism).

    N = length(x);
    % Force column vector and inherit type from input 'x'
    X = x(:);
    
    stages = log2(N);
    
    for s = stages:-1:1
        m = 2^s;
        mh = m/2;
        
        % Twiddle Factors
        % We compute in double, then cast to the input type (e.g., fixed-point)
        W_float = exp(-1j*2*pi*(0:mh-1)/m).';
        W = cast(W_float, 'like', x); 
        
        for k = 1:m:N
            % Butterfly Operations
            a = X(k:k+mh-1);
            b = X(k+mh:k+m-1);
            
            % Radix-2 Butterfly
            X(k:k+mh-1)   = a + b;
            X(k+mh:k+m-1) = (a - b) .* W;
        end
    end
    
    % Bit-Reverse Output
    X = bitrevorder(X);
end