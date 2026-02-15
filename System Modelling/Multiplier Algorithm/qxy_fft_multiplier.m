function Y_hw = qxy_fft_multiplier(A_vec, B_vec, x, y, outType)
%% ---------------------------------------------------------
% Qx.y × Qx.y Wallace Tree Multiplier (Vector Version)
%
% Features:
%   • Supports REAL and COMPLEX inputs
%   • Double / Single → native MATLAB multiply
%   • Fixed → Wallace Tree hardware-style multiply
%
% For complex numbers:
%   (Ar + jAi)(Br + jBi)
%       = (ArBr - AiBi) + j(ArBi + AiBr)
%   → Implemented using FOUR real Wallace multipliers
%
% Inputs:
%   A_vec   : Input vector
%   B_vec   : Input vector
%   x       : Integer bits including sign
%   y       : Fractional bits
%   outType : 'double' | 'single' | 'fixed'
%
% Output:
%   Y_hw    : Output vector
% ---------------------------------------------------------

    if nargin < 5
        outType = 'double';
    end

    % Ensure vectors match
    assert(length(A_vec) == length(B_vec), ...
           'Input vectors must have same length');

    Npts = length(A_vec);

    %% =====================================================
    %  COMPLEX HANDLING (Decomposition into real parts)
    % ======================================================
    %
    % Wallace tree is REAL-only hardware.
    % So if inputs are complex:
    %   → Split into real & imaginary
    %   → Use 4 real multipliers
    %   → Recombine
    %
    if ~isreal(A_vec) || ~isreal(B_vec)

        % Separate components
        Ar = real(A_vec);  Ai = imag(A_vec);
        Br = real(B_vec);  Bi = imag(B_vec);

        % Four real multiplications (reuse SAME function)
        P1 = qxy_fft_multiplier(Ar, Br, x, y, outType); % Ar*Br
        P2 = qxy_fft_multiplier(Ai, Bi, x, y, outType); % Ai*Bi
        P3 = qxy_fft_multiplier(Ar, Bi, x, y, outType); % Ar*Bi
        P4 = qxy_fft_multiplier(Ai, Br, x, y, outType); % Ai*Br

        % Combine results
        Y_real = P1 - P2;
        Y_imag = P3 + P4;

        % Final complex output
        Y_hw = complex(Y_real, Y_imag);
        return;
    end


    %% =====================================================
    %  DOUBLE / SINGLE PATH
    % ======================================================
    if strcmpi(outType,'double')
        Y_hw = double(A_vec .* B_vec);
        return;

    elseif strcmpi(outType,'single')
        Y_hw = single(A_vec .* B_vec);
        return;
    end


    %% =====================================================
    %  FIXED-POINT WALLACE TREE PATH (Qx.y)
    % ======================================================

    WordIn  = x + y;     % Total input bits
    Scale   = 2^y;       % Scaling factor for fraction
    WordOut = 2*WordIn;  % Output width (full precision)
    FracOut = 2*y;       % Fractional bits after multiply

    % Output fixed-point type
    T = numerictype(1, WordOut, FracOut);
    F = fimath('OverflowAction','Wrap', ...
               'RoundingMethod','Floor');

    Y_hw = fi(zeros(Npts,1), T, F);

    %% =====================================================
    %  PROCESS EACH VECTOR ELEMENT
    % ======================================================
    for n = 1:Npts

        A_val = A_vec(n);
        B_val = B_vec(n);

        %% -----------------------------
        % SIGN EXTRACTION
        % -----------------------------
        A_sign = A_val < 0;
        B_sign = B_val < 0;
        result_sign = xor(A_sign, B_sign);

        % Convert to scaled integers (magnitude only)
        A_int = round(abs(A_val) * Scale);
        B_int = round(abs(B_val) * Scale);

        %% -----------------------------
        % PARTIAL PRODUCT GENERATION
        % -----------------------------
        %
        % PP(i,j) = bit_i(A) AND bit_j(B)
        %
        PP = zeros(WordIn, WordIn);

        for i = 1:WordIn
            for j = 1:WordIn
                PP(i,j) = bitand(bitget(A_int,i), ...
                                 bitget(B_int,j));
            end
        end


        %% -----------------------------
        % WALLACE TREE REDUCTION
        % -----------------------------
        %
        % 3:2 compressors (Carry Save Adders)
        % Reduce each column until <= 2 bits remain
        %
        cols = cell(1, WordOut);

        % Arrange partial products into columns
        for i = 1:WordIn
            for j = 1:WordIn
                idx = i + j - 1;
                cols{idx} = [cols{idx}, PP(i,j)];
            end
        end

        % Iterative CSA compression
        while max(cellfun(@length, cols)) > 2

            new_cols = cell(1, WordOut+1);

            for k = 1:WordOut
                bits = cols{k};

                while length(bits) >= 3
                    b1 = bits(end); bits(end)=[];
                    b2 = bits(end); bits(end)=[];
                    b3 = bits(end); bits(end)=[];

                    % 3:2 compressor
                    s = xor(xor(b1,b2),b3);
                    c = (b1&b2) | (b1&b3) | (b2&b3);

                    new_cols{k}   = [new_cols{k}, s];
                    new_cols{k+1} = [new_cols{k+1}, c];
                end

                new_cols{k} = [new_cols{k}, bits];
            end

            cols = new_cols(1:WordOut);
        end


        %% -----------------------------
        % FINAL RIPPLE-CARRY ADDER
        % -----------------------------
        %
        % Convert two rows into final binary result
        %
        sum_row = zeros(1, WordOut);
        carry = 0;

        for k = 1:WordOut
            s = sum(cols{k}) + carry;
            sum_row(k) = mod(s,2);
            carry = floor(s/2);
        end


        %% -----------------------------
        % RECONSTRUCT INTEGER VALUE
        % -----------------------------
        raw_product = 0;

        for k = 1:WordOut
            raw_product = raw_product + ...
                          sum_row(k) * 2^(k-1);
        end


        %% -----------------------------
        % APPLY SIGN (Two's Complement)
        % -----------------------------
        if result_sign
            raw_product = 2^WordOut - raw_product;
        end

        raw_product = mod(raw_product, 2^WordOut);


        %% -----------------------------
        % STORE FIXED-POINT RESULT
        % -----------------------------
        temp_unsigned = fi(raw_product, ...
                           numerictype(0, WordOut, 0));

        out_fmt = numerictype(1, WordOut, FracOut);

        Y_hw(n) = reinterpretcast(temp_unsigned, out_fmt);
    end
end
