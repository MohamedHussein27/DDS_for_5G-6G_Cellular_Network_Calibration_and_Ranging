function Y_hw = qxy_fft_multiplier(A_vec, B_vec, x, y, outType)
%% ---------------------------------------------------------
% Qx.y × Qx.y Wallace Tree Multiplier (Vector Version)
%
% Supports:
%   • Real & Complex inputs
%   • double / single (native MATLAB)
%   • fixed-point (Wallace tree, bit-true)
% ---------------------------------------------------------

    if nargin < 5
        outType = 'double';
    end

    assert(length(A_vec) == length(B_vec), ...
           'Input vectors must have same length');

    Npts = length(A_vec);

    %% =====================================================
    % COMPLEX HANDLING (4 real Wallace multipliers)
    %% =====================================================
    if ~isreal(A_vec) || ~isreal(B_vec)

        Ar = real(A_vec);  Ai = imag(A_vec);
        Br = real(B_vec);  Bi = imag(B_vec);

        P1 = qxy_fft_multiplier(Ar, Br, x, y, outType);
        P2 = qxy_fft_multiplier(Ai, Bi, x, y, outType);
        P3 = qxy_fft_multiplier(Ar, Bi, x, y, outType);
        P4 = qxy_fft_multiplier(Ai, Br, x, y, outType);

        Y_hw = complex(P1 - P2, P3 + P4);
        return;
    end

    %% =====================================================
    % DOUBLE / SINGLE
    %% =====================================================
    if strcmpi(outType,'double')
        Y_hw = double(A_vec .* B_vec);
        return;
    elseif strcmpi(outType,'single')
        Y_hw = single(A_vec .* B_vec);
        return;
    end

    %% =====================================================
    % FIXED-POINT WALLACE TREE
    %% =====================================================
    WordIn  = x + y;
    Scale   = 2^y;
    WordOut = 2 * WordIn;
    FracOut = 2 * y;

    T = numerictype(1, WordOut, FracOut);
    F = fimath('OverflowAction','Wrap', ...
               'RoundingMethod','Floor');

    Y_hw = fi(zeros(Npts,1), T, F);

    %% =====================================================
    % PROCESS VECTOR
    %% =====================================================
    for n = 1:Npts

        A_val = A_vec(n);
        B_val = B_vec(n);

        % ----- Sign handling -----
        A_sign = A_val < 0;
        B_sign = B_val < 0;
        result_sign = xor(A_sign, B_sign);

        % ----- FIX: convert to native integers -----
        A_int = uint64(round(double(abs(A_val)) * Scale));
        B_int = uint64(round(double(abs(B_val)) * Scale));

        % ----- Partial products -----
        PP = zeros(WordIn, WordIn);

        for i = 1:WordIn
            for j = 1:WordIn
                PP(i,j) = bitand(bitget(A_int,i), ...
                                 bitget(B_int,j));
            end
        end

        % ----- Wallace tree reduction -----
        cols = cell(1, WordOut);

        for i = 1:WordIn
            for j = 1:WordIn
                cols{i+j-1} = [cols{i+j-1}, PP(i,j)];
            end
        end

        while max(cellfun(@length, cols)) > 2

            new_cols = cell(1, WordOut+1);

            for k = 1:WordOut
                bits = cols{k};

                while numel(bits) >= 3
                    b1 = bits(end); bits(end)=[];
                    b2 = bits(end); bits(end)=[];
                    b3 = bits(end); bits(end)=[];

                    s = xor(xor(b1,b2),b3);
                    c = (b1&b2) | (b1&b3) | (b2&b3);

                    new_cols{k}   = [new_cols{k}, s];
                    new_cols{k+1} = [new_cols{k+1}, c];
                end

                new_cols{k} = [new_cols{k}, bits];
            end

            cols = new_cols(1:WordOut);
        end

        % ----- Final ripple-carry adder -----
        sum_row = zeros(1, WordOut);
        carry = 0;

        for k = 1:WordOut
            s = sum(cols{k}) + carry;
            sum_row(k) = mod(s,2);
            carry = floor(s/2);
        end

        % ----- Reconstruct integer -----
        raw_product = 0;
        for k = 1:WordOut
            raw_product = raw_product + sum_row(k)*2^(k-1);
        end

        % ----- Apply sign -----
        if result_sign
            raw_product = 2^WordOut - raw_product;
        end
        raw_product = mod(raw_product, 2^WordOut);

        % ----- Store fixed-point result -----
        temp_unsigned = fi(raw_product, numerictype(0,WordOut,0));
        Y_hw(n) = reinterpretcast(temp_unsigned, T);
    end
end
