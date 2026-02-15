function Y_hw = qxy_fft_multiplier(A_vec, B_vec, x, y, outType)
%% ---------------------------------------------------------
% Qx.y × Qx.y Wallace Tree Multiplier (Vector Version)
% - Fixed-point uses Qx.y representation and Wallace Tree
% - Double / Single uses native MATLAB arithmetic
%
% Inputs:
%   A_vec   : Input vector 
%   B_vec   : Input vector 
%   x       : Integer bits including sign (used only for fixed)
%   y       : Fractional bits (used only for fixed)
%   outType : 'double' | 'single' | 'fixed'
%
% Output:
%   Y_hw    : Output vector in requested format
% ---------------------------------------------------------

    if nargin < 5
        outType = 'double';
    end

    assert(length(A_vec) == length(B_vec), 'Input vectors must have same length');
    Npts = length(A_vec);

    % -----------------------------
    % 1. Handle double/single directly
    % -----------------------------
    if strcmpi(outType,'double')
        Y_hw = double(A_vec .* B_vec);
        return;
    elseif strcmpi(outType,'single')
        Y_hw = single(A_vec .* B_vec);
        return;
    end

    % -----------------------------
    % 2. Fixed-point path (Qx.y)
    % -----------------------------
    WordIn = x + y;         % input total bits
    Scale = 2^y;            % fractional scaling
    WordOut = 2*WordIn;     % output width
    FracOut = 2*y;          % fractional bits

    % Preallocate fixed-point output
    T = numerictype(1, WordOut, FracOut);
    F = fimath('OverflowAction','Wrap','RoundingMethod','Floor');
    Y_hw = fi(zeros(Npts,1), T, F);

    % Loop per element
    for n = 1:Npts
        A_val = A_vec(n);
        B_val = B_vec(n);

        % -----------------------------
        % Sign & magnitude
        % -----------------------------
        A_sign = A_val < 0;
        B_sign = B_val < 0;
        result_sign = xor(A_sign, B_sign);

        A_int = round(abs(A_val) * Scale);
        B_int = round(abs(B_val) * Scale);

        % -----------------------------
        % Partial products
        % -----------------------------
        PP = zeros(WordIn, WordIn);
        for i = 1:WordIn
            for j = 1:WordIn
                PP(i,j) = bitand(bitget(A_int,i), bitget(B_int,j));
            end
        end

        % -----------------------------
        % Wallace Tree CSA reduction
        % -----------------------------
        cols = cell(1, WordOut);
        for i = 1:WordIn
            for j = 1:WordIn
                idx = i+j-1;
                cols{idx} = [cols{idx}, PP(i,j)];
            end
        end

        while max(cellfun(@length, cols)) > 2
            new_cols = cell(1, WordOut+1);
            for k = 1:WordOut
                bits = cols{k};
                while length(bits) >= 3
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

        % -----------------------------
        % Final ripple-carry adder
        % -----------------------------
        sum_row = zeros(1, WordOut);
        carry = 0;
        for k = 1:WordOut
            s = sum(cols{k}) + carry;
            sum_row(k) = mod(s,2);
            carry = floor(s/2);
        end

        % -----------------------------
        % Reconstruct integer
        % -----------------------------
        raw_product = 0;
        for k = 1:WordOut
            raw_product = raw_product + sum_row(k) * 2^(k-1);
        end

        % -----------------------------
        % Apply sign
        % -----------------------------
        if result_sign
            raw_product = 2^WordOut - raw_product;
        end
        raw_product = mod(raw_product, 2^WordOut);

        % -----------------------------
        % Store fixed-point output
        % -----------------------------
        temp_unsigned = fi(raw_product, numerictype(0, WordOut, 0));
        out_fmt = numerictype(1, WordOut, FracOut);
        Y_hw(n) = reinterpretcast(temp_unsigned, out_fmt);
    end
end
