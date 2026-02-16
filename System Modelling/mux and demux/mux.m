function [X_out, valid_out] = mux(stream_pos, stream_neg, N)
% HARDWARE_MUX_NOSHIFT Combines streams into Standard FFT Order
%
% Inputs:
%   stream_pos : Data for Positive Frequencies (0 to +Fs/2) -> Indices 1 to N/2
%   stream_neg : Data for Negative Frequencies (-Fs/2 to 0) -> Indices N/2+1 to N
%   N          : Total FFT size
%
% Outputs:
%   X_out      : Combined spectrum ready for direct IFFT (Standard Order)
%   valid_out  : "Data Valid" strobe

    X_out = complex(zeros(N, 1));
    valid_out = false(N, 1);
    
    % Internal Pointers (Hardware Read Addresses)
    ptr_pos = 1;
    ptr_neg = 1;

    % --- HARDWARE PIPELINE SIMULATION ---
    % The counter 'k' represents the memory address writing to the IFFT buffer
    for k = 1:N
        if k <= N/2
            % PHASE 1: POSITIVE FREQUENCIES (0 to +Fs/2)
            % Hardware: Select Input A (Positive Stream)
            X_out(k) = stream_pos(ptr_pos);
            ptr_pos = ptr_pos + 1;
        else
            % PHASE 2: NEGATIVE FREQUENCIES (-Fs/2 to 0)
            % Hardware: Select Input B (Negative Stream)
            X_out(k) = stream_neg(ptr_neg);
            ptr_neg = ptr_neg + 1;
        end
        
        valid_out(k) = true;
    end
end