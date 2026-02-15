function [stream_pos, stream_neg, valid_out] = demux(Y_in, N)

    % Initialize Output Streams (Hardware FIFOs)
    stream_pos = complex(zeros(N/2, 1));
    stream_neg = complex(zeros(N/2, 1));
    valid_out  = false(N, 1);
    
    % Internal Pointers (Write Addresses for the FIFOs)
    ptr_pos = 1;
    ptr_neg = 1;

    % --- HARDWARE PIPELINE SIMULATION ---
    for k = 1:N
        % Read input sample
        sample_in = Y_in(k);
        
        if k <= N/2
            % PHASE 1: POSITIVE FREQUENCIES (0 to +Fs/2)
            % Hardware: Route switch to Output A
            stream_pos(ptr_pos) = sample_in;
            ptr_pos = ptr_pos + 1;
        else
            % PHASE 2: NEGATIVE FREQUENCIES (-Fs/2 to 0)
            % Hardware: Route switch to Output B
            stream_neg(ptr_neg) = sample_in;
            ptr_neg = ptr_neg + 1;
        end
        
        valid_out(k) = true;
    end
end