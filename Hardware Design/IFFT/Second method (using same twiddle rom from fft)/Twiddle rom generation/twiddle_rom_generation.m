%% Corrected Twiddle ROM Generation (Direct Folder Output) - Q2.14
clear; clc;

N_total = 4096;
WL = 16;
scale_factor = 16384; % Q2.14 scaling (16384)

% Get the exact folder where THIS script is currently saved
[script_dir, ~, ~] = fileparts(mfilename('fullpath'));
if isempty(script_dir)
    script_dir = pwd; % Fallback to current directory if not saved
end

fprintf('Generating Q2.14 .mem files directly in: %s\n', script_dir);

for stg = 1:12
    % MAP TO VERILOG STAGES:
    % Stage 1  -> Depth 2048
    % Stage 12 -> Depth 1
    depth = 2^(12 - stg); 
    
    % DIF FFT Phase Calculation
    N_stage = depth * 2;
    k = 0:(depth-1);
    W = exp(-1j * 2 * pi * k / N_stage);
    
    % Fixed-point conversion (Q2.14)
    W_re = round(real(W) * scale_factor);
    W_im = round(imag(W) * scale_factor);
    
    % Saturation limits for 16-bit signed
    W_re = max(min(W_re, 32767), -32768);
    W_im = max(min(W_im, 32767), -32768);
    
    % Convert to 4-digit Hexadecimal strings
    re_hex = dec2hex(mod(W_re, 65536), 4);
    im_hex = dec2hex(mod(W_im, 65536), 4);
    
    % Build file paths in the SAME folder as the script
    f_re_path = fullfile(script_dir, sprintf('twiddle_stage_%d_real.mem', stg));
    f_im_path = fullfile(script_dir, sprintf('twiddle_stage_%d_imag.mem', stg));
    
    % Open and Write Files
    fid_re = fopen(f_re_path, 'w');
    fid_im = fopen(f_im_path, 'w');
    
    if fid_re == -1 || fid_im == -1
        error('Could not open file for writing in: %s', script_dir);
    end
    
    for i = 1:depth
        fprintf(fid_re, '%s\n', re_hex(i,:));
        fprintf(fid_im, '%s\n', im_hex(i,:));
    end
    
    fclose(fid_re); 
    fclose(fid_im);
    fprintf('Stage %2d: Created %4d lines in %s\n', stg, depth, sprintf('twiddle_stage_%d_real.mem', stg));
end

fprintf('\nDone! All 24 Q2.14 files are ready for QuestaSim.\n');