% % Generate Stimulus Data (Raw QAM to Hardware)
% clear; close all; clc;
% 
% [script_dir, ~, ~] = fileparts(mfilename('fullpath'));
% if isempty(script_dir), script_dir = pwd; end
% 
% % 1. Generate QAM Inputs (No FFT/IFFT applied here)
% N = 4096;
% M = 256;
% scale_data = 256; % Q8.8 scaling
% 
% rng(42); 
% data_integers = randi([0, M-1], N, 1);
% 
% if exist('qammod', 'file')
%     ofdm_symbols = qammod(data_integers, M, 'UnitAveragePower', true);
% else
%     [I, Q] = meshgrid(-15:2:15, -15:2:15);
%     const = (I + 1j*Q);
%     const = const(:) / sqrt(mean(abs(const(:)).^2));
%     ofdm_symbols = const(data_integers + 1);
% end
% 
% Save the exact frequency-domain symbols
% mat_file = fullfile(script_dir, 'qam_golden_q88.mat');
% save(mat_file, 'ofdm_symbols');
% 
% % 2. Scale and Quantize
% We divide by 4 to give the unscaled Verilog pipeline the dynamic 
% headroom it needs to grow the signal without clipping the 32767 limit.
% x_hw_in = round((ofdm_symbols / 4) * scale_data); 
% 
% in_real_int = real(x_hw_in);
% in_imag_int = imag(x_hw_in);
% 
% Clamp to 16-bit signed integer limits
% in_real_int = max(min(in_real_int, 32767), -32768);
% in_imag_int = max(min(in_imag_int, 32767), -32768);
% 
% % 3. Write to Text Files
% path_in_re = fullfile(script_dir, 'input_real.txt');
% path_in_im = fullfile(script_dir, 'input_imag.txt');
% 
% fid_in_re = fopen(path_in_re, 'w');
% fid_in_im = fopen(path_in_im, 'w');
% 
% for i = 1:N
%     fprintf(fid_in_re, '%d\n', in_real_int(i));
%     fprintf(fid_in_im, '%d\n', in_imag_int(i));
% end
% 
% fclose(fid_in_re); fclose(fid_in_im);
% fprintf('Done! Raw QAM input text files generated.\n');

%% Generate Stimulus Data (Fixed Scaling)
clear; close all; clc;

[script_dir, ~, ~] = fileparts(mfilename('fullpath'));
if isempty(script_dir), script_dir = pwd; end

%% 1. Generate QAM Inputs
N = 4096;
M = 256;

rng(42); 
data_integers = randi([0, M-1], N, 1);

if exist('qammod', 'file')
    ofdm_symbols = qammod(data_integers, M, 'UnitAveragePower', true);
else
    [I, Q] = meshgrid(-15:2:15, -15:2:15);
    const = (I + 1j*Q);
    const = const(:) / sqrt(mean(abs(const(:)).^2));
    ofdm_symbols = const(data_integers + 1);
end

%% 2. Optimized  Scaling
optimized_scale = 128;

% Apply optimized scale (Removed the * 256 bug)
x_hw_in = round(ofdm_symbols * optimized_scale); 

% Save symbols AND the optimized scale factor for bit-true verification
mat_file = fullfile(script_dir, 'qam_golden_q88.mat');
save(mat_file, 'ofdm_symbols', 'optimized_scale');

in_real_int = real(x_hw_in);
in_imag_int = imag(x_hw_in);

% Clamp to 16-bit signed integer limits for hardware safety
in_real_int = max(min(in_real_int, 32767), -32768);
in_imag_int = max(min(in_imag_int, 32767), -32768);

%% 3. Write to Text Files
path_in_re = fullfile(script_dir, 'input_real.txt');
path_in_im = fullfile(script_dir, 'input_imag.txt');

fid_in_re = fopen(path_in_re, 'w');
fid_in_im = fopen(path_in_im, 'w');

for i = 1:N
    fprintf(fid_in_re, '%d\n', in_real_int(i));
    fprintf(fid_in_im, '%d\n', in_imag_int(i));
end

fclose(fid_in_re); 
fclose(fid_in_im);

fprintf('Done! Optimized input files generated with scale: %.4f\n', optimized_scale);