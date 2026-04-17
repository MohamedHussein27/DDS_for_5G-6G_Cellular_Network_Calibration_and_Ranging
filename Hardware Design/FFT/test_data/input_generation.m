%% Generate Stimulus Data for Verilog FFT Testbench (Inputs Only)
clear; close all; clc;

% Get the exact folder where THIS script is currently saved
[script_dir, ~, ~] = fileparts(mfilename('fullpath'));
if isempty(script_dir)
    script_dir = pwd; % Fallback to current directory
end
fprintf('Files will be saved directly to: %s\n', script_dir);

%% 1. Generate Inputs
fprintf('Generating 256-QAM signal...\n');
N = 4096; 
M = 256;
scale = 16384; % Q2.14 scaling

rng(42); % Fixed seed so your waveforms are consistent across test runs
data_integers = randi([0, M-1], N, 1);

% Robust QAM generation
if exist('qammod', 'file')
    ofdm_symbols = qammod(data_integers, M, 'UnitAveragePower', true);
else
    [I, Q] = meshgrid(-15:2:15, -15:2:15);
    const = (I + 1j*Q); 
    const = const(:) / sqrt(mean(abs(const(:)).^2));
    ofdm_symbols = const(data_integers + 1);
end

% Save the ideal QAM symbols for later MATLAB EVM/BER analysis
mat_file = fullfile(script_dir, 'qam_golden_16b.mat');
save(mat_file, 'ofdm_symbols');
fprintf('Saved QAM symbols to: %s\n', 'qam_golden_16b.mat');

%% 2. Transform and Quantize
x_time = ifft(ofdm_symbols, N);
x_hw_in = round(x_time * scale);

% Separate real and imaginary for hardware text files
in_real_int = real(x_hw_in);
in_imag_int = imag(x_hw_in);

% Clamp to strict 16-bit signed integer limits to prevent Verilog overflow
in_real_int = max(min(in_real_int, 32767), -32768);
in_imag_int = max(min(in_imag_int, 32767), -32768);

%% 3. Write to Text Files (Decimal format for Verilog $fscanf)
path_in_re = fullfile(script_dir, 'input_real.txt');
path_in_im = fullfile(script_dir, 'input_imag.txt');

fid_in_re = fopen(path_in_re, 'w');
fid_in_im = fopen(path_in_im, 'w');

if fid_in_re == -1 || fid_in_im == -1
    error('Failed to open text files. Check folder permissions.');
end

% Write exact integer values line-by-line
for i = 1:N
    fprintf(fid_in_re, '%d\n', in_real_int(i));
    fprintf(fid_in_im, '%d\n', in_imag_int(i));
end

fclose(fid_in_re); 
fclose(fid_in_im);

fprintf('Done! Input text files are ready for your QuestaSim testbench.\n');