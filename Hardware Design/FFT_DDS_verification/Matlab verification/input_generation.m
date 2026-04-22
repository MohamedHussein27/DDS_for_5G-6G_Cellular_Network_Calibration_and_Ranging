%% Prep FFT Input from DDS RTL Output
clear; close all; clc;

N = 4096;
[script_dir, ~, ~] = fileparts(mfilename('fullpath'));
if isempty(script_dir), script_dir = pwd; end

% 1. Read the DDS RTL Output safely
% Assuming the file is named 'rtl_output.txt'
filename = fullfile(script_dir, 'rtl_output.txt');
if ~exist(filename, 'file')
    error('Cannot find rtl_output.txt');
end

% Read the file (skipping any non-numeric text if present)
fid = fopen(filename, 'r');
raw_data = fscanf(fid, '%d');
fclose(fid);

% 2. Extract N points and format for FFT
if length(raw_data) < N
    warning('File has fewer than %d points. Zero-padding.', N);
    input_re = [raw_data; zeros(N - length(raw_data), 1)];
else
    input_re = raw_data(1:N);
end

% Set imaginary part to 0 since DDS gave us a real sine wave
input_im = zeros(N, 1);
input_fixed_complex = complex(input_re, input_im);

% 3. Save for verification later
save(fullfile(script_dir, 'fft_ref_input.mat'), 'input_fixed_complex');

% 4. Write input files for the RTL FFT testbench
fid_re = fopen(fullfile(script_dir, 'input_real.txt'), 'w');
fid_im = fopen(fullfile(script_dir, 'input_imag.txt'), 'w');
for i = 1:N
    fprintf(fid_re, '%d\n', input_re(i));
    fprintf(fid_im, '%d\n', input_im(i));
end
fclose(fid_re); fclose(fid_im);

fprintf('RTL Inputs generated successfully!\n');
fprintf('Peak DDS Amplitude: %d\n', max(abs(input_re)));