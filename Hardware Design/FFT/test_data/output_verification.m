%% Hardware FFT Verification and Analysis Dashboard
% Reads the RTL output, aligns it with the Golden QAM reference, 
% calculates EVM, and plots the constellation and subcarrier errors.

clear; close all; clc;

%% 1. Parameters & Setup
N = 4096;
scale = 16384; % Q2.14 Scaling Factor (2^14)

% Get script directory
[script_dir, ~, ~] = fileparts(mfilename('fullpath'));
if isempty(script_dir), script_dir = pwd; end

%% 2. Load Golden Reference
mat_file = fullfile(script_dir, 'qam_golden_16b.mat');
if ~exist(mat_file, 'file')
    error('Could not find qam_golden_16b.mat. Run the input generation script first.');
end
load(mat_file, 'ofdm_symbols');

%% 3. Read Hardware Output (Decimal Text Files)
path_out_re = fullfile(script_dir, 'output_real.txt');
path_out_im = fullfile(script_dir, 'output_imag.txt');

fid_re = fopen(path_out_re, 'r');
fid_im = fopen(path_out_im, 'r');

if fid_re == -1 || fid_im == -1
    error('Could not find output_real.txt or output_imag.txt. Run QuestaSim first.');
end

% Read signed decimal integers directly
hw_re_dec = fscanf(fid_re, '%d');
hw_im_dec = fscanf(fid_im, '%d');
fclose(fid_re); fclose(fid_im);

% Safety check to ensure exactly 4096 samples were captured
if length(hw_re_dec) ~= N || length(hw_im_dec) ~= N
    warning('Hardware output files do not contain exactly %d lines! Check testbench.', N);
    % Truncate to N just in case to prevent matrix dimension errors
    hw_re_dec = hw_re_dec(1:min(N, end));
    hw_im_dec = hw_im_dec(1:min(N, end));
end

%% 4. Data Recovery & Alignment
% Combine into complex integers
hw_complex_int = complex(hw_re_dec, hw_im_dec);

% Bit-reverse the hardware output (DIF FFT outputs bit-reversed order)
hw_complex_int = bitrevorder(hw_complex_int);

% Convert back to floating point
% NOTE: Because the hardware divides by 2 at each of the 12 stages, the 
% overall gain of the hardware pipeline is (1 / 4096). We apply a dynamic 
% gain alignment here to perfectly map the hardware amplitude to the golden 
% constellation space, isolating pure quantization noise.
hw_float = hw_complex_int / scale;
gain_correction = rms(abs(ofdm_symbols(:))) / rms(abs(hw_float(:)));
X_hw_recovered = hw_float * gain_correction;

%% 5. Calculate Error Vector Magnitude (EVM)
error_vector = X_hw_recovered(:) - ofdm_symbols(:);
error_magnitude = abs(error_vector);

evm_rms_percent = (rms(error_magnitude) / rms(abs(ofdm_symbols(:)))) * 100;
fprintf('Hardware EVM: %.4f %%\n', evm_rms_percent);

%% 6. Plotting: Constellation Diagram
figure('Name', '16-bit Hardware QAM Recovery', 'Position', [100, 100, 600, 600]);
plot(real(ofdm_symbols), imag(ofdm_symbols), 'k+', 'MarkerSize', 8); hold on;
plot(real(X_hw_recovered), imag(X_hw_recovered), 'r.', 'MarkerSize', 8);
title(sprintf('16-bit Recovered 256-QAM\nHardware EVM: %.3f%%', evm_rms_percent));
legend('MATLAB Golden', 'RTL Output', 'Location', 'best');
grid on; 
axis([-1.5 1.5 -1.5 1.5]);

%% 7. Plotting: Subcarriers and Error
figure('Name', 'Hardware vs Reference: Subcarrier Analysis', 'Position', [750, 100, 1000, 800]);

% --- Plot 1: Real Part Comparison ---
subplot(3, 1, 1);
plot(real(ofdm_symbols), 'b', 'LineWidth', 1.5); hold on;
plot(real(X_hw_recovered), 'r--', 'LineWidth', 1);
title('Subcarriers: Real Part (In-Phase)');
ylabel('Amplitude'); xlabel('Subcarrier Index');
legend('MATLAB Golden', 'Verilog RTL', 'Location', 'northeast');
grid on; xlim([1 N]);

% --- Plot 2: Imaginary Part Comparison ---
subplot(3, 1, 2);
plot(imag(ofdm_symbols), 'b', 'LineWidth', 1.5); hold on;
plot(imag(X_hw_recovered), 'r--', 'LineWidth', 1);
title('Subcarriers: Imaginary Part (Quadrature)');
ylabel('Amplitude'); xlabel('Subcarrier Index');
legend('MATLAB Golden', 'Verilog RTL', 'Location', 'northeast');
grid on; xlim([1 N]);

% --- Plot 3: Absolute Error Magnitude ---
subplot(3, 1, 3);
plot(error_magnitude, 'k', 'LineWidth', 1.2);
title(sprintf('Absolute Quantization Error (|Hardware - Golden|) - Max Error: %.4f', max(error_magnitude)));
ylabel('Error Amplitude'); xlabel('Subcarrier Index');
grid on; xlim([1 N]);
ylim([0, max(error_magnitude) * 1.5]); 

sgtitle(sprintf('End-to-End FFT RTL Validation (EVM: %.3f%%)', evm_rms_percent), 'FontSize', 14, 'FontWeight', 'bold');