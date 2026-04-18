%% Hardware IFFT Verification (Time-Domain Output Validation)
clear; close all; clc;

%% 1. Parameters & Setup
N = 4096;
scale_data = 256; % Q8.8 Scaling

[script_dir, ~, ~] = fileparts(mfilename('fullpath'));
if isempty(script_dir), script_dir = pwd; end

%% 2. Load Golden Reference & Compute MATLAB IFFT
mat_file = fullfile(script_dir, 'qam_golden_q88.mat');
load(mat_file, 'ofdm_symbols');

% 1. Apply the same headroom protection (/4) we used on the hardware input
% 2. Compute the exact mathematical IFFT using MATLAB
matlab_time_domain = ifft(ofdm_symbols / 4, N);

%% 3. Read Hardware Output
path_out_re = fullfile(script_dir, 'output_real.txt');
path_out_im = fullfile(script_dir, 'output_imag.txt');
fid_re = fopen(path_out_re, 'r');
fid_im = fopen(path_out_im, 'r');

hw_re_dec = fscanf(fid_re, '%d');
hw_im_dec = fscanf(fid_im, '%d');
fclose(fid_re); fclose(fid_im);

if length(hw_re_dec) ~= N || length(hw_im_dec) ~= N
    hw_re_dec = hw_re_dec(1:min(N, end));
    hw_im_dec = hw_im_dec(1:min(N, end));
end

%% 4. Data Recovery & Alignment
hw_complex_int = complex(hw_re_dec, hw_im_dec);
hw_complex_int = bitrevorder(hw_complex_int);

% Reverse the Q8.8 scaling
hw_float = hw_complex_int / scale_data; 

% --- ALIGNING MATHEMATICS ---
% MATLAB's ifft() formula inherently divides by N.
% Your unscaled Verilog RTL is a pure sum, which DOES NOT divide by N.
% To make the RTL amplitude match MATLAB perfectly, we divide the RTL by N here.
hw_time_recovered = hw_float / N; 

%% 5. Calculate Error Vector Magnitude (Time-Domain)
error_vector = hw_time_recovered(:) - matlab_time_domain(:);
error_magnitude = abs(error_vector);
evm_rms_percent = (rms(error_magnitude) / rms(abs(matlab_time_domain(:)))) * 100;
fprintf('Hardware Time-Domain EVM: %.4f %%\n', evm_rms_percent);

%% 6. Plotting: Time-Domain "Constellation"
figure('Name', 'Hardware vs MATLAB: Time Domain Scatter', 'Position', [100, 100, 600, 600]);
plot(real(matlab_time_domain), imag(matlab_time_domain), 'k+', 'MarkerSize', 8); hold on;
plot(real(hw_time_recovered), imag(hw_time_recovered), 'r.', 'MarkerSize', 8);
title(sprintf('Time-Domain OFDM Signal (Not a Grid!)\nHardware EVM: %.3f%%', evm_rms_percent));
legend('MATLAB Golden IFFT', 'RTL Output', 'Location', 'best');
grid on; 

%% 7. Plotting: Waveform Alignment and Error
figure('Name', 'Hardware vs Reference: Waveform Analysis', 'Position', [750, 100, 1000, 800]);

subplot(3, 1, 1);
plot(real(matlab_time_domain), 'b', 'LineWidth', 1.5); hold on;
plot(real(hw_time_recovered), 'r--', 'LineWidth', 1);
title('Time-Domain Waveform: Real Part (In-Phase)');
grid on; xlim([1 N]);

subplot(3, 1, 2);
plot(imag(matlab_time_domain), 'b', 'LineWidth', 1.5); hold on;
plot(imag(hw_time_recovered), 'r--', 'LineWidth', 1);
title('Time-Domain Waveform: Imaginary Part (Quadrature)');
grid on; xlim([1 N]);

subplot(3, 1, 3);
plot(error_magnitude, 'k', 'LineWidth', 1.2);
title(sprintf('Absolute Time-Domain Error (|RTL - MATLAB|) - Max Error: %.4f', max(error_magnitude)));
grid on; xlim([1 N]); ylim([0, max(error_magnitude) * 1.5]);