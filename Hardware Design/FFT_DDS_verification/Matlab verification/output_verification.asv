%% Verify RTL FFT Output against Fixed-Point MATLAB Model (No-Shift)
clear; close all; clc;

N = 4096;
[script_dir, ~, ~] = fileparts(mfilename('fullpath'));
if isempty(script_dir), script_dir = pwd; end

% --- 1. Load Reference Input ---
if exist(fullfile(script_dir, 'fft_ref_input.mat'), 'file')
    load(fullfile(script_dir, 'fft_ref_input.mat'), 'input_fixed_complex');
else
    error('Reference input not found. Run prep_fft_input.m first.');
end

% --- 2. Read RTL Output Files ---
fid_re = fopen(fullfile(script_dir, 'output_real.txt'), 'r');
fid_im = fopen(fullfile(script_dir, 'output_imag.txt'), 'r');
if fid_re == -1 || fid_im == -1
    error('RTL output files not found.');
end

hw_re = fscanf(fid_re, '%d');
hw_im = fscanf(fid_im, '%d');
fclose(fid_re); fclose(fid_im);

% --- 3. Format & Bit-Reverse ---
N_actual = min([N, length(hw_re), length(hw_im)]);
rtl_out_raw = complex(hw_re(1:N_actual), hw_im(1:N_actual));
rtl_out_natural = bitrevorder(rtl_out_raw); % DIF hardware outputs bit-reversed

% --- 4. Compute MATLAB Reference (Custom Fixed-Point Q2.14 Model) ---
WL = 16; % 16-bit word length
FL = 14; % Q2.14 format (14 fractional bits for twiddle multiplication)

% THE FIX: Scale the raw integers into the Q2.14 fractional domain 
% so MATLAB doesn't saturate the inputs at 1.9999!
input_q_scaled = double(input_fixed_complex) / (2^FL);

% Call your custom bit-accurate FFT design function
% Note: If your custom function outputs bit-reversed data directly (unlike MATLAB's built-in),
% you will need to apply bitrevorder() to matlab_ref as well to match the plotting.
matlab_ref = radix22_dif_fft_fixed(input_q_scaled, WL, FL);

% THE FIX: Convert the Q2.14 output back to raw integers to compare with RTL
% (We handle whether your function outputs a 'fi' object or a standard array)
if isobject(matlab_ref)
    ref_vec = double(matlab_ref.data(:)) * (2^FL);
else
    ref_vec = double(matlab_ref(:)) * (2^FL);
end

% Round just in case of microscopic floating point errors
ref_vec = round(ref_vec);

% --- 5. Performance Metrics Calculation ---

rtl_vec = double(rtl_out_natural(:));
err_vec = ref_vec - rtl_vec;

% Power calculations (Scalar)
sig_pwr   = real(ref_vec' * ref_vec);
noise_pwr = real(err_vec' * err_vec);

% Calculate SQNR and EVM (catching perfect match for Inf SQNR)
if noise_pwr > 0
    sqnr_val = 10 * log10(sig_pwr / noise_pwr);
    evm_rms = (sqrt(mean(abs(err_vec).^2)) / sqrt(mean(abs(ref_vec).^2))) * 100;
else
    sqnr_val = Inf; % Error is exactly 0
    evm_rms = 0.0000; 
end

% --- 6. Final Report ---
fprintf('\n================================================\n');
fprintf('   BIT-ACCURATE FFT VERIFICATION REPORT         \n');
fprintf('================================================\n');
fprintf(' FFT Points:        %d\n', N);
fprintf(' Reference Model:   radix22_dif_fft_fixed (Q2.14)\n');
fprintf(' SQNR:              %8.2f dB\n', sqnr_val);

% Print success message if SQNR is Infinity
if isinf(sqnr_val)
    fprintf('                    --> SUCCESS: PERFECT BIT-ACCURATE MATCH!\n');
end

fprintf(' EVM (RMS):         %8.4f %%\n', evm_rms);
fprintf(' Max Error (LSB):   %8.2f\n', max(abs(err_vec)));
fprintf('================================================\n');

% --- 7. Visualization ---
figure('Name', 'FFT Performance Analysis', 'Position', [100, 100, 900, 400]);

% Subplot 1: Frequency Spectrum
subplot(1,2,1);
Fs = 491.52; % MHz
freq_axis = (0:N-1) * (Fs/N);
plot(freq_axis, 20*log10(abs(rtl_vec)), 'r', 'LineWidth', 1.5); hold on;
plot(freq_axis, 20*log10(abs(ref_vec)), 'b--', 'LineWidth', 1.5);
grid on; axis tight;
title(['Spectrum (SQNR: ', num2str(sqnr_val, '%.2f'), ' dB)']);
xlabel('Frequency (MHz)');
ylabel('Magnitude (dB)');
legend('RTL Hardware', 'MATLAB Ideal');

% Subplot 2: Complex Scatter Plot (Constellation)
subplot(1,2,2);
plot(real(rtl_vec), imag(rtl_vec), 'r.', 'MarkerSize', 8); hold on;
plot(real(ref_vec), imag(ref_vec), 'b+', 'MarkerSize', 4);
grid on; axis square;
title(['Output Scatter (EVM: ', num2str(evm_rms, '%.3f'), '%)']);
xlabel('Real / In-Phase');
ylabel('Imag / Quadrature');
legend('RTL Hardware', 'MATLAB Ideal');

% --- DEBUG PRINT: Let's look at the actual numbers ---
fprintf('\n--- DEBUG: FIRST 5 SAMPLES ---\n');
fprintf('MATLAB Model:\n');
disp(ref_vec(1:5));
fprintf('RTL Hardware:\n');
disp(rtl_vec(1:5));