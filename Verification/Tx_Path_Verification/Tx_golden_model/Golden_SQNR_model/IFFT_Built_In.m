%% --- 1. Parameters ---
N = 4096;

%% --- 2. Load Input Data from MUX ---
disp('Loading floating-point MUX output files...');

% Load the MUX output generated in the previous step
mux_re = load('matlab_mux_re_float.txt');
mux_im = load('matlab_mux_im_float.txt');

% Combine into a complex vector
mux_complex = mux_re + 1j * mux_im;

%% --- 3. Perform the IFFT ---
disp('Executing IFFT...');

% Standard MATLAB IFFT
ifft_out = ifft(mux_complex);

% CRITICAL HARDWARE SCALING: 
% MATLAB automatically divides by N. Hardware typically does not. 
% Multiply by N to match the unscaled integer growth of the RTL.
ifft_unscaled = ifft_out * N;

%% --- 4. Export the Final Floating-Point Outputs ---
disp('Exporting final IFFT results...');

ifft_re_export = real(ifft_unscaled);
ifft_im_export = imag(ifft_unscaled);

% Export to the final text files
writematrix(ifft_re_export, 'matlab_ifft_re_float.txt');
writematrix(ifft_im_export, 'matlab_ifft_im_float.txt');

disp('Success! Final TX datapath outputs saved as floating-point text files.');