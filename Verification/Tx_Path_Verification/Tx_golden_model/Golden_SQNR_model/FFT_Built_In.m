% --- MATLAB FFT Floating-Point Reference ---
clear; clc;

% 1. Load the floating-point DDS output you generated earlier
% (Or you can just run this right after generating 's' in your previous script)
if isfile('matlab_dds_float_out.txt')
    s = load('matlab_dds_float_out.txt');
else
    error('Could not find matlab_dds_float_out.txt. Run the DDS script first!');
end

% 2. Perform the built-in FFT
% s is a 1D vector (e.g., 4096 samples)
fft_out = fft(s);

% 3. Extract Real and Imaginary parts
fft_re = real(fft_out);
fft_im = imag(fft_out);

% 4. Save to text files for Python to read
writematrix(fft_re(:), 'matlab_fft_re_float.txt');
writematrix(fft_im(:), 'matlab_fft_im_float.txt');

disp('Successfully exported FFT Real and Imaginary floating-point vectors.');