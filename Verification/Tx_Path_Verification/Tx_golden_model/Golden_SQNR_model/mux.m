%% --- 1. Parameters ---
N = 4096;
N_OFDM = 2048;
N_ZEROS = 381;
N_RADAR = 1667;

%% --- 2. Load Input Data from Text Files ---
disp('Loading OFDM and FFT text files...');

% Load OFDM Symbols
ofdm_re = load('ofdm_data_float_re.txt');
ofdm_im = load('ofdm_data_float_im.txt');
ofdm_complex = ofdm_re + 1j * ofdm_im;

% Load FFT Output (Chirp)
fft_re = load('matlab_fft_re_float.txt');
fft_im = load('matlab_fft_im_float.txt');
fft_complex = fft_re + 1j * fft_im;

%% --- 3. The MUX (Hardware Assignment) ---
disp('Constructing the MUX frame...');

% Initialize the 4096-point frame with zeros 
% (This perfectly handles the N_ZEROS guardband at bins 2049-2429)
mux_frame = zeros(N, 1); 

% 1. Assign OFDM Symbols (MATLAB indices 1 to 2048)
mux_frame(1 : N_OFDM) = ofdm_complex(1 : N_OFDM);

% 2. Assign Radar/Chirp (MATLAB indices 2430 to 4096)
% CRITICAL: Divide by 128 to match the RTL >> 7 arithmetic shift!
radar_start_idx = N_OFDM + N_ZEROS + 1; % Index 2430
mux_frame(radar_start_idx : end) = fft_complex(1 : N_RADAR) / 128.0;

%% --- 4. Export MUX Output for the IFFT Stage ---
disp('Exporting MUX results...');

mux_re_export = real(mux_frame);
mux_im_export = imag(mux_frame);

writematrix(mux_re_export, 'matlab_mux_re_float.txt');
writematrix(mux_im_export, 'matlab_mux_im_float.txt');

disp('Success! MUX Frame saved as floating-point text files.');