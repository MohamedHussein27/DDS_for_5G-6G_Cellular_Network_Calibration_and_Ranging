% Ideal Mathematical Chirp vs RTL Comparison (Absolute Bit-True ROM Fetch)
clear; clc;

% --- Hardware Parameters ---
Fs = 491.52e6;             % System clock
N_cycles = 4096;           % Number of simulated cycles
T_dur = N_cycles / Fs;     % Total duration in seconds
FTW_step = 426666;         % Hardware FTW step 
FTW_start = 0;             % Hardware FTW start 

% Calculate the exact physical bandwidth reached by the hardware
max_ftw = FTW_start + (FTW_step * N_cycles);
B = ((max_ftw-FTW_start) * Fs) / (2^32); 

% --- 1. Load RTL Output & ROM Data ---
try
    rtl_data = load('rtl_output.txt');
catch
    error('Could not find rtl_output.txt.');
end

try
    % Read the exact hex strings from the Verilog memory file
    fid = fopen('memfile.mem', 'r');
    mem_cells = textscan(fid, '%x'); % Automatically converts Hex to Decimal
    fclose(fid);
    rom_data = mem_cells{1};
catch
    error('Could not load memfile.mem.');
end

% --- 2. Align Pipeline Latency ---
latency = 0; 
rtl_aligned = rtl_data(1+latency : N_cycles+latency)';

% --- 3. Generate Ideal Discrete-Time Model (BIT-TRUE DIRECT FETCH) ---
n = 0:(N_cycles-1); 
t = n / Fs;

% Discrete Phase Accumulation (32-bit wrapping)
ideal_discrete_phase_word = mod((FTW_start .* n) + (FTW_step * (n .* (n - 1)) / 2), 2^32);

% ---- EMULATE PHASE_ACC.v ----
% Truncate to ADDRESS_WIDTH = 16 bits (drop the bottom 16 bits)
truncated_phase_16b = floor(ideal_discrete_phase_word / 65536); 

% ---- EMULATE Quadrant Mapper & LUT ----
ideal_out = zeros(size(truncated_phase_16b));

for k = 1:length(truncated_phase_16b)
    current_phase = truncated_phase_16b(k);
    
    % Extract top 2 bits for quadrant, lower 14 bits for the mapped address
    quadrant = floor(current_phase / 16384); 
    addr = mod(current_phase, 16384); 
    
    % Quadrant Mapping Logic (Emulating first_quad_address.v)
    if quadrant == 0     % 1st quad
        mapped_addr = addr;
        neg_flag = 0;
    elseif quadrant == 1 % 2nd quad
        mapped_addr = 16383 - addr;
        neg_flag = 0;
    elseif quadrant == 2 % 3rd quad
        mapped_addr = addr; 
        neg_flag = 1;
    else                 % 4th quad
        mapped_addr = 16383 - addr;
        neg_flag = 1;
    end
    
    % ---- EMULATE LUT.v (Direct ROM Read) ----
    % MATLAB arrays are 1-indexed, so we add 1 to the mapped_addr
    lut_amplitude = double(rom_data(mapped_addr + 1));
    
    % ---- EMULATE negative_mux.v ----
    if neg_flag == 1
        ideal_out(k) = -lut_amplitude;
    else
        ideal_out(k) = lut_amplitude;
    end
end

% --- 4. Plotting Time Domain ---
figure('Name', 'DDS Output: RTL vs Ideal');

% Subplot 1: The Waveforms
subplot(2,1,1);
plot(t * 1e6, ideal_out, 'b', 'LineWidth', 1.5); hold on;
plot(t * 1e6, rtl_aligned, 'r--', 'LineWidth', 1.5);
title('LFM Chirp: RTL vs Bit-True Hardware Emulator');
xlabel('Time (\mus)');
ylabel('Amplitude');
legend('Ideal Hardware Emulator', 'Verilog RTL (Aligned)');
grid on;

% Subplot 2: The Error (Should be entirely flat at 0)
subplot(2,1,2);
noise_signal = ideal_out - rtl_aligned;
plot(t * 1e6, noise_signal, 'k');
title('Hardware Error (Verilog vs MATLAB)');
xlabel('Time (\mus)');
ylabel('Amplitude Error (LSBs)');
ylim([-5 5]); % Zoomed in to prove absolute zero
grid on;

% --- 5. Signal Quality Metrics ---
signal_power = rms(ideal_out)^2;
noise_power = rms(noise_signal)^2;

% Calculate SQNR safely to output 'Inf'
if noise_power == 0
    sqnr_db = Inf;
else
    sqnr_db = 10 * log10(double(signal_power) / double(noise_power));
end

fprintf('Calculated Physical Bandwidth: %.2f MHz\n', B / 1e6);
fprintf('Measured Hardware SQNR: %.2f dB\n', sqnr_db);

% --- 6. Frequency Spectrum (FFT) ---
Y = fft(rtl_aligned);
P2 = abs(Y / N_cycles);
P1 = P2(1:N_cycles/2+1);
P1(2:end-1) = 2 * P1(2:end-1); 

max_amp = 127;
P1_dBFS = 20 * log10(P1 / max_amp + eps); 

f = Fs * (0:(N_cycles/2)) / N_cycles;

figure('Name', 'RTL Output Frequency Spectrum');
plot(f / 1e6, P1_dBFS, 'b', 'LineWidth', 1.2); hold on;
title('Single-Sided Amplitude Spectrum of RTL Chirp');
xlabel('Frequency (MHz)');
ylabel('Magnitude (dBFS)');
ylim([-100 0]); 
grid on;

xline(B / 1e6, 'r--', 'Calculated Max Bandwidth', 'LabelHorizontalAlignment', 'left', 'LineWidth', 1.5);