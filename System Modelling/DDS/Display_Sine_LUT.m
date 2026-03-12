% =========================================================
% QUARTER-WAVE SINE LUT VIEWER (N/4 + 1 Rule)
% =========================================================
clc; clear; close all;

LUT_bits = 8;               % Total full-wave address bits
Fractional_Bits = 6;        % Example Q-format (Q2.6) used in your test
Word_Length = 8;            % Total bits per word (1 Sign + 1 Int + 6 Frac)

% 1. Calculate Quarter-Wave Sizes
Q_LUT_bits = LUT_bits - 2; 
Q_LUT_size = 2^Q_LUT_bits;  % 64
addr_ROM = 0 : Q_LUT_size;  % Indices 0 to 64 (65 points total)

% 2. Generate the True Mathematical Values
phase_rads = 2*pi*addr_ROM / (2^LUT_bits);
sine_LUT_double = sin(phase_rads);

% 3. Generate Fixed-Point Values (Hardware Representation)
% Using a signed fixed-point type with specified word/fractional lengths
T = numerictype(1, Word_Length, Fractional_Bits);
F = fimath('RoundingMethod', 'Nearest', 'OverflowAction', 'Saturate');
sine_LUT_fixed = fi(sine_LUT_double, T, F);

% =========================================================
% PRINT THE ROM CONTENTS
% =========================================================
fprintf('--- QUARTER-WAVE SINE LUT CONTENTS (LUT_bits = %d) ---\n', LUT_bits);
fprintf('Total Entries: %d (Indices 0 to %d)\n', length(addr_ROM), Q_LUT_size);
fprintf('Format: Q%d.%d\n\n', Word_Length - Fractional_Bits - 1, Fractional_Bits);

fprintf(' Index | Phase (deg) | Ideal Value | Fixed-Point | HEX Value\n');
fprintf('--------------------------------------------------------------\n');

for k = 1:length(addr_ROM)
    idx   = addr_ROM(k);
    deg   = phase_rads(k) * (180/pi);
    ideal = sine_LUT_double(k);
    fixed = double(sine_LUT_fixed(k));
    hex_val = hex(sine_LUT_fixed(k)); % <--- Renamed variable here
    
    fprintf('  %3d  |   %6.2f    |  %8.5f   |  %8.5f   |   %s\n', ...
            idx, deg, ideal, fixed, hex_val); % <--- And here
end

% =========================================================
% VISUALIZE THE MEMORY
% =========================================================
figure('Name', 'ROM Contents', 'Color', 'w', 'Position', [200, 200, 700, 400]);
stem(addr_ROM, sine_LUT_double, 'k', 'LineWidth', 1.5, 'MarkerFaceColor', 'k'); hold on;
plot(addr_ROM, double(sine_LUT_fixed), 'r--', 'LineWidth', 1.5);
grid on;
title('Quarter-Wave LUT Values Stored in Memory');
xlabel('ROM Address (Index)');
ylabel('Amplitude');
legend('Ideal Mathematical Value', sprintf('Fixed-Point (Q2.%d)', Fractional_Bits), 'Location', 'SouthEast');
xlim([0 Q_LUT_size]);