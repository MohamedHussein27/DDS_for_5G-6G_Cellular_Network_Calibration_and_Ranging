%% Radix-2² DIF FFT: Fixed-Point Workflow & SQNR Analysis
% 6G ISAC System | NFFT = 4096
clear; close all; clc;

%% 1. Parameters
NFFT = 4096;
num_seeds = 50;                  
WL = 24;                         % Word Length 
FL = 13;                         % Fractional Length
do_instrumentation = true;       % Set TRUE to update range proposals

%% 2. Generate ONE Representative Signal (Required for Instrumentation)
% We need this logic here so we don't feed "zeros" to the tool.
rng(42);
qamOrder = 256;
numBits = log2(qamOrder);
bits = randi([0 1], NFFT*numBits, 1);
X_qam = qammod(bits, qamOrder, 'InputType', 'bit', 'UnitAveragePower', true); %true make normalization on the power
x_float = ifft(X_qam, NFFT);

% Backoff & Scaling
rms_val = rms(x_float); 
target_rms = 0.25; 
scaling_factor = target_rms / rms_val;
x_float = x_float * scaling_factor;

% ADC Clipping
x_real = max(min(real(x_float), 0.999), -0.999);
x_imag = max(min(imag(x_float), 0.999), -0.999);
x_sample = complex(x_real, x_imag); % <--- USE THIS REAL DATA

%% 3. Instrumentation Step (Step 2 & 3)
if do_instrumentation
    fprintf('----------------------------------------\n');
    fprintf('Step 2: Running Instrumentation...\n');
    
    % A) Build the instrumented MEX file
    % We use x_sample in the args so it knows the size and type
    buildInstrumentedMex radix22_dif_fft -o radix22_dif_fft_mex ...
        -args {x_sample} -histogram;
    
    % B) RUN the instrumented MEX file with REAL DATA
    % This is the fix: We pass x_sample instead of dummy_input
    fprintf('Executing compiled code with REAL data to gather stats...\n');
    radix22_dif_fft_mex(x_sample);
    
    % C) SHOW the results
    fprintf('Opening Instrumentation Report Window...\n');
    showInstrumentationResults radix22_dif_fft_mex;
    
    fprintf('Instrumentation Done. Proceeding to SQNR Analysis...\n');
    fprintf('----------------------------------------\n');
end

%% 4. Main SQNR Analysis Loop (Step 4 & Verification)
SQNR_results = zeros(num_seeds, 1);

fprintf('========================================\n');
fprintf('Radix-2² DIF FFT Fixed-Point Analysis\n');
fprintf('========================================\n');
fprintf('NFFT: %d | Seeds: %d\n', NFFT, num_seeds);
fprintf('Fixed-Point Mode: WL=%d, FL=%d\n', WL, FL);
fprintf('----------------------------------------\n');

% Load Type Table
types = get_fft_types('Fixed', WL, FL);
F = fimath('RoundingMethod', 'Nearest', 'OverflowAction', 'Saturate', ...
           'ProductMode', 'FullPrecision', 'SumMode', 'FullPrecision');

for seed = 1:num_seeds
    rng(seed);
    
    % Generate Random Signal (Unique per seed)
    bits = randi([0 1], NFFT*numBits, 1);
    X_qam = qammod(bits, qamOrder, 'InputType', 'bit', 'UnitAveragePower', true);
    x_curr = ifft(X_qam, NFFT);
    
    % Scaling & Clipping
    rms_val = rms(x_curr); 
    scaling_factor = 0.25 / rms_val;
    x_curr = x_curr * scaling_factor;
    
    x_real = max(min(real(x_curr), 0.999), -0.999);
    x_imag = max(min(imag(x_curr), 0.999), -0.999);
    x_curr = complex(x_real, x_imag);
    
    % --- Fixed Point Execution ---
    x_fixed = fi(x_curr, types.Input, F);
    X_ref = fft(x_curr);
    X_dut_fixed = radix22_dif_fft(x_fixed);
    
    % --- SQNR ---
    SQNR_results(seed) = calculate_SQNR(X_ref, double(X_dut_fixed));
    
    if mod(seed, 10) == 0
        fprintf('Completed seed %d/%d, SQNR = %.2f dB\n', seed, num_seeds, SQNR_results(seed));
    end
end

%% 5. Statistical Analysis
fprintf('\n========================================\n');
fprintf('Final Results\n');
fprintf('Mean SQNR:   %.2f dB\n', mean(SQNR_results));
fprintf('Min SQNR:    %.2f dB\n', min(SQNR_results));

%% 6. Plotting Results
figure('Position', [50, 50, 1400, 900]);

% Subplot 1: SQNR vs Seed
subplot(2, 3, 1);
plot(1:num_seeds, SQNR_results, 'b-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'b');
grid on; title(sprintf('SQNR vs Seed (WL=%d, FL=%d)', WL, FL));
yline(mean(SQNR_results), 'r--', 'LineWidth', 2);

% Subplot 2: Histogram
subplot(2, 3, 2);
histogram(SQNR_results, 25, 'FaceColor', [0.3, 0.6, 0.9]);
grid on; title('SQNR Distribution');

% Subplot 3: Magnitude Comparison
subplot(2, 3, 3);
freq_axis = (0:NFFT-1) / NFFT;
plot(freq_axis(1:NFFT/2), 20*log10(abs(X_ref(1:NFFT/2))), 'b-', 'LineWidth', 2);
hold on;
plot(freq_axis(1:NFFT/2), 20*log10(abs(double(X_dut_fixed(1:NFFT/2)))), 'r--', 'LineWidth', 1.5);
grid on; title('Magnitude: Float vs Fixed'); legend('Ref', 'Fixed');

% Subplot 4: Error Magnitude
subplot(2, 3, 4);
error_sig = X_ref - double(X_dut_fixed);
semilogy(freq_axis(1:NFFT/2), abs(error_sig(1:NFFT/2)), 'k-');
grid on; title('Error Magnitude');

% Subplot 5: Phase
subplot(2, 3, 5);
plot(freq_axis(1:NFFT/2), angle(X_ref(1:NFFT/2)), 'b.', 'MarkerSize', 3);
hold on;
plot(freq_axis(1:NFFT/2), angle(double(X_dut_fixed(1:NFFT/2))), 'r.', 'MarkerSize', 3);
grid on; title('Phase');

% Subplot 6: Relative Error
subplot(2, 3, 6);
rel_error = abs(error_sig) ./ (abs(X_ref) + eps);
semilogy(freq_axis(1:NFFT/2), rel_error(1:NFFT/2), 'm-');
grid on; title('Relative Error');

sgtitle(['Fixed-Point Analysis (WL=' num2str(WL) ', FL=' num2str(FL) ')'], 'FontSize', 15);

function sqnr_db = calculate_SQNR(signal_ref, signal_test)
    signal_power = sum(abs(signal_ref).^2);
    noise_power = sum(abs(signal_ref - signal_test).^2);
    if noise_power == 0
        sqnr_db = Inf;
    else
        sqnr_db = 10 * log10(signal_power / noise_power);
    end
end