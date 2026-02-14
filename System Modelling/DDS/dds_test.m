clear; clc; close all;
% =========================================================
% 1. DESIGN PARAMETERS
% =========================================================
Fs       = 500e6;        % System clock (Hz)
T_dur    = 4.096e-6;     % Signal duration (s)
Nacc     = 24;           % Phase accumulator bits
LUT_bits = 10;           % LUT address bits
nSeeds   = 60;           % Number of random tests
SQNR_fixed  = zeros(nSeeds, 1);
Ns = round(T_dur * Fs);  % Number of samples
t  = (0:Ns-1) / Fs;      % Time vector

% -----------------------------------------------------
% CONFIGURATION
% -----------------------------------------------------
current_mode = 'fixed'; % Options: 'fixed', 'double', 'single'
fprintf('Testing Mode: %s\n', current_mode);

% Define Q-Formats to sweep (Only used if mode is 'fixed')
if strcmp(current_mode, 'fixed')
    Q_list = [14 12 10 8 6]; 
else
    Q_list = [0]; % Dummy value for non-fixed modes
end

avg_sqnr_results = zeros(length(Q_list), 1);
min_sqnr_results = zeros(length(Q_list), 1);

% Global variable for your hardware core
global DDS_FRAC_BITS; 

% Storage for results
error_log = zeros(nSeeds, length(Q_list));      % MAE against Ideal
sqnr_log  = zeros(nSeeds, length(Q_list));      % SQNR against Single
sqnr_ideal_log = zeros(nSeeds, length(Q_list)); % SQNR against Ideal

% Storage for plotting (Seed #1 only)
seed1_data = struct(); 

% =========================================================
% 2. TEST LOOPS (SWAPPED: Q IS OUTER, SEED IS INNER)
% =========================================================
% --- OUTER LOOP: Iterate through Q-Formats first ---
for q_idx = 1:length(Q_list)
    
    % Set Global Q-Format for this entire pass
    if strcmp(current_mode, 'fixed')
        DDS_FRAC_BITS = Q_list(q_idx);
        fprintf('Testing Q2.%d ...\n', DDS_FRAC_BITS);
    else
        fprintf('Testing Standard Precision...\n');
    end
    seed_sqnrs = zeros(nSeeds, 1);
    
    % --- INNER LOOP: Iterate through Random Seeds ---
    for seed = 1:nSeeds
        
        % CRITICAL: Reset RNG to 'seed' so that Seed #1 is 
        % the EXACT SAME chirp for Q15, Q13, Q11, etc.
        rng(seed); 
        
        % -----------------------------------------------------
        % A. CHIRP PARAMETERS (Your Fixed Spec)
        % -----------------------------------------------------
        f0 = 0;
        B  = 200e6;
        
        % Instantaneous frequency law
        f_inst = f0 + B * (t / T_dur);
        
        % Nyquist Check
        if max(f_inst) >= Fs/2
            f_inst(f_inst >= Fs/2) = (Fs/2) - 1e3;
        end
        
        % Calculate Tuning Word (M)
        M = round(f_inst * (2^Nacc) / Fs);
        
        % -----------------------------------------------------
        % B. IDEAL REFERENCE (GROUND TRUTH)
        % -----------------------------------------------------
        ideal_phase = zeros(1, Ns);
        ideal_phase(2:end) = cumsum(double(M(1:end-1)));
        yExpected = sin(2*pi * ideal_phase / 2^Nacc);
        
        % -----------------------------------------------------
        % C. SINGLE PRECISION REFERENCE (For SQNR)
        % -----------------------------------------------------
        dds_ref_single = double(dds_core(M, Nacc, LUT_bits, 'single'));
        
        % Use MEAN for power calculations to be consistent
        P_signal_ref   = mean(dds_ref_single.^2); 
        P_signal_ideal = mean(yExpected.^2); 
        
        % -----------------------------------------------------
        % D. CORE EXECUTION & METRICS
        % -----------------------------------------------------
        
        % 1. Run Core (DDS_FRAC_BITS is already set by outer loop)
        dds_out = dds_core(M, Nacc, LUT_bits, current_mode);
        dds_out_dbl = double(dds_out);
        
        % 2. Compute MAE
        error_log(seed, q_idx) = mean(abs(dds_out_dbl - yExpected));
        
        % 3. Compute SQNR vs SINGLE
        noise_signal = dds_out_dbl - dds_ref_single;
        P_noise = mean(noise_signal.^2);
        if P_noise == 0, P_noise = eps; end
        sqnr_log(seed, q_idx) = 10 * log10(P_signal_ref / P_noise);
        
        % 4. Compute SQNR vs IDEAL
        noise_ideal = dds_out_dbl - yExpected;
        P_noise_ideal = mean(noise_ideal.^2);
        if P_noise_ideal == 0, P_noise_ideal = eps; end
        sqnr_ideal_log(seed, q_idx) = 10 * log10(P_signal_ideal / P_noise_ideal);
        
        % 5. Capture Data for Seed #1 Plotting
        if seed == 1
            seed1_data.time = t;
            seed1_data.ideal = yExpected;
            seed1_data.f_inst = f_inst;
            seed1_data.M = M;
            
            % Save specific waveform for this Q-format
            seed1_data.waves{q_idx} = dds_out;
            seed1_data.labels{q_idx} = ['Q2.' num2str(Q_list(q_idx))];
        end
        SQNR_fixed(seed) = calculate_SQNR(yExpected, dds_out_dbl);
    end
    
    %% Plot Fixed-Point SQNR Results
    figure('Position', [100, 100, 1200, 500]);
    
    subplot(1,2,1);
    plot(1:nSeeds, SQNR_fixed, 'm-o', ...
        'LineWidth', 1.5, 'MarkerSize', 5, 'MarkerFaceColor', 'm');
    grid on;
    xlabel('Seed Number', 'FontSize', 11, 'FontWeight', 'bold');
    ylabel('SQNR (dB)', 'FontSize', 11, 'FontWeight', 'bold');
    
    % Mean line
    hold on;
    yline(mean(SQNR_fixed), 'k--', 'LineWidth', 2, ...
        'Label', sprintf('Mean = %.2f dB', mean(SQNR_fixed)), ...
        'LabelHorizontalAlignment', 'left');
    hold off;
    
    fprintf('Mean SQNR (fixed): %.2f dB\n', mean(SQNR_fixed));
    fprintf('Min  SQNR (fixed): %.2f dB\n', min(SQNR_fixed));
    fprintf('Max  SQNR (fixed): %.2f dB\n', max(SQNR_fixed));
    fprintf('Std  SQNR (fixed): %.2f dB\n', std(SQNR_fixed));
    
    % Record Statistics
    avg_sqnr_results(q_idx) = mean(SQNR_fixed);
    min_sqnr_results(q_idx) = min(SQNR_fixed);
end

%% 3. Plotting the Trade-off Curve
figure('Position', [100, 100, 900, 600]);
% Main Curve
plot(Q_list, avg_sqnr_results, 'b-o', 'LineWidth', 2, 'MarkerFaceColor', 'b');
hold on;
plot(Q_list, min_sqnr_results, 'r--', 'LineWidth', 1.5);
% Annotations
grid on;
xlabel('Fractional Length (bits)');
ylabel('SQNR (dB)');
legend('Mean SQNR', 'Worst-Case SQNR', 'Location', 'Best');
fprintf('Simulation Complete.\n');

% =========================================================
% 3. RESULTS VISUALIZATION 
% =========================================================
% --- FIGURE 1: Error Statistics ---
figure('Name', 'Error Analysis', 'Color', 'w');
if strcmp(current_mode, 'fixed')
    avg_errors = mean(error_log, 1);
    bar(categorical(seed1_data.labels), avg_errors);
    xlabel('Fixed-Point Format');
    ylabel('Mean Absolute Error (Avg over 60 Seeds)');
    title(['Impact of Bit Depth on Accuracy (' current_mode ')']);
    grid on;
else
    plot(1:nSeeds, error_log(:,1), '.-');
    xlabel('Seed Index'); ylabel('MAE');
    title(['Error per Seed (' current_mode ')']);
end

% --- FIGURE 2: Waveform Comparison ---
figure('Name', 'Waveform Comparison', 'Color', 'w');
subplot(2,1,1);
plot(seed1_data.time*1e6, seed1_data.ideal, 'k', 'LineWidth', 1.5); hold on;
if strcmp(current_mode, 'fixed')
    plot(seed1_data.time*1e6, double(seed1_data.waves{1}), '--'); 
    if length(Q_list) > 1
        plot(seed1_data.time*1e6, double(seed1_data.waves{end}), ':'); 
        legend('Ideal', seed1_data.labels{1}, seed1_data.labels{end});
    else
        legend('Ideal', seed1_data.labels{1});
    end
else
    plot(seed1_data.time*1e6, double(seed1_data.waves{1}), '--');
    legend('Ideal', current_mode);
end
title('Waveform: Ideal vs DDS Output (Seed #1)');
xlabel('Time (\mus)'); ylabel('Amplitude');
grid on; xlim([0 2]); 

subplot(2,1,2);
f_dds_recalc = double(seed1_data.M) * Fs / 2^Nacc;
plot(seed1_data.time*1e6, seed1_data.f_inst/1e6, 'k', 'LineWidth', 2); hold on;
plot(seed1_data.time*1e6, f_dds_recalc/1e6, '--r');
legend('Target Chirp', 'DDS M-Word');
title('Instantaneous Frequency Tracking');
xlabel('Time (\mus)'); ylabel('Frequency (MHz)');
grid on;

% =========================================================
% 4. TEXT SUMMARY
% =========================================================
fprintf('\n===========================================\n');
fprintf('           RESULTS SUMMARY                 \n');
fprintf('===========================================\n');
fprintf('Testing Mode: %s\n', current_mode);
if strcmp(current_mode, 'fixed')
    for q_idx = 1:length(Q_list)
        current_errors = error_log(:, q_idx);
        fprintf('\n--- Q2.%d Format ---\n', Q_list(q_idx));
        fprintf('Average Error: %e\n', mean(current_errors));
        fprintf('Max Error    : %e\n', max(current_errors));
        fprintf('Min Error    : %e\n', min(current_errors));
    end
else
    current_errors = error_log(:, 1);
    fprintf('\n--- Standard Precision Analysis ---\n');
    fprintf('Average Error: %e\n', mean(current_errors));
    fprintf('Max Error    : %e\n', max(current_errors));
    fprintf('Min Error    : %e\n', min(current_errors));
end
fprintf('\n===========================================\n');

% =========================================================
% 5. FIGURE: ERROR VS SEED
% =========================================================
figure('Name', 'DDS Accuracy vs Seed', 'Color', 'w');
plot(1:nSeeds, error_log, '.-', 'LineWidth', 1.5);
grid on;
xlabel('Test Seed'); ylabel('Mean Absolute Error');
title(['DDS Chirp Accuracy - ' current_mode ' Mode']);
subtitle('Error due to Phase Truncation and LUT Quantization Only');
if strcmp(current_mode, 'fixed')
    legends = cell(1, length(Q_list));
    for k = 1:length(Q_list)
        legends{k} = ['Q2.' num2str(Q_list(k))];
    end
    legend(legends, 'Location', 'best');
else
    legend('Double/Single Precision');
end
xlim([1 nSeeds]);

% =========================================================
% 6. FIGURE: SQNR VS BITS (Ref: Single Precision)
% =========================================================
if strcmp(current_mode, 'fixed')
    figure('Name', 'SQNR Analysis (HW Noise Only)', 'Color', 'w');
    avg_sqnr = mean(sqnr_log, 1);
    plot(Q_list, avg_sqnr, '-o', 'LineWidth', 2, 'MarkerFaceColor', 'b'); hold on;
    theory_line = 6.02 * Q_list;
    offset = avg_sqnr(1) - theory_line(1); 
    plot(Q_list, theory_line + offset, '--r', 'LineWidth', 1.5);
    grid on;
    xlabel('Fractional Bits (Q2.x)'); ylabel('SQNR (dB)');
    title('DDS Output SQNR vs. Bit Depth');
    subtitle('Reference: Single Precision (Shows Quantization Noise Only)');
    legend('Measured SQNR', 'Theoretical 6dB/bit Trend', 'Location', 'best');
    xticks(sort(Q_list)); 
end

% =========================================================
% 7. FIGURE: SQNR VS BITS (Ref: Ideal Mathematical Model)
% =========================================================
if strcmp(current_mode, 'fixed')
    figure('Name', 'Total SQNR vs Ideal', 'Color', 'w');
    avg_sqnr_ideal = mean(sqnr_ideal_log, 1);
    plot(Q_list, avg_sqnr_ideal, '-s', 'LineWidth', 2, 'MarkerFaceColor', 'm', 'Color', 'm'); hold on;
    theory_line = 6.02 * Q_list;
    offset = avg_sqnr_ideal(end) - theory_line(end); 
    plot(Q_list, theory_line + offset, '--k', 'LineWidth', 1.5);
    grid on;
    xlabel('Fractional Bits (Q2.x)'); ylabel('SQNR (dB)');
    title('Total SQNR vs. Bit Depth');
    subtitle('Reference: Ideal Math Model (Shows Phase Truncation Ceiling)');
    legend('Total SQNR (Quantization + Phase Noise)', 'Theoretical 6dB/bit Trend', 'Location', 'best');
    xticks(sort(Q_list)); 
end

% =========================================================
% 8. FIGURE: POWER SPECTRAL DENSITY (PSD)
% =========================================================
figure('Name', 'Power Spectral Density', 'Color', 'w', 'Position', [150, 150, 900, 500]);
hold on;

% Plot Ideal PSD using a rectangular window (best for chirps)
[pxx_ideal, f_psd] = periodogram(seed1_data.ideal, rectwin(length(seed1_data.ideal)), length(seed1_data.ideal), Fs);
plot(f_psd/1e6, 10*log10(pxx_ideal), 'k', 'LineWidth', 2, 'DisplayName', 'Ideal Chirp');

% Plot Quantized DDS PSDs
colors = lines(length(Q_list));
if strcmp(current_mode, 'fixed')
    for q_idx = 1:length(Q_list)
        [pxx_dds, ~] = periodogram(double(seed1_data.waves{q_idx}), rectwin(length(seed1_data.waves{q_idx})), length(seed1_data.waves{q_idx}), Fs);
        plot(f_psd/1e6, 10*log10(pxx_dds), 'Color', colors(q_idx,:), ...
            'LineWidth', 1, 'DisplayName', seed1_data.labels{q_idx});
    end
else
    [pxx_dds, ~] = periodogram(double(seed1_data.waves{1}), rectwin(length(seed1_data.waves{1})), length(seed1_data.waves{1}), Fs);
    plot(f_psd/1e6, 10*log10(pxx_dds), 'b--', 'LineWidth', 1, 'DisplayName', current_mode);
end

grid on;
xlabel('Frequency (MHz)', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Power/Frequency (dB/Hz)', 'FontSize', 11, 'FontWeight', 'bold');
title('Power Spectral Density (PSD) of DDS Chirp Output', 'FontSize', 12, 'FontWeight', 'bold');
subtitle('Evaluated using Seed #1');
legend('Location', 'best');
xlim([0 Fs/2/1e6]); % Plot up to Nyquist frequency (250 MHz)
ylim([-120 max(10*log10(pxx_ideal))+10]); % Set a clean noise floor view
hold off;

%% SQNR Calculation Function
function sqnr_db = calculate_SQNR(signal_ref, signal_test)
    % Calculate Signal-to-Quantization Noise Ratio
    signal_power = sum(abs(signal_ref).^2);
    noise_power = sum(abs(signal_ref - signal_test).^2);
    if noise_power == 0
        sqnr_db = Inf;
    else
        sqnr_db = 10 * log10(signal_power / noise_power);
    end
end