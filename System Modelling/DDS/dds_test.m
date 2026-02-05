clear; clc; close all;

% =========================================================
% 1. DESIGN PARAMETERS
% =========================================================
Fs       = 500e6;        % System clock (Hz)
T_dur    = 20e-6;        % Signal duration (s) (Shortened slightly for speed)
Nacc     = 32;           % Phase accumulator bits
LUT_bits = 16;           % LUT address bits
nSeeds   = 60;           % Number of random tests

Ns = round(T_dur * Fs);  % Number of samples
t  = (0:Ns-1) / Fs;      % Time vector

% -----------------------------------------------------
% CONFIGURATION
% -----------------------------------------------------
current_mode = 'fixed'; % Options: 'fixed', 'double', 'single'
fprintf('Testing Mode: %s\n', current_mode);

% Define Q-Formats to sweep (Only used if mode is 'fixed')
if strcmp(current_mode, 'fixed')
    Q_list = [15 13 11 9]; 
else
    Q_list = [0]; % Dummy value for non-fixed modes
end

% Global variable for your hardware core
global DDS_FRAC_BITS; 

% Storage for results
% Rows = Seeds, Cols = Q_formats
error_log = zeros(nSeeds, length(Q_list)); 

% Storage for plotting (We only save Seed #1 waveforms to save RAM)
seed1_data = struct(); 

% =========================================================
% 2. RANDOMIZED DDS TEST LOOP
% =========================================================
fprintf('Running %d simulations...\n', nSeeds);

for seed = 1:nSeeds
    rng(seed);
    
    % -----------------------------------------------------
    % A. RANDOMIZED CHIRP PARAMETERS
    % -----------------------------------------------------
    f0_min = 1e6;   f0_max = 20e6;
    B_min  = 5e6;   B_max  = 80e6;
    
    f0 = f0_min + (f0_max - f0_min) * rand;   
    B  = B_min  + (B_max  - B_min)  * rand;   
    
    % Instantaneous frequency law (linear chirp)
    f_inst = f0 + B * (t / T_dur);
    
    % Nyquist Check
    if max(f_inst) >= Fs/2
        warning('Seed %d exceeded Nyquist. Clamping.', seed);
        f_inst(f_inst >= Fs/2) = (Fs/2) - 1e3;
    end
    
    % Calculate Tuning Word (M) for the whole vector
    M = round(f_inst * (2^Nacc) / Fs);

    % -----------------------------------------------------
    % B. IDEAL REFERENCE (GROUND TRUTH)
    % -----------------------------------------------------
    ideal_phase = zeros(1, Ns);
    ideal_phase(2:end) = cumsum(double(M(1:end-1)));
    yExpected = sin(2*pi * ideal_phase / 2^Nacc);
    
    % -----------------------------------------------------
    % C. DDS CORE EXECUTION & ERROR CALCULATION
    % -----------------------------------------------------
    
    % Loop through Q formats (if fixed) or run once (if double/single)
    for q_idx = 1:length(Q_list)
        
        % Set the global Q-format if in fixed mode
        if strcmp(current_mode, 'fixed')
            DDS_FRAC_BITS = Q_list(q_idx);
        end
        
        % Call your core
        % Note: Ensure your core accepts 'current_mode' correctly
        dds_out = dds_core(M, Nacc, LUT_bits, current_mode);
        
        % Compute Mean Absolute Error (MAE) for this specific Q-format
        current_error = mean(abs(double(dds_out) - yExpected));
        error_log(seed, q_idx) = current_error;
        
        % -------------------------------------------------
        % D. CAPTURE DATA FOR PLOTTING (SEED 1 ONLY)
        % -------------------------------------------------
        if seed == 1
            % Save ideal reference once
            if q_idx == 1
                seed1_data.time = t;
                seed1_data.ideal = yExpected;
                seed1_data.f_inst = f_inst;
                seed1_data.M = M;
            end
            % Save the waveform for this Q format
            seed1_data.waves{q_idx} = dds_out;
            seed1_data.labels{q_idx} = ['Q1.' num2str(Q_list(q_idx))];
        end
    end
end

fprintf('Simulation Complete.\n');

% =========================================================
% 3. RESULTS VISUALIZATION 
% =========================================================

% --- FIGURE 1: Error Statistics Comparison ---
figure('Name', 'Error Analysis', 'Color', 'w');
if strcmp(current_mode, 'fixed')
    % Calculate average error across all seeds for each Q format
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

% --- FIGURE 2: Waveform Comparison (Seed #1) ---
figure('Name', 'Waveform Comparison', 'Color', 'w');
subplot(2,1,1);
plot(seed1_data.time*1e6, seed1_data.ideal, 'k', 'LineWidth', 1.5); hold on;
% Plot the lowest and highest precision for comparison
if strcmp(current_mode, 'fixed')
    plot(seed1_data.time*1e6, double(seed1_data.waves{1}), '--'); 
    plot(seed1_data.time*1e6, double(seed1_data.waves{end}), ':'); 
    legend('Ideal', seed1_data.labels{1}, seed1_data.labels{end});
else
    plot(seed1_data.time*1e6, double(seed1_data.waves{1}), '--');
    legend('Ideal', current_mode);
end
title('Waveform: Ideal vs DDS Output (Seed #1)');
xlabel('Time (\mus)'); ylabel('Amplitude');
grid on;
xlim([0 2]); % Zoom in on the start to see details

% --- FIGURE 3: Chirp Frequency Verification ---
subplot(2,1,2);
f_dds_recalc = double(seed1_data.M) * Fs / 2^Nacc;
plot(seed1_data.time*1e6, seed1_data.f_inst/1e6, 'k', 'LineWidth', 2); hold on;
plot(seed1_data.time*1e6, f_dds_recalc/1e6, '--r');
legend('Target Chirp', 'DDS M-Word');
title('Instantaneous Frequency Tracking');
xlabel('Time (\mus)'); ylabel('Frequency (MHz)');
grid on;

% =========================================================
% 4. TEXT SUMMARY OUTPUT
% =========================================================
fprintf('\n===========================================\n');
fprintf('           RESULTS SUMMARY                 \n');
fprintf('===========================================\n');
fprintf('Testing Mode: %s\n', current_mode);

if strcmp(current_mode, 'fixed')
    % LOOP through each Q format and print stats
    for q_idx = 1:length(Q_list)
        % Extract the error column for this specific Q format across all seeds
        current_errors = error_log(:, q_idx);
        
        fprintf('\n--- Q1.%d Format ---\n', Q_list(q_idx));
        fprintf('Average Error: %e\n', mean(current_errors));
        fprintf('Max Error    : %e\n', max(current_errors));
        fprintf('Min Error    : %e\n', min(current_errors));
    end
else
    % SINGLE/DOUBLE mode (only one column of errors)
    current_errors = error_log(:, 1);
    
    fprintf('\n--- Standard Precision Analysis ---\n');
    fprintf('Average Error: %e\n', mean(current_errors));
    fprintf('Max Error    : %e\n', max(current_errors));
    fprintf('Min Error    : %e\n', min(current_errors));
end
fprintf('\n===========================================\n');

% =========================================================
% 5. FIGURE: ERROR VS SEED (COMPARISON PLOT)
% =========================================================
figure('Name', 'DDS Accuracy vs Seed', 'Color', 'w');

% Plot all columns of error_log (each column is a different Q format)
% Columns correspond to the Q_list order
p = plot(1:nSeeds, error_log, '.-', 'LineWidth', 1.5);

grid on;
xlabel('Test Seed');
ylabel('Mean Absolute Error');
title(['DDS Chirp Accuracy - ' current_mode ' Mode']);
subtitle('Error due to Phase Truncation and LUT Quantization Only');

% Dynamic Legend Generation
if strcmp(current_mode, 'fixed')
    % Create legend labels based on your Q_list
    legends = cell(1, length(Q_list));
    for k = 1:length(Q_list)
        legends{k} = ['Q1.' num2str(Q_list(k))];
    end
    legend(legends, 'Location', 'best');
    
    % Optional: Make the lines distinct colors/styles if needed
    % (MATLAB does this automatically for multiple columns)
else
    legend('Double/Single Precision');
end

% Set axis limits to look nice (optional)
xlim([1 nSeeds]);