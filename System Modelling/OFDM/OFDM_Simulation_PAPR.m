%% OFDM: Three Channel Conditions - Rayleigh only, AWGN only, Both
clear; clc; close all;
%% ---------------- Parameters ----------------
NFFT = 4096;
SCS = 120e3;                 % Subcarrier spacing [Hz]
NumSCinRB = 12;
NusedFFT = 264 * NumSCinRB;  % active subcarriers
NGuardSC = NFFT - NusedFFT;
NSymbols_data = 100;         % number of data symbols
numPreamble = 3;             % number of identical preamble OFDM symbols (will average)
totalSymbols = numPreamble + NSymbols_data; % OFDM modulator NumSymbols
ModOrder = 256;              % 256-QAM
% Channel / simulation
pathDelays = [0 0.3e-8 0.8e-8 2.5e-8];  % seconds
pathPowerdB = [0 -3 -6 -9];
maxDopplerHz = 20;            % set >0 for mobility tests

% SNR sweep for AWGN cases
SNR_dB_vec = 0:2:60;

% Pilot / estimation settings
pilotSpacingSC = 8;          % pilots every 8 subcarriers
pilotSpacingSym = 2;         % pilots every 2 OFDM symbols (in data block)
pilotBoostdB = 6;            % pilot amplitude boost
pilotVal = 1+0j;

%% ---------------- Derived & auto-CP ----------------
fs = NFFT * SCS;
Ts = 1/fs;
maxDelay = max(pathDelays);
cpLen = ceil(maxDelay / Ts) + 4;  % small safety margin

fprintf('cpLen = %d samples (maxDelay=%g s)\n', cpLen, maxDelay);

%% ---------------- OFDM objects ----------------
ofdmmod = comm.OFDMModulator('FFTLength', NFFT, ...
    'NumSymbols', totalSymbols, ...
    'NumGuardBandCarriers', [NGuardSC/2; NGuardSC/2], ...
    'CyclicPrefixLength', cpLen, ...
    'InsertDCNull', false);
ofdmdemod = comm.OFDMDemodulator('FFTLength', NFFT, ...
    'NumSymbols', totalSymbols, ...
    'NumGuardBandCarriers', [NGuardSC/2; NGuardSC/2], ...
    'CyclicPrefixLength', cpLen, ...
    'RemoveDCCarrier', false);
ofdmInfo = info(ofdmmod);
numSC = ofdmInfo.DataInputSize(1);
numSym = ofdmInfo.DataInputSize(2);
fprintf('OFDM grid: %d subcarriers x %d symbols (preamble=%d, data=%d)\n', numSC, numSym, numPreamble, NSymbols_data);

%% ---------------- Common Data Generation ----------------
% Generate random data (same for all simulations)
rng(42); % Set seed for reproducibility
msg = randi([0 ModOrder-1], numSC, NSymbols_data);
dataQAM = qammod(msg, ModOrder, 'UnitAveragePower', true);

% Pilot mask definition over data block (size: numSC x NSymbols_data)
pilotMask = false(numSC, NSymbols_data);
pilotMask(1:pilotSpacingSC:end, 1:pilotSpacingSym:end) = true;

% Build txGrid: columns 1..numPreamble are preamble, columns numPreamble+1..end are data/pilots
txGrid = zeros(numSC, numSym);
preambleRef = ones(numSC, numPreamble);  % known preamble pattern (unit amplitude)
txGrid(:, 1:numPreamble) = preambleRef;
txGrid(:, numPreamble+1:end) = dataQAM;  % initial data

% Insert boosted pilots into txGrid at correct positions using sub2ind
[pRows_rel, pCols_rel] = find(pilotMask);                % indices relative to data block
pCols_full = pCols_rel + numPreamble;                    % convert to full-frame columns
pilotBoost = 10^(pilotBoostdB/20);
linearIdxPilots = sub2ind(size(txGrid), pRows_rel, pCols_full);
txGrid(linearIdxPilots) = pilotBoost * pilotVal;         % place boosted pilots (overwrites data at those positions)

txDataGrid = txGrid(:, numPreamble+1:end); % Transmitted data block (with pilots)
dataMask = ~pilotMask;  % relative to data block
refGrid_data_only = txDataGrid(dataMask); % Transmitted data symbols (for EVM)
msg_data_only = msg(dataMask); % Transmitted data integers (for BER)
k = log2(ModOrder); % Bits per symbol
txBits = int2bit(msg_data_only, k);

%% ---------------- OFDM modulate ----------------
txSig = ofdmmod(txGrid);

% Calculate samples per symbol for plotting
samples_per_symbol = NFFT + cpLen;

% Channel Object for Rayleigh cases
rayleighChan = comm.RayleighChannel('SampleRate', fs, ...
    'PathDelays', pathDelays, 'AveragePathGains', pathPowerdB, ...
    'NormalizePathGains', true, 'MaximumDopplerShift', maxDopplerHz);

%% ---------------- Three Simulation Cases ----------------
channel_types = {'Rayleigh Only', 'AWGN Only', 'Rayleigh + AWGN'};
colors = {'r', 'g', 'b'};
markers = {'o', 's', 'd'};

% Pre-allocate results storage
results = struct();

% For Rayleigh Only case, we'll do SNR sweep too
for case_idx = 1:3
    channel_type = channel_types{case_idx};
    fprintf('\n=== Running Simulation: %s ===\n', channel_type);
    
    switch case_idx
        case 1 % Rayleigh Only - NOW WITH SNR SWEEP
            % --- SNR sweep with Rayleigh only ---
            reset(rayleighChan);
            rxSig_faded = rayleighChan(txSig);
            
            [results(case_idx).rmsEVM_vec_pilot, results(case_idx).rmsEVM_vec_preamble, ...
             results(case_idx).rmsEVM_vec_noEQ, results(case_idx).ber_vec_pilot, ...
             results(case_idx).rxSig, results(case_idx).rxGrid, ...
             results(case_idx).H_est_pilot_grid, results(case_idx).eqPilotVec, ...
             results(case_idx).eqPreambleVec, results(case_idx).rxVec_noEQ, ...
             results(case_idx).last_EVM] = ...
                processOFDM_SNR(rxSig_faded, txGrid, txDataGrid, refGrid_data_only, ...
                              txBits, ModOrder, dataMask, pilotMask, numPreamble, ...
                              preambleRef, pilotBoost, pilotVal, ofdmdemod, SNR_dB_vec, fs, 'rayleigh');
            
            results(case_idx).SNR_dB_vec = SNR_dB_vec;
            
        case 2 % AWGN Only
            % --- SNR sweep with AWGN only ---
            [results(case_idx).rmsEVM_vec_pilot, results(case_idx).rmsEVM_vec_preamble, ...
             results(case_idx).rmsEVM_vec_noEQ, results(case_idx).ber_vec_pilot, ...
             results(case_idx).rxSig, results(case_idx).rxGrid, ...
             results(case_idx).H_est_pilot_grid, results(case_idx).eqPilotVec, ...
             results(case_idx).eqPreambleVec, results(case_idx).rxVec_noEQ, ...
             results(case_idx).last_EVM] = ...
                processOFDM_SNR(txSig, txGrid, txDataGrid, refGrid_data_only, ...
                              txBits, ModOrder, dataMask, pilotMask, numPreamble, ...
                              preambleRef, pilotBoost, pilotVal, ofdmdemod, SNR_dB_vec, fs, 'awgn');
            
            results(case_idx).SNR_dB_vec = SNR_dB_vec;
            
        case 3 % Rayleigh + AWGN
            % --- SNR sweep with Rayleigh + AWGN ---
            reset(rayleighChan);
            rxSig_faded = rayleighChan(txSig);
            
            [results(case_idx).rmsEVM_vec_pilot, results(case_idx).rmsEVM_vec_preamble, ...
             results(case_idx).rmsEVM_vec_noEQ, results(case_idx).ber_vec_pilot, ...
             results(case_idx).rxSig, results(case_idx).rxGrid, ...
             results(case_idx).H_est_pilot_grid, results(case_idx).eqPilotVec, ...
             results(case_idx).eqPreambleVec, results(case_idx).rxVec_noEQ, ...
             results(case_idx).last_EVM] = ...
                processOFDM_SNR(rxSig_faded, txGrid, txDataGrid, refGrid_data_only, ...
                              txBits, ModOrder, dataMask, pilotMask, numPreamble, ...
                              preambleRef, pilotBoost, pilotVal, ofdmdemod, SNR_dB_vec, fs, 'rayleigh_awgn');
            
            results(case_idx).SNR_dB_vec = SNR_dB_vec;
    end
    
    % Calculate PAPR for transmitted and received signals
    results(case_idx).PAPR_tx = calculatePAPR(txSig);
    results(case_idx).PAPR_rx = calculatePAPR(results(case_idx).rxSig);
    
    % Store additional parameters for plotting
    results(case_idx).fs = fs;
    results(case_idx).numPreamble = numPreamble;
    results(case_idx).SCS = SCS;
    results(case_idx).numSC = numSC;
    results(case_idx).samples_per_symbol = samples_per_symbol;
    
    % Generate figures for this case - pass txSig and fs
    generateFigures(results(case_idx), case_idx, channel_type, colors{case_idx}, markers{case_idx}, txSig, fs);
end

%% ---------------- Combined Comparison Figures ----------------
fprintf('\n=== Generating Combined Comparison Figures ===\n');

% Figure 19: Combined EVM vs SNR
figure(19);
hold on;
for case_idx = 1:3
    % All cases now have SNR sweep - plot vs SNR
    plot(results(case_idx).SNR_dB_vec, results(case_idx).rmsEVM_vec_pilot, ...
         [colors{case_idx} markers{case_idx} '-'], 'LineWidth', 1.5, ...
         'DisplayName', [channel_types{case_idx} ' (Pilot MMSE)']);
    plot(results(case_idx).SNR_dB_vec, results(case_idx).rmsEVM_vec_preamble, ...
         [colors{case_idx} markers{case_idx} '--'], 'LineWidth', 1.5, ...
         'DisplayName', [channel_types{case_idx} ' (Preamble MMSE)']);
end
hold off;
title('Combined: EVM vs. SNR - All Channel Types');
xlabel('SNR (dB)');
ylabel('RMS EVM (%)');
legend('show', 'Location', 'northeast');
grid on;
set(gca, 'YScale', 'log');

% Figure 20: Combined BER vs SNR
figure(20);
hold on;
for case_idx = 1:3
    % All cases now have SNR sweep - plot vs SNR
    semilogy(results(case_idx).SNR_dB_vec, results(case_idx).ber_vec_pilot, ...
             [colors{case_idx} markers{case_idx} '-'], 'LineWidth', 1.5, ...
             'DisplayName', [channel_types{case_idx} ' (Pilot MMSE)']);
end
hold off;
title('Combined: BER vs. SNR - All Channel Types');
xlabel('SNR (dB)');
ylabel('Bit Error Rate (BER)');
legend('show', 'Location', 'southwest');
grid on;

% Figure 21: Combined PAPR Comparison
figure(21);
PAPR_tx = [results.PAPR_tx];
PAPR_rx = [results.PAPR_rx];
bar_data = [PAPR_tx; PAPR_rx]';
bar_handle = bar(bar_data);
set(gca, 'XTickLabel', channel_types);
ylabel('PAPR (dB)');
title('PAPR Comparison: Transmitted vs Received Signals');
legend('Transmitted', 'Received', 'Location', 'northeast');
grid on;

% Add PAPR values on top of bars
for i = 1:length(channel_types)
    text(i-0.2, PAPR_tx(i)+0.1, sprintf('%.2f dB', PAPR_tx(i)), 'HorizontalAlignment', 'center');
    text(i+0.2, PAPR_rx(i)+0.1, sprintf('%.2f dB', PAPR_rx(i)), 'HorizontalAlignment', 'center');
end

fprintf('\n=== All simulations complete! Generated 21 figures total ===\n');
fprintf('Figures 1-6: Rayleigh Only\n');
fprintf('Figures 7-12: AWGN Only\n');
fprintf('Figures 13-18: Rayleigh + AWGN\n');
fprintf('Figures 19-20: Combined Comparison\n');
fprintf('Figure 21: PAPR Comparison\n');

%% ---------------- PAPR Calculation Function ----------------
function papr = calculatePAPR(signal)
    % Calculate Peak-to-Average Power Ratio (PAPR) in dB
    power = abs(signal).^2;
    peak_power = max(power);
    avg_power = mean(power);
    papr = 10*log10(peak_power / avg_power);
end

%% ---------------- Support Functions ----------------
function [rmsEVM_pilot, rmsEVM_preamble, rmsEVM_noEQ, ber_pilot, ...
          rxSig, rxGrid, H_est_pilot_grid, eqPilotVec, eqPreambleVec, rxVec_noEQ] = ...
    processOFDM(rxSig, txGrid, txDataGrid, refGrid_data_only, txBits, ModOrder, ...
               dataMask, pilotMask, numPreamble, preambleRef, pilotBoost, pilotVal, ofdmdemod, noiseVar, fs)
    
    % Demodulate
    rxGrid = ofdmdemod(rxSig);
    
    % Separate preamble RX (average the numPreamble columns)
    rxPreamble = rxGrid(:, 1:numPreamble);
    rxPreambleAvg = mean(rxPreamble, 2);
    
    % Received data area
    rxDataGrid = rxGrid(:, numPreamble+1:end);
    
    % Channel estimation and equalization
    [H_est_pilot_grid, eqPilotVec, eqPreambleVec, rxVec_noEQ] = ...
        channelEstimationAndEqualization(rxGrid, rxDataGrid, txGrid, dataMask, ...
                                       pilotMask, numPreamble, preambleRef, ...
                                       rxPreambleAvg, pilotBoost, pilotVal, noiseVar, fs);
    
    % EVM Measurement
    evmObj = comm.EVM('AveragingDimensions', [1 2]);
    rmsEVM_noEQ = evmObj(refGrid_data_only, rxVec_noEQ);
    rmsEVM_pilot = evmObj(refGrid_data_only, eqPilotVec);
    rmsEVM_preamble = evmObj(refGrid_data_only, eqPreambleVec);
    
    % BER Calculation
    rxMsg = qamdemod(eqPilotVec, ModOrder, 'UnitAveragePower', true);
    k = log2(ModOrder);
    rxBits = int2bit(rxMsg, k);
    berObj = comm.ErrorRate;
    ber_results = berObj(txBits, rxBits);
    ber_pilot = ber_results(1);
end

function [rmsEVM_vec_pilot, rmsEVM_vec_preamble, rmsEVM_vec_noEQ, ber_vec_pilot, ...
          rxSig_last, rxGrid_last, H_est_pilot_grid_last, eqPilotVec_last, ...
          eqPreambleVec_last, rxVec_noEQ_last, last_EVM] = ...
    processOFDM_SNR(txSig_faded, txGrid, txDataGrid, refGrid_data_only, txBits, ModOrder, ...
                   dataMask, pilotMask, numPreamble, preambleRef, pilotBoost, pilotVal, ofdmdemod, SNR_dB_vec, fs, channel_type)
    
    % Pre-allocate arrays
    rmsEVM_vec_pilot = zeros(size(SNR_dB_vec));
    rmsEVM_vec_preamble = zeros(size(SNR_dB_vec));
    rmsEVM_vec_noEQ = zeros(size(SNR_dB_vec));
    ber_vec_pilot = zeros(size(SNR_dB_vec));
    
    evmObj = comm.EVM('AveragingDimensions', [1 2]);
    berObj = comm.ErrorRate;
    
    for snr_idx = 1:length(SNR_dB_vec)
        current_SNR_dB = SNR_dB_vec(snr_idx);
        
        % Add AWGN based on channel type
        if strcmp(channel_type, 'rayleigh')
            % For Rayleigh Only - use very high SNR to simulate no AWGN
            rxSig = awgn(txSig_faded, 100, 'measured');
        elseif strcmp(channel_type, 'awgn')
            % For AWGN Only - use the specified SNR
            rxSig = awgn(txSig_faded, current_SNR_dB, 'measured');
        else % rayleigh_awgn
            % For Rayleigh + AWGN - use the specified SNR
            rxSig = awgn(txSig_faded, current_SNR_dB, 'measured');
        end
        
        % Demodulate
        rxGrid = ofdmdemod(rxSig);
        
        % Separate preamble RX
        rxPreamble = rxGrid(:, 1:numPreamble);
        rxPreambleAvg = mean(rxPreamble, 2);
        
        % Received data area
        rxDataGrid = rxGrid(:, numPreamble+1:end);
        
        % Estimate noise variance
        signalPower = mean(abs(txDataGrid(dataMask)).^2, 'all');
        noiseVar = signalPower / (10^(current_SNR_dB/10));
        
        % Channel estimation and equalization
        [H_est_pilot_grid, eqPilotVec, eqPreambleVec, rxVec_noEQ] = ...
            channelEstimationAndEqualization(rxGrid, rxDataGrid, txGrid, dataMask, ...
                                           pilotMask, numPreamble, preambleRef, ...
                                           rxPreambleAvg, pilotBoost, pilotVal, noiseVar, fs);
        
        % EVM Measurement
        rmsEVM_vec_noEQ(snr_idx) = evmObj(refGrid_data_only, rxVec_noEQ);
        rmsEVM_vec_pilot(snr_idx) = evmObj(refGrid_data_only, eqPilotVec);
        rmsEVM_vec_preamble(snr_idx) = evmObj(refGrid_data_only, eqPreambleVec);
        
        % BER Calculation
        rxMsg = qamdemod(eqPilotVec, ModOrder, 'UnitAveragePower', true);
        k = log2(ModOrder);
        rxBits = int2bit(rxMsg, k);
        ber_results = berObj(txBits, rxBits);
        ber_vec_pilot(snr_idx) = ber_results(1);
        reset(berObj);
        
        % Store last iteration results (for highest SNR)
        if snr_idx == length(SNR_dB_vec)
            rxSig_last = rxSig;
            rxGrid_last = rxGrid;
            H_est_pilot_grid_last = H_est_pilot_grid;
            eqPilotVec_last = eqPilotVec;
            eqPreambleVec_last = eqPreambleVec;
            rxVec_noEQ_last = rxVec_noEQ;
            
            % Store EVM values for constellation plots in a separate structure
            last_EVM.rmsEVM_pilot = rmsEVM_vec_pilot(snr_idx);
            last_EVM.rmsEVM_preamble = rmsEVM_vec_preamble(snr_idx);
            last_EVM.rmsEVM_noEQ = rmsEVM_vec_noEQ(snr_idx);
        end
    end
end

function [H_est_pilot_grid, eqPilotVec, eqPreambleVec, rxVec_noEQ] = ...
    channelEstimationAndEqualization(rxGrid, rxDataGrid, txGrid, dataMask, ...
                                   pilotMask, numPreamble, preambleRef, ...
                                   rxPreambleAvg, pilotBoost, pilotVal, noiseVar, fs)
    
    numSC = size(rxGrid, 1);
    numSym = size(rxGrid, 2);
    NSymbols_data = size(rxDataGrid, 2);
    
    % Pilot-based LS & interpolation
    [pRows_rel, pCols_rel] = find(pilotMask);
    pCols_full = pCols_rel + numPreamble;
    Yp = rxGrid(sub2ind(size(rxGrid), pRows_rel, pCols_full));
    Xp = pilotBoost * pilotVal * ones(size(Yp));
    Hp = Yp ./ (Xp + eps);
    
    xPil = double(pCols_full(:));
    yPil = double(pRows_rel(:));
    pvals_re = real(Hp(:));
    pvals_im = imag(Hp(:));
    
    [symGrid, scGrid] = meshgrid(1:numSym, 1:numSC);
    F_re = scatteredInterpolant(xPil, yPil, pvals_re, 'linear', 'nearest');
    F_im = scatteredInterpolant(xPil, yPil, pvals_im, 'linear', 'nearest');
    
    H_est_pilot_grid = F_re(symGrid, scGrid) + 1j * F_im(symGrid, scGrid);
    H_est_pilot_grid = fillmissing(H_est_pilot_grid, 'nearest');
    
    % Smoothing
    k_smooth = [1 4 6 4 1]; k_smooth = k_smooth/sum(k_smooth);
    H_est_pilot_grid = conv2(k_smooth, k_smooth, H_est_pilot_grid, 'same');
    
    % Signal power for MMSE
    txDataGrid = txGrid(:, numPreamble+1:end);
    signalPower = mean(abs(txDataGrid(dataMask)).^2, 'all');
    
    % MMSE equalizer (pilot-based)
    H_use = H_est_pilot_grid;
    mmseWeight = conj(H_use) ./ (abs(H_use).^2 + noiseVar / (signalPower + eps));
    eqGrid_pilotMMSE = mmseWeight .* rxGrid;
    eqData_pilotMMSE = eqGrid_pilotMMSE(:, numPreamble+1:end);
    
    % CPE correction (pilot-based)
    eqData_pilotMMSE_CPE = eqData_pilotMMSE;
    for symRel = 1:NSymbols_data
        if any(pilotMask(:, symRel))
            rows = find(pilotMask(:, symRel));
            colFull = symRel + numPreamble;
            receivedPilotsEq = eqGrid_pilotMMSE(rows, colFull);
            ph = angle(mean(receivedPilotsEq ./ (pilotBoost * pilotVal + eps)));
            eqData_pilotMMSE_CPE(:, symRel) = eqData_pilotMMSE_CPE(:, symRel) * exp(-1j * ph);
        end
    end
    
    % Preamble-based MMSE
    H_preamble = rxPreambleAvg ./ (mean(preambleRef,2) + eps);
    H_preamble_grid = repmat(H_preamble, 1, numSym);
    mmseWeight_preamble = conj(H_preamble_grid) ./ (abs(H_preamble_grid).^2 + noiseVar / (signalPower + eps));
    eqGrid_preambleMMSE = mmseWeight_preamble .* rxGrid;
    eqData_preambleMMSE = eqGrid_preambleMMSE(:, numPreamble+1:end);
    
    % CPE correction (preamble-based)
    eqData_preambleMMSE_CPE = eqData_preambleMMSE;
    for symRel = 1:NSymbols_data
        if any(pilotMask(:, symRel))
            rows = find(pilotMask(:, symRel));
            colFull = symRel + numPreamble;
            receivedPilotsEq = eqGrid_preambleMMSE(rows, colFull);
            ph = angle(mean(receivedPilotsEq ./ (pilotBoost * pilotVal + eps)));
            eqData_preambleMMSE_CPE(:, symRel) = eqData_preambleMMSE_CPE(:, symRel) * exp(-1j * ph);
        end
    end
    
    % Extract data vectors
    rxVec_noEQ = rxDataGrid(dataMask);
    eqPilotVec = eqData_pilotMMSE_CPE(dataMask);
    eqPreambleVec = eqData_preambleMMSE_CPE(dataMask);
end

function generateFigures(results, case_idx, channel_type, color, marker, txSig, fs)
    base_fig = (case_idx-1)*6;
    
    % Figure 1: Time Domain and PSD
    figure(base_fig + 1);
    
    % Subplot 1: Time Domain of TRANSMITTED signal (first 3 symbols)
    subplot(2,1,1);
    samples_per_symbol = results.samples_per_symbol;
    num_samples_to_plot = 3 * samples_per_symbol;
    plot(real(txSig(1:min(num_samples_to_plot, length(txSig)))));
    title([channel_type ' - Time Domain Transmitted Signal (First 3 OFDM Symbols)']);
    xlabel('Samples'); ylabel('Amplitude'); grid on;
    
    % Add vertical lines to show symbol boundaries
    hold on;
    for i = 1:2
        xline(i * samples_per_symbol, 'r--', 'LineWidth', 1);
    end
    hold off;
    legend('Signal', 'Symbol Boundaries', 'Location', 'best');
    
    % Subplot 2: PSD of transmitted signal (before channel effects)
    subplot(2,1,2);
    pwelch(txSig, [], [], [], fs, 'centered', 'power');
    title([channel_type ' - PSD of Transmitted OFDM Signal (Before Channel)']);
    
    % Figure 2: EVM vs SNR
    figure(base_fig + 2);
    plot(results.SNR_dB_vec, results.rmsEVM_vec_pilot, 'b-o', 'LineWidth', 1.5, 'DisplayName', 'Pilot MMSE + CPE');
    hold on;
    plot(results.SNR_dB_vec, results.rmsEVM_vec_preamble, 'r-s', 'LineWidth', 1.5, 'DisplayName', 'Preamble MMSE + CPE');
    plot(results.SNR_dB_vec, results.rmsEVM_vec_noEQ, 'g-^', 'LineWidth', 1.5, 'DisplayName', 'No Equalization');
    hold off;
    xlabel('SNR (dB)');
    title([channel_type ' - EVM Performance']);
    ylabel('RMS EVM (%)'); 
    legend('show', 'Location', 'northeast');
    grid on;
    set(gca, 'YScale', 'log');
    
    % Figure 3: BER vs SNR
    figure(base_fig + 3);
    semilogy(results.SNR_dB_vec, results.ber_vec_pilot, [color marker '-'], 'LineWidth', 1.5, 'DisplayName', 'Pilot MMSE + CPE');
    xlabel('SNR (dB)');
    title([channel_type ' - BER Performance']);
    ylabel('Bit Error Rate (BER)'); 
    legend('show', 'Location', 'southwest'); 
    grid on;
    
    % Figure 4: Constellation Diagrams
    figure(base_fig + 4);
    nplot = min(6000, numel(results.rxVec_noEQ));
    sel = randperm(numel(results.rxVec_noEQ), nplot);
    
    % Get EVM values from the last_EVM structure
    evm_pilot = results.last_EVM.rmsEVM_pilot;
    evm_preamble = results.last_EVM.rmsEVM_preamble;
    evm_noEQ = results.last_EVM.rmsEVM_noEQ;
    
    subplot(1,3,1);
    plot(real(results.rxVec_noEQ(sel)), imag(results.rxVec_noEQ(sel)), '.');
    title(sprintf('%s - No Equalization\nEVM: %.3f%%', channel_type, evm_noEQ)); 
    axis square; grid on;
    
    subplot(1,3,2);
    plot(real(results.eqPilotVec(sel)), imag(results.eqPilotVec(sel)), '.');
    title(sprintf('%s - Pilot MMSE+CPE\nEVM: %.3f%%', channel_type, evm_pilot)); 
    axis square; grid on;
    
    subplot(1,3,3);
    plot(real(results.eqPreambleVec(sel)), imag(results.eqPreambleVec(sel)), '.');
    title(sprintf('%s - Preamble MMSE+CPE\nEVM: %.3f%%', channel_type, evm_preamble)); 
    axis square; grid on;
    
    % Figure 5: Channel Analysis
    figure(base_fig + 5);
    
    % Subplot 1: Power Delay Profile
    subplot(2,1,1);
    pathDelays = [0 0.3e-8 0.8e-8 2.5e-8];
    pathPowerdB = [0 -3 -6 -9];
    stem(pathDelays * 1e6, pathPowerdB, 'LineWidth', 1.5);
    title([channel_type ' - Power Delay Profile']);
    xlabel('Delay (microseconds)'); ylabel('Power (dB)'); grid on;
    
    % Subplot 2: Channel Frequency Response
    subplot(2,1,2);
    if ~isempty(results.H_est_pilot_grid)
        % Use the stored numPreamble from results
        numPreamble_used = results.numPreamble;
        H_freq_slice = results.H_est_pilot_grid(:, numPreamble_used + 1); 
        
        % Create frequency axis
        numSC = results.numSC;
        SCS = results.SCS;
        freq_axis_MHz = ((1:numSC) - numSC/2) * (SCS / 1e6);
        
        plot(freq_axis_MHz, 20*log10(abs(H_freq_slice)));
        title([channel_type ' - Channel Frequency Response']);
        xlabel('Frequency (MHz)'); ylabel('Gain (dB)'); grid on;
        xlim([min(freq_axis_MHz), max(freq_axis_MHz)]);
    else
        text(0.5, 0.5, 'No channel estimate available', 'HorizontalAlignment', 'center');
        title([channel_type ' - Channel Frequency Response']);
        xlim([0 1]); ylim([0 1]);
    end
    
    % Figure 6: PAPR Analysis - NEW FIGURE
    figure(base_fig + 6);
    
    % Subplot 1: PAPR Comparison
    subplot(2,1,1);
    papr_values = [results.PAPR_tx, results.PAPR_rx];
    bar(papr_values);
    set(gca, 'XTickLabel', {'Transmitted', 'Received'});
    ylabel('PAPR (dB)');
    title([channel_type ' - PAPR Comparison']);
    grid on;
    
    % Add PAPR values on top of bars
    text(1, results.PAPR_tx+0.1, sprintf('%.2f dB', results.PAPR_tx), 'HorizontalAlignment', 'center');
    text(2, results.PAPR_rx+0.1, sprintf('%.2f dB', results.PAPR_rx), 'HorizontalAlignment', 'center');
    
    % Subplot 2: CCDF of PAPR
    subplot(2,1,2);
    [ccdf_tx, papr_levels_tx] = calculateCCDF(txSig);
    [ccdf_rx, papr_levels_rx] = calculateCCDF(results.rxSig);
    
    semilogy(papr_levels_tx, ccdf_tx, 'b-', 'LineWidth', 2, 'DisplayName', 'Transmitted');
    hold on;
    semilogy(papr_levels_rx, ccdf_rx, 'r-', 'LineWidth', 2, 'DisplayName', 'Received');
    hold off;
    xlabel('PAPR (dB)');
    ylabel('CCDF (Probability)');
    title([channel_type ' - PAPR CCDF']);
    legend('show', 'Location', 'southwest');
    grid on;
    xlim([0, max([papr_levels_tx, papr_levels_rx])]);
end

%% ---------------- CCDF Calculation Function ----------------
function [ccdf, papr_levels] = calculateCCDF(signal)
    % Calculate Complementary Cumulative Distribution Function (CCDF) of PAPR
    num_points = 1000;
    papr_levels = linspace(0, 20, num_points);
    ccdf = zeros(size(papr_levels));
    
    % Calculate instantaneous PAPR for the entire signal
    power = abs(signal).^2;
    avg_power = mean(power);
    instant_papr = 10*log10(power / avg_power);
    
    for i = 1:length(papr_levels)
        ccdf(i) = mean(instant_papr > papr_levels(i));
    end
end