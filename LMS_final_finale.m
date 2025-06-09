% 3D SNR Analysis for LMS Filter with Different Filter Orders and Step Sizes
% Using WAV files for both desired (“clean”) signal and wind noise
clear all; close all; clc;

%% 1. Specify and read WAV files
desiredFilename   = "D:\MSRIT\Mini Project\Data sets\Desired Signals\piano_2_Cn_n_m_34.wav";    % e.g., clean speech or tone
windNoiseFilename = "D:\MSRIT\Mini Project\Data sets\Wind noises\033_009.wav";  % e.g., recorded wind noise

[d_raw,  fs_d ] = audioread(desiredFilename);    
[wn_raw, fs_wn] = audioread(windNoiseFilename);  

% Convert to mono if stereo
if size(d_raw,2) > 1
    d_raw = mean(d_raw, 2);
end
if size(wn_raw,2) > 1
    wn_raw = mean(wn_raw, 2);
end

% Ensure sampling rates match
if fs_d ~= fs_wn
    error('Sampling rates do not match: %d Hz vs %d Hz. Please resample one file.', fs_d, fs_wn);
end
fs = fs_d;

% Truncate both signals to the same length N (use the shorter length)
len_desired = length(d_raw);
len_noise   = length(wn_raw);
N = min(len_desired, len_noise);

d  = d_raw(1:N);
wn = wn_raw(1:N);

t = (0:N-1)/fs;  % Time vector

%% 2. Normalize wind noise (optional: scale for target SNR)
wn = wn / max(abs(wn));  

% % If you want a specific input SNR, uncomment and set targetSNR_dB
% targetSNR_dB = 0;  % desired input SNR in dB
% signalPower = mean(d.^2);
% noisePower  = mean(wn.^2);
% currentSNR  = 10*log10(signalPower / noisePower);
% scaleFactor = 10^((currentSNR - targetSNR_dB)/20);
% wn = wn * scaleFactor;

noisy_signal = d + wn;

%% 3. Compute Input SNR
signal_power = mean(d.^2);
noise_power  = mean(wn.^2);
SNR_input    = 10 * log10(signal_power / noise_power);

%% 4. Define Parameter Grid for LMS
filter_orders = [4, 6, 8, 12, 14, 18, 20, 24, 32, 64];
step_sizes    = [0.0001, 0.0005, 0.001, 0.005, 0.01, 0.02, 0.05, 0.1];

% Initialize result matrices
SNR_results             = zeros(length(filter_orders), length(step_sizes));
SNR_improvement_results = zeros(length(filter_orders), length(step_sizes));
MSE_results             = zeros(length(filter_orders), length(step_sizes));

fprintf('Starting LMS SNR Analysis with WAV inputs...\n');
total_iterations = length(filter_orders) * length(step_sizes);
current_iteration = 0;

for i = 1:length(filter_orders)
    for j = 1:length(step_sizes)
        current_iteration = current_iteration + 1;
        if mod(current_iteration, 10) == 0 || current_iteration == total_iterations
            fprintf('Progress: %d/%d (%.1f%%)\n', current_iteration, total_iterations, ...
                    100*current_iteration/total_iterations);
        end
        
        filterLength = filter_orders(i);
        mu           = step_sizes(j);
        
        try
            lmsFilter = dsp.LMSFilter('Length', filterLength, 'Method', 'LMS', 'StepSize', mu);
            [~, filtered_signal, ~] = lmsFilter(wn, noisy_signal);
            
            % Compute output SNR
            residual_noise_power = mean((d - filtered_signal).^2);
            SNR_output = 10 * log10(signal_power / residual_noise_power);
            SNR_improvement = SNR_output - SNR_input;
            
            % Compute MSE
            MSE = mean((d - filtered_signal).^2);
            
            % Store results
            SNR_results(i, j)             = SNR_output;
            SNR_improvement_results(i, j) = SNR_improvement;
            MSE_results(i, j)             = MSE;
            
        catch
            fprintf('Warning: LMS unstable at Order=%d, Step=%.4f\n', filterLength, mu);
            SNR_results(i, j)             = SNR_input;
            SNR_improvement_results(i, j) = 0;
            MSE_results(i, j)             = inf;
        end
    end
end

% Find best SNR improvement
[max_improvement, max_idx] = max(SNR_improvement_results(:));
[max_i, max_j] = ind2sub(size(SNR_improvement_results), max_idx);
optimal_order = filter_orders(max_i);
optimal_step  = step_sizes(max_j);

%% 5. Plot SNR Improvement Surface
[X, Y] = meshgrid(step_sizes, filter_orders);

figure;
surf(X, Y, SNR_improvement_results, 'EdgeColor', 'interp');
colormap jet;
colorbar;
xlabel('Step Size (\mu)'); 
ylabel('Filter Order'); 
zlabel('SNR Improvement (dB)');
title('3D SNR Improvement Surface for LMS (WAV inputs)');
set(gca, 'XScale', 'log');
view(45, 30);
hold on;
plot3(optimal_step, optimal_order, max_improvement, 'ro');
text(optimal_step, optimal_order, max_improvement + 0.5, ...
     sprintf('Max SNR: %.2f dB\nOrder: %d, \\mu: %.4f', ...
     max_improvement, optimal_order, optimal_step), 'Color', 'k');

%% 6. Display SNR Result Matrix as Table
disp(' ');
disp('=== SNR Output Matrix (Rows: Filter Orders, Columns: Step Sizes) ===');

SNR_matrix = array2table(SNR_results, ...
    'VariableNames', strcat('mu_', strrep(string(step_sizes), '.', '_')));
SNR_matrix.FilterOrder = filter_orders';
SNR_matrix = movevars(SNR_matrix, 'FilterOrder', 'Before', 1);

disp(SNR_matrix);

%% 7. Summary
disp(' ');
disp('=== LMS Filter Analysis Results ===');
fprintf('Input SNR: %.2f dB\n', SNR_input);
fprintf('Max SNR Improvement: %.2f dB at Order=%d, Step=%.4f\n', ...
        max_improvement, optimal_order, optimal_step);

fprintf('\nPerformance by Filter Order:\n');
for i = 1:length(filter_orders)
    avg_imp = mean(SNR_improvement_results(i, :));
    fprintf('  Order %2d: Avg SNR Improvement = %.2f dB\n', filter_orders(i), avg_imp);
end

fprintf('\nPerformance by Step Size:\n');
for j = 1:length(step_sizes)
    avg_imp = mean(SNR_improvement_results(:, j));
    fprintf('  Step %.4f: Avg SNR Improvement = %.2f dB\n', step_sizes(j), avg_imp);
end