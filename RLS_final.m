% 3D SNR Analysis for RLS Filter with Different Filter Orders and Forgetting Factors
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

% Truncate both signals to the same length N (use the shorter of the two)
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

%% 4. Define Parameter Grid for RLS
filter_orders      = [4, 8, 12, 16, 20, 24, 28, 32];
forgetting_factors = [0.90, 0.92, 0.95, 0.97, 0.98, 0.99, 0.995, 0.999];

% Initialize result matrices
SNR_results               = zeros(length(filter_orders), length(forgetting_factors));
SNR_improvement_results   = zeros(length(filter_orders), length(forgetting_factors));
MSE_results               = zeros(length(filter_orders), length(forgetting_factors));
convergence_samples       = zeros(length(filter_orders), length(forgetting_factors));

fprintf('Starting 3D RLS SNR Analysis with WAV inputs...\n');
total_iterations = length(filter_orders) * length(forgetting_factors);
current_iteration = 0;

for i = 1:length(filter_orders)
    for j = 1:length(forgetting_factors)
        current_iteration = current_iteration + 1;
        if mod(current_iteration, 8) == 0 || current_iteration == total_iterations
            fprintf('Progress: %d/%d (%.1f%%)\n', current_iteration, total_iterations, 100*current_iteration/total_iterations);
        end
        
        filterLength = filter_orders(i);
        lambda       = forgetting_factors(j);
        
        try
            rlsFilter = dsp.RLSFilter('Length', filterLength, 'ForgettingFactor', lambda);
            
            % Run RLS on the entire length N
            [noise_estimate, filtered_signal] = rlsFilter(wn, noisy_signal);
            
            % Compute residual noise power and output SNR
            residual_noise_power = mean((d - filtered_signal).^2);
            SNR_output = 10 * log10(signal_power / residual_noise_power);
            SNR_improvement = SNR_output - SNR_input;
            
            % Compute MSE over all samples
            MSE = mean((d - filtered_signal).^2);
            
            % Determine convergence: look at moving‐average squared error
            squared_error   = (d - filtered_signal).^2;
            learning_curve  = movmean(squared_error, 20);
            error_derivative = abs(diff(learning_curve));
            threshold = 0.01 * mean(learning_curve);
            indices = find(error_derivative < threshold);
            conv_samples = N;
            if ~isempty(indices)
                conv_samples = indices(1) + 20;
            end
            
            % Store results
            SNR_results(i, j)             = SNR_output;
            SNR_improvement_results(i, j) = SNR_improvement;
            MSE_results(i, j)             = MSE;
            convergence_samples(i, j)     = conv_samples;
            
        catch
            fprintf('Warning: RLS unstable at Order=%d, Lambda=%.3f\n', filterLength, lambda);
            SNR_results(i, j)             = SNR_input;
            SNR_improvement_results(i, j) = 0;
            MSE_results(i, j)             = inf;
            convergence_samples(i, j)     = N;
        end
    end
end

% Create meshgrid for plotting
[X, Y] = meshgrid(forgetting_factors, filter_orders);

% Find best SNR improvement
[max_improvement, max_idx] = max(SNR_improvement_results(:));
[max_i, max_j] = ind2sub(size(SNR_improvement_results), max_idx);
optimal_order   = filter_orders(max_i);
optimal_lambda  = forgetting_factors(max_j);

% Find fastest convergence
[min_conv, min_idx] = min(convergence_samples(:));
[conv_i, conv_j] = ind2sub(size(convergence_samples), min_idx);
best_conv_order   = filter_orders(conv_i);
best_conv_lambda  = forgetting_factors(conv_j);

%% 5. Plotting Results

% 5.1 SNR Improvement Surface
figure;
surf(X, Y, SNR_improvement_results, 'EdgeColor', 'interp');
colormap jet;
colorbar;
xlabel('\lambda'); 
ylabel('Filter Order'); 
zlabel('SNR Improvement (dB)');
title('3D Surface: SNR Improvement vs. Order & Forgetting Factor');
hold on;
plot3(optimal_lambda, optimal_order, max_improvement, 'ro', 'MarkerFaceColor', 'r');
text(optimal_lambda, optimal_order, max_improvement + 0.5, ...
     sprintf('Max SNR: %.2f dB\nOrder: %d, \\lambda: %.3f', ...
     max_improvement, optimal_order, optimal_lambda), 'Color', 'k');

% 5.2 Convergence Samples Surface
figure;
surf(X, Y, convergence_samples, 'EdgeColor', 'none');
colormap turbo;
colorbar;
xlabel('\lambda'); 
ylabel('Filter Order'); 
zlabel('Convergence Samples');
title('3D Surface: Convergence Samples vs. Order & Forgetting Factor');
hold on;
plot3(best_conv_lambda, best_conv_order, min_conv, 'ko', 'MarkerFaceColor', 'g');
text(best_conv_lambda, best_conv_order, min_conv + 10, ...
     sprintf('Fastest\nOrder: %d, \\lambda: %.3f\nSamples: %d', ...
     best_conv_order, best_conv_lambda, min_conv), 'Color', 'k');

%% 6. Display Best‐Parameter Summary
fprintf('\n=== Best Parameters ===\n');
fprintf('Max SNR Improvement: %.2f dB at Order = %d, Lambda = %.3f\n', ...
        max_improvement, optimal_order, optimal_lambda);
fprintf('Fastest Convergence: %d samples at Order = %d, Lambda = %.3f\n', ...
        min_conv, best_conv_order, best_conv_lambda);

%% 7. Create and Display Combined Result Matrix
CombinedMatrix = cell(length(filter_orders) + 1, length(forgetting_factors) + 1);
CombinedMatrix{1, 1} = 'Order/\lambda';

for j = 1:length(forgetting_factors)
    CombinedMatrix{1, j+1} = sprintf('%.3f', forgetting_factors(j));
end

for i = 1:length(filter_orders)
    CombinedMatrix{i+1, 1} = filter_orders(i);
    for j = 1:length(forgetting_factors)
        CombinedMatrix{i+1, j+1} = sprintf('%.2f dB | %d', ...
                                           SNR_results(i, j), ...
                                           convergence_samples(i, j));
    end
end

fprintf('\n=== SNR and Convergence Matrix ===\n');
disp(CombinedMatrix);