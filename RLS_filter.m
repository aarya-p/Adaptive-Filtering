% RLS Filter Parameter Optimization with 3D Analysis
% Testing different filter orders and forgetting factors
clear all; close all; clc;

% Parameters
fs = 16000;  % Sampling frequency (Hz)
N = 800;    % Number of samples
t = (0:N-1)/fs; % Time vector

% Generate base wind noise using Weibull distribution
rng(12);
shape = 2.0; % Shape parameter (k)
scale = 1.0; % Scale parameter (lambda)
wind_noise = wblrnd(scale, shape, 1, N);

% Normalize wind noise
wind_noise = wind_noise / max(abs(wind_noise));
wind_noise = wind_noise'; % Ensure column vector

% Generate clean signal (440Hz sine wave)
fm = 440;  % Signal frequency
Am = 1;     % Signal amplitude
clean_signal = Am * sin(2*pi*fm*t)'; % Clean signal (column vector)

% Create noisy signal
noisy_signal = clean_signal + wind_noise; % Noisy signal

% Calculate input SNR
signal_power = mean(clean_signal.^2);
noise_power = mean(wind_noise.^2);
SNR_input = 10 * log10(signal_power / noise_power);

fprintf('=== RLS FILTER PARAMETER OPTIMIZATION ===\n');
fprintf('Input SNR: %.2f dB\n', SNR_input);
fprintf('Signal length: %d samples\n', N);
fprintf('Testing different filter orders and forgetting factors...\n\n');

%% Parameter Ranges for Testing
filter_orders = [1,2,3, 4,5, 6,7, 8,9, 10,11, 12,13,14,15 ];  % Different filter lengths
forgetting_factors = [0.95, 0.97, 0.98, 0.99, 0.995, 0.998, 0.999];  % Different forgetting factors
initialInvCov = 1000;  % Keep constant

% Initialize result matrices
num_orders = length(filter_orders);
num_factors = length(forgetting_factors);
SNR_improvements = zeros(num_orders, num_factors);
convergence_rates = zeros(num_orders, num_factors);
MSE_results = zeros(num_orders, num_factors);

% Progress tracking
total_tests = num_orders * num_factors;
test_count = 0;

%% Main Testing Loop
fprintf('Progress: ');
for i = 1:num_orders
    for j = 1:num_factors
        test_count = test_count + 1;
        
        % Current parameters
        filterLength = filter_orders(i);
        forgettingFactor = forgetting_factors(j);
        
        try
            % Create RLS filter
            rlsFilter = dsp.RLSFilter('Length', filterLength, ...
                                     'ForgettingFactor', forgettingFactor, ...
                                     'InitialInverseCovariance', initialInvCov);
            
            % Apply RLS filter
            [noise_estimate, error_signal] = rlsFilter(wind_noise, noisy_signal);
            filtered_signal = noisy_signal - noise_estimate;
            
            % Calculate performance metrics
            residual_noise_power = mean((clean_signal - filtered_signal).^2);
            residual_noise_power = max(residual_noise_power, eps);
            
            SNR_output = 10 * log10(signal_power / residual_noise_power);
            SNR_improvement = SNR_output - SNR_input;
            MSE = mean((clean_signal - filtered_signal).^2);
            
            % Convergence analysis
            squared_error = (clean_signal - filtered_signal).^2;
            window_size = min(50, N/10);
            learning_curve = movmean(squared_error, window_size);
            
            % Multiple convergence detection methods
            convergence_samples = N; % Default to end if no convergence detected
            
            % Method 1: Derivative-based
            if length(learning_curve) > 10
                error_derivative = abs(diff(learning_curve));
                convergence_threshold = 0.01 * mean(learning_curve(end-min(50,length(learning_curve)/4):end));
                convergence_indices = find(error_derivative < convergence_threshold);
                
                if ~isempty(convergence_indices) && convergence_indices(1) > window_size
                    convergence_samples = convergence_indices(1);
                end
            end
            
            % Method 2: Stability-based (more robust)
            if length(learning_curve) > 100 && convergence_samples == N
                final_segment = learning_curve(end-min(50,length(learning_curve)/5):end);
                final_mse = mean(final_segment);
                mse_std = std(final_segment);
                stability_threshold = final_mse + 2*mse_std;
                stable_indices = find(learning_curve <= stability_threshold);
                
                if ~isempty(stable_indices) && stable_indices(1) > window_size/2
                    convergence_samples = stable_indices(1);
                end
            end
            
            % Method 3: Improvement rate threshold
            if convergence_samples == N && length(learning_curve) > 20
                improvement_rate = -diff(learning_curve);
                smoothed_rate = movmean(improvement_rate, min(10, length(improvement_rate)/4));
                rate_threshold = 0.001 * mean(abs(smoothed_rate(1:min(50,end))));
                slow_improvement = find(abs(smoothed_rate) < rate_threshold);
                if ~isempty(slow_improvement)
                    convergence_samples = slow_improvement(1) + 1;
                end
            end
            
            % Ensure convergence sample is within bounds
            convergence_samples = min(convergence_samples, N);
            
            % Calculate convergence rate (samples per filter length - normalized)
            convergence_rate = convergence_samples / filterLength;
            
            % Store results
            SNR_improvements(i, j) = SNR_improvement;
            convergence_rates(i, j) = convergence_rate;
            MSE_results(i, j) = MSE;
            
        catch ME
            % Handle any errors gracefully
            fprintf('\nError with Order=%d, Factor=%.4f: %s\n', filterLength, forgettingFactor, ME.message);
            SNR_improvements(i, j) = NaN;
            convergence_rates(i, j) = NaN;
            MSE_results(i, j) = NaN;
        end
        
        % Progress indicator
        if mod(test_count, 8) == 0
            fprintf('%.0f%% ', 100*test_count/total_tests);
        end
    end
end
fprintf('\nTesting completed!\n\n');

%% Find Best Configurations
% Remove NaN values for analysis
valid_mask = ~isnan(SNR_improvements) & ~isnan(convergence_rates);
valid_snr = SNR_improvements(valid_mask);
valid_conv = convergence_rates(valid_mask);

if any(valid_mask(:))
    % Find best SNR improvement
    [max_snr, max_snr_idx] = max(SNR_improvements(:));
    [best_snr_order_idx, best_snr_factor_idx] = ind2sub(size(SNR_improvements), max_snr_idx);
    
    % Find fastest convergence (minimum convergence rate)
    [min_conv, min_conv_idx] = min(convergence_rates(:));
    [best_conv_order_idx, best_conv_factor_idx] = ind2sub(size(convergence_rates), min_conv_idx);
    
    % Find best overall performance (combination of SNR and convergence)
    % Normalize both metrics and create combined score
    snr_normalized = (SNR_improvements - min(valid_snr)) / (max(valid_snr) - min(valid_snr));
    conv_normalized = (max(valid_conv) - convergence_rates) / (max(valid_conv) - min(valid_conv)); % Invert for convergence
    combined_score = 0.7 * snr_normalized + 0.3 * conv_normalized; % Weight SNR more heavily
    
    [max_combined, max_combined_idx] = max(combined_score(:));
    [best_combined_order_idx, best_combined_factor_idx] = ind2sub(size(combined_score), max_combined_idx);
    
    fprintf('=== OPTIMIZATION RESULTS ===\n');
    fprintf('Best SNR Improvement: %.2f dB\n', max_snr);
    fprintf('  Filter Order: %d, Forgetting Factor: %.4f\n', ...
            filter_orders(best_snr_order_idx), forgetting_factors(best_snr_factor_idx));
    fprintf('  Convergence Rate: %.1f samples/tap\n', convergence_rates(best_snr_order_idx, best_snr_factor_idx));
    
    fprintf('\nFastest Convergence: %.1f samples/tap\n', min_conv);
    fprintf('  Filter Order: %d, Forgetting Factor: %.4f\n', ...
            filter_orders(best_conv_order_idx), forgetting_factors(best_conv_factor_idx));
    fprintf('  SNR Improvement: %.2f dB\n', SNR_improvements(best_conv_order_idx, best_conv_factor_idx));
    
    fprintf('\nBest Overall Performance (70%% SNR + 30%% Convergence):\n');
    fprintf('  Filter Order: %d, Forgetting Factor: %.4f\n', ...
            filter_orders(best_combined_order_idx), forgetting_factors(best_combined_factor_idx));
    fprintf('  SNR Improvement: %.2f dB\n', SNR_improvements(best_combined_order_idx, best_combined_factor_idx));
    fprintf('  Convergence Rate: %.1f samples/tap\n', convergence_rates(best_combined_order_idx, best_combined_factor_idx));
else
    fprintf('No valid results found!\n');
    return;
end

%% 3D Visualization
% Create meshgrid for 3D plotting
[X, Y] = meshgrid(forgetting_factors, filter_orders);

% Figure 1: SNR Improvement Surface
figure('Name', '3D Analysis: SNR Improvement', 'Position', [100, 100, 1200, 500]);

subplot(1, 2, 1);
surf(X, Y, SNR_improvements, 'EdgeColor', 'none', 'FaceAlpha', 0.8);
hold on;

% Mark best points
scatter3(forgetting_factors(best_snr_factor_idx), filter_orders(best_snr_order_idx), ...
         max_snr, 200, 'r', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 2);
scatter3(forgetting_factors(best_combined_factor_idx), filter_orders(best_combined_order_idx), ...
         SNR_improvements(best_combined_order_idx, best_combined_factor_idx), ...
         150, 'g', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 2);

xlabel('Forgetting Factor');
ylabel('Filter Order');
zlabel('SNR Improvement (dB)');
title('SNR Improvement vs Filter Parameters');
colorbar;
grid on;
legend('SNR Surface', 'Best SNR', 'Best Overall', 'Location', 'best');
view(45, 30);

% Figure 2: Convergence Rate Surface
subplot(1, 2, 2);
surf(X, Y, convergence_rates, 'EdgeColor', 'none', 'FaceAlpha', 0.8);
hold on;

% Mark best points
scatter3(forgetting_factors(best_conv_factor_idx), filter_orders(best_conv_order_idx), ...
         min_conv, 200, 'b', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 2);
scatter3(forgetting_factors(best_combined_factor_idx), filter_orders(best_combined_order_idx), ...
         convergence_rates(best_combined_order_idx, best_combined_factor_idx), ...
         150, 'g', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 2);

xlabel('Forgetting Factor');
ylabel('Filter Order');
zlabel('Convergence Rate (samples/tap)');
title('Convergence Rate vs Filter Parameters');
colorbar;
grid on;
legend('Convergence Surface', 'Fastest Convergence', 'Best Overall', 'Location', 'best');
view(45, 30);

% Figure 3: Combined Performance Analysis
figure('Name', '3D Analysis: Combined Performance', 'Position', [200, 200, 1200, 900]);

% Combined score surface
subplot(2, 2, 1);
surf(X, Y, combined_score, 'EdgeColor', 'none', 'FaceAlpha', 0.8);
hold on;
scatter3(forgetting_factors(best_combined_factor_idx), filter_orders(best_combined_order_idx), ...
         max_combined, 200, 'r', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 2);
xlabel('Forgetting Factor');
ylabel('Filter Order');
zlabel('Combined Score');
title('Combined Performance Score');
colorbar;
grid on;
view(45, 30);

% Contour plots for detailed analysis
subplot(2, 2, 2);
contour(X, Y, SNR_improvements, 20, 'ShowText', 'on');
hold on;
plot(forgetting_factors(best_snr_factor_idx), filter_orders(best_snr_order_idx), ...
     'ro', 'MarkerSize', 10, 'LineWidth', 3);
plot(forgetting_factors(best_combined_factor_idx), filter_orders(best_combined_order_idx), ...
     'gs', 'MarkerSize', 10, 'LineWidth', 3);
xlabel('Forgetting Factor');
ylabel('Filter Order');
title('SNR Improvement Contours (dB)');
legend('SNR Contours', 'Best SNR', 'Best Overall', 'Location', 'best');
grid on;

subplot(2, 2, 3);
contour(X, Y, convergence_rates, 20, 'ShowText', 'on');
hold on;
plot(forgetting_factors(best_conv_factor_idx), filter_orders(best_conv_order_idx), ...
     'bo', 'MarkerSize', 10, 'LineWidth', 3);
plot(forgetting_factors(best_combined_factor_idx), filter_orders(best_combined_order_idx), ...
     'gs', 'MarkerSize', 10, 'LineWidth', 3);
xlabel('Forgetting Factor');
ylabel('Filter Order');
title('Convergence Rate Contours (samples/tap)');
legend('Conv. Contours', 'Fastest Conv.', 'Best Overall', 'Location', 'best');
grid on;

% Performance trade-off scatter plot
subplot(2, 2, 4);
scatter(SNR_improvements(:), convergence_rates(:), 50, combined_score(:), 'filled');
hold on;
scatter(max_snr, convergence_rates(best_snr_order_idx, best_snr_factor_idx), ...
        200, 'r', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 2);
scatter(SNR_improvements(best_conv_order_idx, best_conv_factor_idx), min_conv, ...
        200, 'b', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 2);
scatter(SNR_improvements(best_combined_order_idx, best_combined_factor_idx), ...
        convergence_rates(best_combined_order_idx, best_combined_factor_idx), ...
        200, 'g', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 2);
xlabel('SNR Improvement (dB)');
ylabel('Convergence Rate (samples/tap)');
title('Performance Trade-off Analysis');
colorbar;
legend('All Configs', 'Best SNR', 'Fastest Conv.', 'Best Overall', 'Location', 'best');
grid on;

%% Performance Summary Table
fprintf('\n=== PERFORMANCE SUMMARY TABLE ===\n');
fprintf('Order\tFactor\tSNR(dB)\tConv.Rate\tMSE\t\tScore\n');
fprintf('-----\t------\t-------\t---------\t-------\t\t-----\n');

% Show top 5 configurations
[~, sorted_idx] = sort(combined_score(:), 'descend');
for k = 1:min(5, length(sorted_idx))
    [row, col] = ind2sub(size(combined_score), sorted_idx(k));
    if ~isnan(combined_score(row, col))
        fprintf('%d\t%.4f\t%.2f\t\t%.1f\t\t%.6f\t%.3f\n', ...
                filter_orders(row), forgetting_factors(col), ...
                SNR_improvements(row, col), convergence_rates(row, col), ...
                MSE_results(row, col), combined_score(row, col));
    end
end

%% Combined Results Matrix
fprintf('\n=== COMBINED RESULTS MATRIX ===\n');
fprintf('Matrix format: SNR_improvement(dB) / Convergence_rate(samples/tap)\n');
fprintf('Rows: Filter Orders, Columns: Forgetting Factors\n\n');

% Create header with forgetting factors
fprintf('\t\t');
for j = 1:num_factors
    fprintf('λ=%.3f\t\t', forgetting_factors(j));
end
fprintf('\n');

% Print separator line
fprintf('\t\t');
for j = 1:num_factors
    fprintf('--------\t\t');
end
fprintf('\n');

% Print each row (filter order) with combined SNR/Convergence data
for i = 1:num_orders
    fprintf('Order=%d\t\t', filter_orders(i));
    for j = 1:num_factors
        if ~isnan(SNR_improvements(i, j)) && ~isnan(convergence_rates(i, j))
            % Mark best configurations with special symbols
            marker = '';
            if i == best_snr_order_idx && j == best_snr_factor_idx
                marker = '*';  % Best SNR
            elseif i == best_conv_order_idx && j == best_conv_factor_idx
                marker = '+';  % Fastest convergence
            elseif i == best_combined_order_idx && j == best_combined_factor_idx
                marker = '◆';  % Best overall
            end
            
            fprintf('%.2f/%.1f%s\t\t', SNR_improvements(i, j), convergence_rates(i, j), marker);
        else
            fprintf('NaN/NaN\t\t');
        end
    end
    fprintf('\n');
end

fprintf('\nLegend:\n');
fprintf('* = Best SNR improvement\n');
fprintf('+ = Fastest convergence\n');
fprintf('◆ = Best overall performance\n');
fprintf('Format: SNR_improvement(dB) / Convergence_rate(samples/tap)\n');

%% Create Visual Matrix Plot
figure('Name', 'Combined Performance Matrix Visualization', 'Position', [400, 100, 1400, 800]);

% Create combined text matrix for visualization
combined_text = cell(num_orders, num_factors);
color_matrix_snr = SNR_improvements;
color_matrix_conv = convergence_rates;

for i = 1:num_orders
    for j = 1:num_factors
        if ~isnan(SNR_improvements(i, j)) && ~isnan(convergence_rates(i, j))
            % Create text with both values
            marker = '';
            if i == best_snr_order_idx && j == best_snr_factor_idx
                marker = ' *';
            elseif i == best_conv_order_idx && j == best_conv_factor_idx
                marker = ' +';
            elseif i == best_combined_order_idx && j == best_combined_factor_idx
                marker = ' ◆';
            end
            combined_text{i, j} = sprintf('%.2f/%.1f%s', SNR_improvements(i, j), convergence_rates(i, j), marker);
        else
            combined_text{i, j} = 'NaN';
        end
    end
end

% Subplot 1: SNR Improvement Matrix
subplot(2, 2, 1);
imagesc(SNR_improvements);
colorbar;
title('SNR Improvement Matrix (dB)');
xlabel('Forgetting Factor Index');
ylabel('Filter Order Index');

% Add text annotations
for i = 1:num_orders
    for j = 1:num_factors
        if ~isnan(SNR_improvements(i, j))
            text(j, i, sprintf('%.2f', SNR_improvements(i, j)), ...
                 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                 'FontSize', 8, 'Color', 'white', 'FontWeight', 'bold');
        end
    end
end

% Set axis labels
set(gca, 'XTick', 1:num_factors, 'XTickLabel', arrayfun(@(x) sprintf('%.3f', x), forgetting_factors, 'UniformOutput', false));
set(gca, 'YTick', 1:num_orders, 'YTickLabel', arrayfun(@(x) sprintf('%d', x), filter_orders, 'UniformOutput', false));
xtickangle(45);

% Subplot 2: Convergence Rate Matrix
subplot(2, 2, 2);
imagesc(convergence_rates);
colorbar;
title('Convergence Rate Matrix (samples/tap)');
xlabel('Forgetting Factor Index');
ylabel('Filter Order Index');

% Add text annotations
for i = 1:num_orders
    for j = 1:num_factors
        if ~isnan(convergence_rates(i, j))
            text(j, i, sprintf('%.1f', convergence_rates(i, j)), ...
                 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                 'FontSize', 8, 'Color', 'white', 'FontWeight', 'bold');
        end
    end
end

set(gca, 'XTick', 1:num_factors, 'XTickLabel', arrayfun(@(x) sprintf('%.3f', x), forgetting_factors, 'UniformOutput', false));
set(gca, 'YTick', 1:num_orders, 'YTickLabel', arrayfun(@(x) sprintf('%d', x), filter_orders, 'UniformOutput', false));
xtickangle(45);

% Subplot 3: Combined Score Matrix
subplot(2, 2, 3);
imagesc(combined_score);
colorbar;
title('Combined Performance Score Matrix');
xlabel('Forgetting Factor Index');
ylabel('Filter Order Index');

% Add text annotations and mark best points
for i = 1:num_orders
    for j = 1:num_factors
        if ~isnan(combined_score(i, j))
            % Determine text color based on background
            if combined_score(i, j) > 0.5
                text_color = 'white';
            else
                text_color = 'black';
            end
            
            marker = '';
            if i == best_combined_order_idx && j == best_combined_factor_idx
                marker = ' ◆';
                text_color = 'red';
            end
            
            text(j, i, sprintf('%.3f%s', combined_score(i, j), marker), ...
                 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                 'FontSize', 8, 'Color', text_color, 'FontWeight', 'bold');
        end
    end
end

set(gca, 'XTick', 1:num_factors, 'XTickLabel', arrayfun(@(x) sprintf('%.3f', x), forgetting_factors, 'UniformOutput', false));
set(gca, 'YTick', 1:num_orders, 'YTickLabel', arrayfun(@(x) sprintf('%d', x), filter_orders, 'UniformOutput', false));
xtickangle(45);

% Subplot 4: Combined Text Matrix (SNR/Convergence)
subplot(2, 2, 4);
% Create a custom visualization for the combined text
axis([0.5 num_factors+0.5 0.5 num_orders+0.5]);
set(gca, 'YDir', 'reverse');

% Draw grid
for i = 1:num_orders+1
    line([0.5 num_factors+0.5], [i-0.5 i-0.5], 'Color', 'k', 'LineWidth', 0.5);
end
for j = 1:num_factors+1
    line([j-0.5 j-0.5], [0.5 num_orders+0.5], 'Color', 'k', 'LineWidth', 0.5);
end

% Add combined text
for i = 1:num_orders
    for j = 1:num_factors
        if ~isnan(SNR_improvements(i, j)) && ~isnan(convergence_rates(i, j))
            % Color code based on combined score
            if ~isnan(combined_score(i, j))
                if combined_score(i, j) > 0.7
                    bg_color = [0.2 0.8 0.2]; % Green for high score
                elseif combined_score(i, j) > 0.4
                    bg_color = [1 1 0.2]; % Yellow for medium score
                else
                    bg_color = [1 0.7 0.7]; % Light red for low score
                end
                
                % Draw background rectangle
                rectangle('Position', [j-0.4, i-0.4, 0.8, 0.8], ...
                         'FaceColor', bg_color, 'EdgeColor', 'none');
            end
            
            % Add text
            marker = '';
            text_color = 'black';
            if i == best_snr_order_idx && j == best_snr_factor_idx
                marker = ' *';
                text_color = 'red';
            elseif i == best_conv_order_idx && j == best_conv_factor_idx
                marker = ' +';
                text_color = 'blue';
            elseif i == best_combined_order_idx && j == best_combined_factor_idx
                marker = ' ◆';
                text_color = 'red';
            end
            
            text(j, i, sprintf('%.2f/%.1f%s', SNR_improvements(i, j), convergence_rates(i, j), marker), ...
                 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                 'FontSize', 7, 'Color', text_color, 'FontWeight', 'bold');
        else
            text(j, i, 'NaN', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                 'FontSize', 8, 'Color', 'gray');
        end
    end
end

title('Combined Matrix: SNR(dB)/Conv.Rate(samp/tap)');
xlabel('Forgetting Factor');
ylabel('Filter Order');

% Set axis labels
set(gca, 'XTick', 1:num_factors, 'XTickLabel', arrayfun(@(x) sprintf('%.3f', x), forgetting_factors, 'UniformOutput', false));
set(gca, 'YTick', 1:num_orders, 'YTickLabel', arrayfun(@(x) sprintf('%d', x), filter_orders, 'UniformOutput', false));
xtickangle(45);

% Add legend
text(num_factors+0.7, 1, '* Best SNR', 'Color', 'red', 'FontSize', 8);
text(num_factors+0.7, 2, '+ Fastest Conv.', 'Color', 'blue', 'FontSize', 8);
text(num_factors+0.7, 3, '◆ Best Overall', 'Color', 'red', 'FontSize', 8);
text(num_factors+0.7, 5, 'Green: High Score', 'BackgroundColor', [0.2 0.8 0.2], 'FontSize', 8);
text(num_factors+0.7, 6, 'Yellow: Med Score', 'BackgroundColor', [1 1 0.2], 'FontSize', 8);
text(num_factors+0.7, 7, 'Pink: Low Score', 'BackgroundColor', [1 0.7 0.7], 'FontSize', 8);

%% Detailed Analysis of Best Configuration
fprintf('\n=== DETAILED ANALYSIS OF BEST CONFIGURATION ===\n');
best_order = filter_orders(best_combined_order_idx);
best_factor = forgetting_factors(best_combined_factor_idx);

% Re-run with best parameters for detailed analysis
rlsFilter_best = dsp.RLSFilter('Length', best_order, ...
                              'ForgettingFactor', best_factor, ...
                              'InitialInverseCovariance', initialInvCov);

[noise_estimate_best, ~] = rlsFilter_best(wind_noise, noisy_signal);
filtered_signal_best = noisy_signal - noise_estimate_best;

% Detailed performance metrics
residual_noise_power_best = mean((clean_signal - filtered_signal_best).^2);
SNR_output_best = 10 * log10(signal_power / residual_noise_power_best);
SNR_improvement_best = SNR_output_best - SNR_input;
MSE_best = mean((clean_signal - filtered_signal_best).^2);

fprintf('Best Configuration: Order=%d, Factor=%.4f\n', best_order, best_factor);
fprintf('Input SNR: %.2f dB\n', SNR_input);
fprintf('Output SNR: %.2f dB\n', SNR_output_best);
fprintf('SNR Improvement: %.2f dB\n', SNR_improvement_best);
fprintf('MSE: %.6f\n', MSE_best);
fprintf('Convergence Rate: %.1f samples/tap\n', convergence_rates(best_combined_order_idx, best_combined_factor_idx));

% Signal comparison plot
figure('Name', 'Best Configuration Results', 'Position', [300, 300, 1000, 600]);
subplot(2,1,1);
plot(t, clean_signal, 'g-', 'LineWidth', 2);
hold on;
plot(t, noisy_signal, 'r:', 'LineWidth', 1);
plot(t, filtered_signal_best, 'b-', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Amplitude');
title(sprintf('Best RLS Configuration: Order=%d, Factor=%.4f', best_order, best_factor));
legend('Clean Signal', 'Noisy Signal', 'Filtered Signal', 'Location', 'best');
grid on;

% MSE evolution for best configuration
squared_error_best = (clean_signal - filtered_signal_best).^2;
learning_curve_best = movmean(squared_error_best, min(50, N/10));

subplot(2,1,2);
semilogy(1:length(learning_curve_best), learning_curve_best, 'b-', 'LineWidth', 2);
hold on;
semilogy([1, length(learning_curve_best)], [MSE_best, MSE_best], 'g--', 'LineWidth', 1.5);
convergence_best = convergence_rates(best_combined_order_idx, best_combined_factor_idx) * best_order;
if convergence_best <= length(learning_curve_best)
    semilogy(convergence_best, learning_curve_best(round(convergence_best)), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
    legend('MSE Evolution', 'Final MSE', 'Convergence Point', 'Location', 'best');
else
    legend('MSE Evolution', 'Final MSE', 'Location', 'best');
end
xlabel('Sample');
ylabel('MSE (log scale)');
title('MSE Convergence for Best Configuration');
grid on;

fprintf('\n=== ANALYSIS COMPLETE ===\n');
fprintf('The 3D plots show the relationship between filter parameters and performance.\n');
fprintf('Red markers indicate best SNR, blue markers indicate fastest convergence,\n');
fprintf('and green markers indicate the best overall balance.\n');
