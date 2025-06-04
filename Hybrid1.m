% Dynamic Adaptive Filtering with Real-time LMS-RLS Switching
% Continuously compares and switches between LMS and RLS based on SNR performance
clear all; close all; clc;

%% Signal Parameters
fs = 50000; % Sampling frequency (Hz) - reduced for better visualization
N = 800; % Number of samples - reduced for clearer plots
t = (0:N-1)/fs; % Time vector

%% Clean Signal Generation - Audible Sine Wave
% Signal parameters
freq_signal = 200;   % Signal frequency (Hz) - A4 musical note, audible range
Am = 1;             % Signal amplitude

% Clean signal - Sine wave in audible range
d = Am * sin(2*pi*freq_signal*t)'; % Clean signal (column vector)
message = d; % Keep reference for plots

%% Wind Noise Generation using Weibull Distribution
shape = 2.0; % Shape parameter (k)
scale = 1.0; % Scale parameter (lambda)
wind_noise = wblrnd(scale, shape, 1, N);

% Normalize wind noise
wind_noise = wind_noise / max(abs(wind_noise));

% Add wind noise to clean signal
x = d + wind_noise'; % Noisy signal

%% Filter Setup
filterLength = 32;  % Filter length
mu_lms = 0.001;     % LMS step size
forgetting_factor = 0.99; % RLS forgetting factor

% Initialize both filters
lmsFilter = dsp.LMSFilter('Length', filterLength, 'Method', 'LMS', 'StepSize', mu_lms);
rlsFilter = dsp.RLSFilter('Length', filterLength, 'ForgettingFactor', forgetting_factor);

%% Initial Evaluation Phase (First 200 samples)
fprintf('=== Dynamic Adaptive Filtering with Real-time Switching ===\n');
fprintf('Phase 1: Initial evaluation of both filters...\n');

eval_length = min(200, N);
eval_samples = 1:eval_length;

% Test both filters on initial samples
[y_lms_init, e_lms_init, ~] = lmsFilter(wind_noise(eval_samples)', x(eval_samples));
[y_rls_init, e_rls_init] = rlsFilter(wind_noise(eval_samples)', x(eval_samples));

% Calculate initial SNR for both filters
signal_power = mean(d.^2);
lms_init_noise = mean((d(eval_samples) - e_lms_init).^2);
rls_init_noise = mean((d(eval_samples) - e_rls_init).^2);

SNR_lms_init = 10 * log10(signal_power / lms_init_noise);
SNR_rls_init = 10 * log10(signal_power / rls_init_noise);

% Choose initial filter based on performance
if SNR_lms_init >= SNR_rls_init
    current_filter = 'LMS';
    fprintf('✓ Initial selection: LMS (SNR: %.2f dB vs RLS: %.2f dB)\n', SNR_lms_init, SNR_rls_init);
else
    current_filter = 'RLS';
    fprintf('✓ Initial selection: RLS (SNR: %.2f dB vs LMS: %.2f dB)\n', SNR_rls_init, SNR_lms_init);
end

%% Dynamic Switching Phase - Sample-by-sample comparison
fprintf('Phase 2: Real-time performance monitoring and switching...\n');

% Reset filters for full processing
lmsFilter.reset();
rlsFilter.reset();

% Initialize arrays for storing results
adaptive_output = zeros(N, 1);
filter_choice = cell(N, 1);
snr_lms_track = zeros(N, 1);
snr_rls_track = zeros(N, 1);
switch_points = [];

% Parameters for dynamic switching
comparison_window = 50; % Window size for SNR comparison
switch_threshold = 0.5; % Minimum SNR improvement to trigger switch (dB)
min_samples_before_switch = 100; % Minimum samples before allowing first switch

% Process each sample
for i = 1:N
    % Get current sample
    current_noise = wind_noise(i);
    current_input = x(i);
    
    % Process with both filters
    [y_lms_sample, e_lms_sample, ~] = lmsFilter(current_noise, current_input);
    [y_rls_sample, e_rls_sample] = rlsFilter(current_noise, current_input);
    
    % Calculate running SNR every 'comparison_window' samples
    if i >= comparison_window
        % Calculate SNR for recent window
        window_start = max(1, i - comparison_window + 1);
        window_end = i;
        
        % Get stored outputs for window calculation
        if i == comparison_window
            % For first calculation, we need to run filters on the window
            lmsFilter_temp = dsp.LMSFilter('Length', filterLength, 'Method', 'LMS', 'StepSize', mu_lms);
            rlsFilter_temp = dsp.RLSFilter('Length', filterLength, 'ForgettingFactor', forgetting_factor);
            
            [~, e_lms_window, ~] = lmsFilter_temp(wind_noise(1:i)', x(1:i));
            [~, e_rls_window] = rlsFilter_temp(wind_noise(1:i)', x(1:i));
        else
            % Use incremental calculation for efficiency
            e_lms_window = e_lms_sample;  % Current LMS output
            e_rls_window = e_rls_sample;  % Current RLS output
        end
        
        % Calculate window SNR (simplified - using current sample performance)
        lms_window_noise = (d(i) - e_lms_sample)^2;
        rls_window_noise = (d(i) - e_rls_sample)^2;
        
        if lms_window_noise > 0
            snr_lms_current = 10 * log10(signal_power / lms_window_noise);
        else
            snr_lms_current = 100; % Very high SNR for perfect reconstruction
        end
        
        if rls_window_noise > 0
            snr_rls_current = 10 * log10(signal_power / rls_window_noise);
        else
            snr_rls_current = 100; % Very high SNR for perfect reconstruction
        end
        
        snr_lms_track(i) = snr_lms_current;
        snr_rls_track(i) = snr_rls_current;
        
        % Switching decision (only after minimum samples)
        if i >= min_samples_before_switch
            if strcmp(current_filter, 'LMS') && snr_rls_current > snr_lms_current + switch_threshold
                current_filter = 'RLS';
                switch_points = [switch_points, i];
                fprintf('→ Switched to RLS at sample %d (%.3fs): RLS=%.2fdB vs LMS=%.2fdB\n', ...
                        i, t(i), snr_rls_current, snr_lms_current);
            elseif strcmp(current_filter, 'RLS') && snr_lms_current > snr_rls_current + switch_threshold
                current_filter = 'LMS';
                switch_points = [switch_points, i];
                fprintf('→ Switched to LMS at sample %d (%.3fs): LMS=%.2fdB vs RLS=%.2fdB\n', ...
                        i, t(i), snr_lms_current, snr_rls_current);
            end
        end
    else
        % For initial samples before comparison window
        snr_lms_track(i) = SNR_lms_init;
        snr_rls_track(i) = SNR_rls_init;
    end
    
    % Select output based on current filter choice
    if strcmp(current_filter, 'LMS')
        adaptive_output(i) = e_lms_sample;
    else
        adaptive_output(i) = e_rls_sample;
    end
    
    filter_choice{i} = current_filter;
end

%% Performance Analysis
fprintf('\nProcessing complete. Total switches: %d\n', length(switch_points));

% Calculate overall performance metrics
noise_power = mean(wind_noise.^2);
input_snr = 10 * log10(signal_power / noise_power);

final_residual_noise = mean((d - adaptive_output).^2);
final_snr = 10 * log10(signal_power / final_residual_noise);
snr_improvement = final_snr - input_snr;

% Calculate filter usage statistics
lms_usage = sum(strcmp(filter_choice, 'LMS')) / N * 100;
rls_usage = sum(strcmp(filter_choice, 'RLS')) / N * 100;

%% Results Display
fprintf('\n=== Performance Metrics ===\n');
fprintf('Signal Power: %.4f\n', signal_power);
fprintf('Noise Power: %.4f\n', noise_power);
fprintf('Input SNR: %.2f dB\n', input_snr);
fprintf('Final Output SNR: %.2f dB\n', final_snr);
fprintf('SNR Improvement: %.2f dB\n', snr_improvement);
fprintf('\n--- Filter Usage Statistics ---\n');
fprintf('LMS Usage: %.1f%% of samples\n', lms_usage);
fprintf('RLS Usage: %.1f%% of samples\n', rls_usage);
fprintf('Total Switches: %d\n', length(switch_points));
fprintf('Average time between switches: %.3f seconds\n', mean(diff([0, switch_points]) / fs));

%% Visualization
% Main signals comparison
figure('Name', 'Dynamic Adaptive Filtering Results');
subplot(4,1,1);
plot(t, wind_noise);
title(sprintf('Wind Noise (Weibull) - Power: %.4f', noise_power));
xlabel('Time (s)'); ylabel('Amplitude'); grid on;

subplot(4,1,2);
plot(t, d);
title(sprintf('Clean Sine Signal (%.0f Hz) - Power: %.4f', freq_signal, signal_power));
xlabel('Time (s)'); ylabel('Amplitude'); grid on;

subplot(4,1,3);
plot(t, x);
title(sprintf('Noisy Sine Signal - Input SNR: %.2f dB', input_snr));
xlabel('Time (s)'); ylabel('Amplitude'); grid on;

subplot(4,1,4);
plot(t, adaptive_output, 'b');
% Mark switch points
title(sprintf('Dynamic Adaptive Output - SNR: %.2f dB (Improvement: %.2f dB)', final_snr, snr_improvement));
xlabel('Time (s)'); ylabel('Amplitude'); grid on;

% Filter selection timeline
figure('Name', 'Filter Selection and Performance Timeline');
subplot(3,1,1);
filter_numeric = zeros(N,1);
for i = 1:N
    if strcmp(filter_choice{i}, 'LMS')
        filter_numeric(i) = 1;
    else
        filter_numeric(i) = 2;
    end
end
plot(t, filter_numeric, 'b' );% hold on;
% for sp = switch_points
%     xline(t(sp), 'r--', 'Switch', 'LineWidth', 1);
% end
ylim([0.5, 2.5]);
yticks([1, 2]);
yticklabels({'LMS', 'RLS'});
title('Active Filter Selection Over Time');
xlabel('Time (s)'); ylabel('Filter Type'); grid on;

subplot(3,1,2);
plot(t, snr_lms_track, 'b', 'DisplayName', 'LMS SNR'); hold on;
plot(t, snr_rls_track, 'r', 'DisplayName', 'RLS SNR');
% for sp = switch_points
%     xline(t(sp), 'k--', 'Switch', 'LineWidth', 1);
% end
title('Real-time SNR Comparison');
xlabel('Time (s)'); ylabel('SNR (dB)'); 
legend('show'); grid on;

subplot(3,1,3);
plot(t, (d - adaptive_output).^2, 'g'); %hold on;
% for sp = switch_points
%     xline(t(sp), 'r--', 'Switch', 'LineWidth', 1);
% end
title('Instantaneous Squared Error');
xlabel('Time (s)'); ylabel('Squared Error'); grid on;

% Frequency domain analysis
figure('Name', 'Frequency Domain Analysis');
Y_noise = fft(wind_noise);
Y_clean = fft(d);
Y_adaptive = fft(adaptive_output);
f = (0:N-1)*(fs/N);

subplot(3,1,1);
plot(f(1:N/2), abs(Y_noise(1:N/2))/N);
title('Wind Noise Frequency Spectrum');
xlabel('Frequency (Hz)'); ylabel('Magnitude'); grid on;

subplot(3,1,2);
plot(f(1:N/2), abs(Y_clean(1:N/2))/N);
title('Clean Sine Signal Frequency Spectrum');
xlabel('Frequency (Hz)'); ylabel('Magnitude'); grid on;

subplot(3,1,3);
plot(f(1:N/2), abs(Y_adaptive(1:N/2))/N);
title('Dynamic Adaptive Filtered Signal Frequency Spectrum');
xlabel('Time (s)'); ylabel('Magnitude'); grid on;

% Performance summary
fprintf('\n=== Dynamic Filtering Analysis ===\n');
if snr_improvement > 2
    fprintf('✓ Dynamic adaptive filtering successful\n');
    fprintf('  Significant SNR improvement: %.2f dB\n', snr_improvement);
    fprintf('  Effective filter switching strategy\n');
elseif snr_improvement > 0
    fprintf('✓ Dynamic adaptive filtering working\n');
    fprintf('  Moderate SNR improvement: %.2f dB\n', snr_improvement);
else
    fprintf('✗ Dynamic filtering needs optimization\n');
    fprintf('  Consider adjusting:\n');
    fprintf('  - Switch threshold (current: %.1f dB)\n', switch_threshold);
    fprintf('  - Comparison window (current: %d samples)\n', comparison_window);
    fprintf('  - Filter parameters\n');
end

fprintf('\n=== Analysis Complete ===\n');
fprintf('Dynamic adaptive filtering completed with %d filter switches\n', length(switch_points));
fprintf('Primary filter used: %s (%.1f%% of time)\n', ...
        char(filter_choice{end}), max(lms_usage, rls_usage));
fprintf('Check the generated plots for detailed visual analysis.\n');
