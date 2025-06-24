% Dynamic Adaptive Filtering with Real‐time LMS‐RLS Switching
% Using WAV files for both desired ("clean") signal and wind noise.

clear all; close all; clc;

%% 1. Specify and read WAV files
desiredFilename   = "D:\MSRIT\Mini Project\Data sets\Desired Signals\piano_2_An_j_p_72.wav";    % e.g., clean speech or tone
windNoiseFilename = "D:\MSRIT\Mini Project\Data sets\Wind noises\000_002.wav";  % e.g., recorded wind noise

[d, fs_d]  = audioread(desiredFilename);    % d: desired signal
[wn, fs_wn] = audioread(windNoiseFilename);  % wn: wind noise

% If stereo, convert to mono by averaging channels
if size(d,2) > 1
    d = mean(d, 2);
end
if size(wn,2) > 1
    wn = mean(wn, 2);
end

% Ensure sampling rates match
if fs_d ~= fs_wn
    error('Sampling rates do not match: %d Hz vs %d Hz. Please resample one file.', fs_d, fs_wn);
end
fs = fs_d/4;  % common sampling frequency

% Truncate (or zero‐pad) both signals to the same length N
len_desired = length(d);
len_noise   = length(wn);
N = min(len_desired, len_noise);

% Option 1: Truncate both to the shorter length:
d  = d(1:N);
wn = wn(1:N);

% % Option 2: If you prefer to keep the entire desired file, zero‐pad the shorter:
% if len_desired > len_noise
%     wn = [wn; zeros(len_desired - len_noise, 1)];
%     N = len_desired;
% elseif len_noise > len_desired
%     d = [d; zeros(len_noise - len_desired, 1)];
%     N = len_noise;
% end

t = (0:N-1) / fs;  % Time vector

%% 2. Form the noisy input
% If your wind‐noise WAV is already scaled, you can skip normalization.
% Otherwise, normalize wn to peak = 1 so that you can adjust SNR manually.
wn = wn / max(abs(wn));  

% Adjust wind‐noise amplitude if you want a specific input SNR:
% targetSNR_dB = 0;  % e.g., 0 dB input SNR
% signalPower = mean(d.^2);
% noisePower  = mean(wn.^2);
% currentSNR  = 10*log10(signalPower / noisePower);
% scaleFactor = 10^((currentSNR - targetSNR_dB)/20);
% wn = wn * scaleFactor;

x = d + wn;  % noisy observation

%% 3. Filter Setup
filterLength_LMS = 1;
filterLength_RLS     = 3;       % Filter length (number of taps)
mu_lms           = 0.1;    % LMS step size
forgetting_factor = 0.999;    % RLS forgetting factor

% Create DSP filter objects
lmsFilter = dsp.LMSFilter('Length', filterLength_LMS, ...
                          'Method', 'LMS', ...
                          'StepSize', mu_lms);
rlsFilter = dsp.RLSFilter('Length', filterLength_RLS, ...
                          'ForgettingFactor', forgetting_factor);

%% 4. Initial Evaluation Phase (First 200 samples)
fprintf('=== Dynamic Adaptive Filtering with Real‐time Switching ===\n');
fprintf('Phase 1: Initial evaluation of both filters...\n');

eval_length     = min(200, N);
eval_samples    = 1:eval_length;

% Run LMS and RLS on the first eval_length samples
[y_lms_init, e_lms_init, ~] = lmsFilter( wn(eval_samples), x(eval_samples) );
[y_rls_init, e_rls_init]     = rlsFilter( wn(eval_samples), x(eval_samples) );

% Compute actual desired signal power (over entire d)
signal_power = mean(d.^2);

% Estimate noise power after each filter
lms_init_noise = mean( (d(eval_samples) - e_lms_init).^2 );
rls_init_noise = mean( (d(eval_samples) - e_rls_init).^2 );

SNR_lms_init = 10 * log10(signal_power / lms_init_noise);
SNR_rls_init = 10 * log10(signal_power / rls_init_noise);

% Pick initial filter
if SNR_lms_init >= SNR_rls_init
    current_filter = 'LMS';
    fprintf('✓ Initial selection: LMS (SNR: %.2f dB vs RLS: %.2f dB)\n', ...
            SNR_lms_init, SNR_rls_init);
else
    current_filter = 'RLS';
    fprintf('✓ Initial selection: RLS (SNR: %.2f dB vs LMS: %.2f dB)\n', ...
            SNR_rls_init, SNR_lms_init);
end

%% 5. Dynamic Switching Phase
fprintf('Phase 2: Real‐time performance monitoring and switching...\n');

% Reset filters to clear internal states before full‐length run
reset(lmsFilter);
reset(rlsFilter);

adaptive_output = zeros(N, 1);
filter_choice   = cell(N, 1);
snr_lms_track   = zeros(N, 1);
snr_rls_track   = zeros(N, 1);
switch_points   = [];
lms_output = zeros(N,1);
rls_output = zeros(N,1);

% Parameters for dynamic switching
comparison_window        = 50;   % #samples over which to compare SNR (simplified)
switch_threshold         = 0.5;  % dB improvement needed to switch
min_samples_before_switch = 100; % wait at least this many samples before switching

for i = 1:N
    % Current samples
    current_noise = wn(i);
    current_input = x(i);
    
    % Filter outputs and error signals
    [y_lms_sample, e_lms_sample, ~] = lmsFilter(current_noise, current_input);
    [y_rls_sample, e_rls_sample]     = rlsFilter(current_noise, current_input);

    % Store each filter's output
    lms_output(i) = e_lms_sample;
    rls_output(i) = e_rls_sample;
    
    % Once we have at least comparison_window samples, compute a "running" SNR
    if i >= comparison_window
        % Simplified: use instantaneous squared error as "noise" measure
        lms_window_noise = (d(i) - e_lms_sample)^2;
        rls_window_noise = (d(i) - e_rls_sample)^2;
        
        if lms_window_noise > 0
            snr_lms_current = 10 * log10(signal_power / lms_window_noise);
        else
            snr_lms_current = 100; % very high if no error
        end
        
        if rls_window_noise > 0
            snr_rls_current = 10 * log10(signal_power / rls_window_noise);
        else
            snr_rls_current = 100;
        end
        
        snr_lms_track(i) = snr_lms_current;
        snr_rls_track(i) = snr_rls_current;
        
        if i >= min_samples_before_switch
            if strcmp(current_filter, 'LMS') && (snr_rls_current > snr_lms_current + switch_threshold)
                current_filter = 'RLS';
                switch_points = [switch_points, i];
                fprintf('→ Switched to RLS at sample %d (%.3fs): RLS=%.2fdB vs LMS=%.2fdB\n', ...
                        i, t(i), snr_rls_current, snr_lms_current);
            elseif strcmp(current_filter, 'RLS') && (snr_lms_current > snr_rls_current + switch_threshold)
                current_filter = 'LMS';
                switch_points = [switch_points, i];
                fprintf('→ Switched to LMS at sample %d (%.3fs): LMS=%.2fdB vs RLS=%.2fdB\n', ...
                        i, t(i), snr_lms_current, snr_rls_current);
            end
        end
    else
        % Before the first comparison, just hold initial SNRs
        snr_lms_track(i) = SNR_lms_init;
        snr_rls_track(i) = SNR_rls_init;
    end
    
    % Choose the output from the currently selected filter
    if strcmp(current_filter, 'LMS')
        adaptive_output(i) = e_lms_sample;
    else
        adaptive_output(i) = e_rls_sample;
    end
    
    filter_choice{i} = current_filter;
end

%% 6. Performance Analysis
fprintf('\nProcessing complete. Total switches: %d\n', length(switch_points));

noise_power  = mean(wn.^2);
input_snr    = 10 * log10(signal_power / noise_power);

% Calculate final SNR for all three outputs
final_residual_noise_adaptive = mean( (d - adaptive_output).^2 );
final_residual_noise_lms = mean( (d - lms_output).^2 );
final_residual_noise_rls = mean( (d - rls_output).^2 );

final_snr_adaptive = 10 * log10(signal_power / final_residual_noise_adaptive);
final_snr_lms = 10 * log10(signal_power / final_residual_noise_lms);
final_snr_rls = 10 * log10(signal_power / final_residual_noise_rls);

snr_improvement_adaptive = final_snr_adaptive - input_snr;
snr_improvement_lms = final_snr_lms - input_snr;
snr_improvement_rls = final_snr_rls - input_snr;

%% 6.1. Convergence Analysis (LMS characteristic)
% Calculate Mean Squared Error (MSE) over time for convergence analysis
window_size = 100;
mse_lms = zeros(N, 1);
mse_adaptive = zeros(N, 1);
mse_rls = zeros(N,1);
for i = 1:N
    start_idx = max(1, i - window_size + 1);
    end_idx = i;
    
    % Calculate windowed MSE
    mse_lms(i) = mean((d(start_idx:end_idx) - lms_output(start_idx:end_idx)).^2);
    mse_adaptive(i) = mean((d(start_idx:end_idx) - adaptive_output(start_idx:end_idx)).^2);
    mse_rls(i) = mean((d(start_idx:end_idx) - rls_output(start_idx:end_idx)).^2);
end

% Convergence metrics
convergence_threshold = 0.1 * mse_lms(end); % 10% of final MSE
lms_convergence_time = find(mse_lms <= (mse_lms(end) + convergence_threshold), 1);
adaptive_convergence_time = find(mse_adaptive <= (mse_adaptive(end) + convergence_threshold), 1);

if isempty(lms_convergence_time), lms_convergence_time = N; end
if isempty(adaptive_convergence_time), adaptive_convergence_time = N; end

% Convergence rate (slope of MSE decrease in early samples)
early_samples = min(500, floor(N/4));
if early_samples > 50
    p_lms = polyfit(1:early_samples, log10(mse_lms(1:early_samples))', 1);
    p_adaptive = polyfit(1:early_samples, log10(mse_adaptive(1:early_samples))', 1);
    lms_convergence_rate = -p_lms(1); % Negative slope means faster convergence
    adaptive_convergence_rate = -p_adaptive(1);
else
    lms_convergence_rate = 0;
    adaptive_convergence_rate = 0;
end

lms_usage = sum(strcmp(filter_choice, 'LMS')) / N * 100;
rls_usage = sum(strcmp(filter_choice, 'RLS')) / N * 100;

fprintf('\n=== Performance Metrics ===\n');
fprintf('Signal Power: %.4f\n', signal_power);
fprintf('Noise Power: %.4f\n', noise_power);
fprintf('Input SNR: %.2f dB\n', input_snr);
fprintf('\n--- Final SNR Comparison ---\n');
fprintf('Adaptive Filter SNR: %.2f dB (Improvement: %.2f dB)\n', final_snr_adaptive, snr_improvement_adaptive);
fprintf('LMS Filter SNR: %.2f dB (Improvement: %.2f dB)\n', final_snr_lms, snr_improvement_lms);
fprintf('RLS Filter SNR: %.2f dB (Improvement: %.2f dB)\n', final_snr_rls, snr_improvement_rls);


fprintf('\n--- Combined Performance Analysis ---\n');
fprintf('✓ Adaptive Filter Advantages:\n');
fprintf('  • Fast Convergence (LMS trait): ');
if adaptive_convergence_rate >= lms_convergence_rate * 0.8
    fprintf('ACHIEVED - %.1f%% of LMS convergence speed\n', (adaptive_convergence_rate/lms_convergence_rate)*100);
else
    fprintf('Partial - %.1f%% of LMS convergence speed\n', (adaptive_convergence_rate/lms_convergence_rate)*100);
end

fprintf('\n--- Filter Usage Statistics ---\n');
fprintf('LMS Usage: %.1f%% of samples\n', lms_usage);
fprintf('RLS Usage: %.1f%% of samples\n', rls_usage);
fprintf('Total Switches: %d\n', length(switch_points));
if length(switch_points) > 1
    fprintf('Average time between switches: %.3f seconds\n', mean(diff([0, switch_points]) / fs));
end

%% 7. Visualization
% (You can comment out any figure sections if you don't need plots.)

figure('Name', 'Input Signals');
subplot(2,1,1);
plot(t, wn);
title(sprintf('Wind Noise (from WAV) – Power: %.4f', noise_power));
xlabel('Time (s)'); ylabel('Amplitude'); grid on;

subplot(2,1,2);
plot(t, d);
title(sprintf('Desired Signal (from WAV) – Power: %.4f', signal_power));
xlabel('Time (s)'); ylabel('Amplitude'); grid on;

figure('Name',"Noisy signal and Filtered Signal")
subplot(2,1,1);
plot(t, x);
title(sprintf('Noisy Input: Desired + Wind – Input SNR: %.2f dB', input_snr));
xlabel('Time (s)'); ylabel('Amplitude'); grid on;

subplot(2,1,2);
plot(t, adaptive_output, 'b');
title(sprintf('Adaptive Output – Final SNR: %.2f dB (Improvement: %.2f dB)', final_snr_adaptive, snr_improvement_adaptive));
xlabel('Time (s)'); ylabel('Amplitude'); grid on;


figure('Name', 'Convergence Analysis - LMS vs Adaptive');
semilogy(t, mse_lms, 'r', 'LineWidth', 1, 'DisplayName', 'LMS MSE'); hold on;
semilogy(t, mse_rls, 'g', 'LineWidth', 1, 'DisplayName', 'RLS MSE');
semilogy(t, mse_adaptive, 'b', 'LineWidth', 1, 'DisplayName', 'Adaptive MSE');
xlabel('Time (s)'); ylabel('Mean Squared Error (log scale)'); 
title('Convergence Comparison: LMS vs Adaptive Filter');
legend('show'); grid on; hold off;

% subplot(2,1,2);
% plot(t(1:min(1000,N)), mse_lms(1:min(1000,N)), 'r', 'LineWidth', 1.5, 'DisplayName', 'LMS MSE'); hold on;
% plot(t(1:min(1000,N)), mse_adaptive(1:min(1000,N)), 'b', 'LineWidth', 1.5, 'DisplayName', 'Adaptive MSE');
% xlabel('Time (s)'); ylabel('Mean Squared Error'); 
% title('Early Convergence Detail (First 1000 samples)');
% legend('show'); grid on; hold off;

figure('Name', 'Filter Outputs Comparison');
subplot(3,1,1);
plot(t, lms_output, 'r');
title(sprintf('LMS Filter Output – SNR: %.2f dB (Improvement: %.2f dB)', final_snr_lms, snr_improvement_lms));
xlabel('Time (s)'); ylabel('Amplitude'); grid on;

subplot(3,1,2);
plot(t, rls_output, 'g');
title(sprintf('RLS Filter Output – SNR: %.2f dB (Improvement: %.2f dB)', final_snr_rls, snr_improvement_rls));
xlabel('Time (s)'); ylabel('Amplitude'); grid on;

subplot(3,1,3);
plot(t, adaptive_output, 'b');
title(sprintf('Adaptive Filter Output – SNR: %.2f dB (Improvement: %.2f dB)', final_snr_adaptive, snr_improvement_adaptive));
xlabel('Time (s)'); ylabel('Amplitude'); grid on;

% subplot(4,1,4);
% plot(t, lms_output, 'r', 'DisplayName', 'LMS'); hold on;
% plot(t, rls_output, 'g', 'DisplayName', 'RLS');
% plot(t, adaptive_output, 'b', 'DisplayName', 'Adaptive');
% legend('show');
% title('All Filter Outputs Comparison');
% xlabel('Time (s)'); ylabel('Amplitude'); grid on; hold off;

figure('Name', 'Filter Selection and Performance');
subplot(3,1,1);
filter_numeric = zeros(N,1);
for i = 1:N
    if strcmp(filter_choice{i}, 'LMS'), filter_numeric(i) = 1;
    else                         filter_numeric(i) = 2; end
end
plot(t, filter_numeric, 'LineWidth', 1.2);
ylim([0.5, 2.5]);
yticks([1, 2]);
yticklabels({'LMS','RLS'});
title('Active Filter Selection Over Time');
xlabel('Time (s)'); ylabel('Filter Type'); grid on;

subplot(3,1,2);
plot(t, snr_lms_track, 'b', 'DisplayName','LMS SNR'); hold on;
plot(t, snr_rls_track, 'r', 'DisplayName','RLS SNR');
legend('show');
title('Real‐time SNR Comparison');
xlabel('Time (s)'); ylabel('SNR (dB)'); grid on;

subplot(3,1,3);
plot(t, (d - adaptive_output).^2, 'g');
title('Instantaneous Squared Error');
xlabel('Time (s)'); ylabel('Squared Error'); grid on;

figure('Name', 'Frequency Domain Analysis');
Y_noise    = fft(wn);
Y_desired  = fft(d);
Y_adaptive = fft(adaptive_output);
f = (0:N-1) * (fs / N);

subplot(3,1,1);
plot(f(1:floor(N/2)), abs(Y_noise(1:floor(N/2)))/N);
title('Wind Noise Spectrum');
xlabel('Frequency (Hz)'); ylabel('Magnitude'); grid on;

subplot(3,1,2);
plot(f(1:floor(N/2)), abs(Y_desired(1:floor(N/2)))/N);
title('Desired Signal Spectrum');
xlabel('Frequency (Hz)'); ylabel('Magnitude'); grid on;

subplot(3,1,3);
plot(f(1:floor(N/2)), abs(Y_adaptive(1:floor(N/2)))/N);
title('Adaptive Filtered Output Spectrum');
xlabel('Frequency (Hz)'); ylabel('Magnitude'); grid on;

%% 8. Final Remarks
if snr_improvement_adaptive > 2
    fprintf('\n✓ Dynamic adaptive filtering successful – SNR improved by %.2f dB\n', snr_improvement_adaptive);
elseif snr_improvement_adaptive > 0
    fprintf('\n✓ Some SNR improvement: %.2f dB – consider tweaking parameters for more gain\n', snr_improvement_adaptive);
else
    fprintf('\n✗ No SNR improvement – you may need to adjust step sizes, thresholds, or check WAV files\n');
end

fprintf('\n=== Analysis Complete ===\n');
fprintf('Total filter switches: %d\n', length(switch_points));
fprintf('Primary filter used at end: %s\n', filter_choice{end});
fprintf('LMS used %.1f%%, RLS used %.1f%% of the time\n', lms_usage, rls_usage);
fprintf('\n=== SNR Performance Summary ===\n');
fprintf('Best performing filter: ');
if final_snr_adaptive >= max(final_snr_lms, final_snr_rls)
    fprintf('Adaptive Filter (%.2f dB)\n', final_snr_adaptive);
elseif final_snr_lms >= final_snr_rls
    fprintf('LMS Filter (%.2f dB)\n', final_snr_lms);
else
    fprintf('RLS Filter (%.2f dB)\n', final_snr_rls);
end