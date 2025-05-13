% PRSA with increasing SNR over time
clear; clc; close all;

% Sampling parameters
fs = 25;  % Sampling frequency (Hz)
T = 100;  % Total duration (s)
t = 0:1/fs:T; % Time vector

% Generate signal with increasing SNR
base_signal = sin(2*pi*0.2*t) + 0.5*sin(2*pi*0.05*t); % Sinusoidal components
noise = randn(size(t)); % Initial white noise
snr_increase = linspace(0.5, 5, length(t)); % Increasing SNR
signal = base_signal + noise ./ snr_increase; % Composite signal

% Define PRSA parameters
L = 50;  % Half window length
tau = 5; % Delay parameter
theta = 0.2; % Threshold for anchor point selection

% Step 1: Detect anchor points
anchor_points = find(signal(1:end-tau) - signal(tau+1:end) > theta);
N = length(anchor_points);

% Step 2: Apply PRSA on different time segments
num_segments = 3; % Divide signal into parts to analyze PRSA evolution
segment_length = floor(length(t) / num_segments);
prsa_curves = zeros(num_segments, 2*L+1);
prsa_variance = zeros(num_segments, 1); % Array to store variance of PRSA curves

for seg = 1:num_segments
    % Define time segment
    seg_start = (seg-1) * segment_length + 1;
    seg_end = seg * segment_length;
    
    % Find anchor points in this segment
    segment_anchors = anchor_points(anchor_points >= seg_start & anchor_points <= seg_end);
    Nseg = length(segment_anchors);
    
    % Extract windows
    prsa_matrix = nan(Nseg, 2*L+1);
    for i = 1:Nseg
        idx = segment_anchors(i);
        if idx - L > 0 && idx + L <= length(signal)
            prsa_matrix(i, :) = signal(idx - L : idx + L);
        end
    end
    
    % Compute PRSA curve for this segment
    prsa_curves(seg, :) = nanmean(prsa_matrix, 1);
    
    % Calculate variance of PRSA curve for this segment
    prsa_variance(seg) = var(prsa_curves(seg, :), 'omitnan');
end

% Plot results in one figure with two subplots
figure;

% Plot original signal with SNR evolution
subplot(2,1,1);
plot(t, signal, 'k'); hold on;
xlabel('Time (s)');
ylabel('Signal Amplitude');
title('Original Signal with Increasing SNR');
grid on;

% Plot PRSA curves for different time segments
subplot(2,1,2);
k = -L:L;
plot(k/fs, prsa_curves(1, :), 'r', 'LineWidth', 2); hold on;
plot(k/fs, prsa_curves(2, :), 'g', 'LineWidth', 2);
plot(k/fs, prsa_curves(3, :), 'b', 'LineWidth', 2);

% Annotate variance for each PRSA curve
for seg = 1:num_segments
    text(0, prsa_curves(seg, L+1), sprintf('Var = %.3f', prsa_variance(seg)), ...
        'Color', 'black', 'FontSize', 10, 'FontWeight', 'bold');
end

xlabel('Time Relative to Anchor Point (s)');
ylabel('PRSA Amplitude');
title('PRSA Curves for Different Time Segments');
legend('Segment 1 (SNR=0dB)', 'Segment 2 (SNR=5dB)', 'Segment 3 (SNR=10dB)');
grid on;
