% Permutation Entropy with increasing SNR from -5 dB to 40 dB
clear; clc; close all;

% Sampling parameters
fs = 25;  % Sampling frequency (Hz)
T = 100;  % Total duration (s) - Signal is now 100 seconds long
t = 0:1/fs:T; % Time vector

% Generate base signal (sinusoidal components)
base_signal = sin(2*pi*0.2*t) + 0.5*sin(2*pi*0.05*t); % Sinusoidal components

% Define SNR in dB from -5 to 40 over time
snr_dB = linspace(-10, 30, length(t));  % Linear increase from -5 dB to 40 dB

% Convert SNR from dB to linear scale
snr_linear = 10.^(snr_dB / 10);

% Generate noise and adjust based on the SNR
noise = randn(size(t));  % White noise
signal = base_signal + noise ./ snr_linear;  % Signal with increasing SNR

% Define Permutation Entropy parameters
m = 3;  % Embedding dimension
tau = 5; % Time delay

% Step 1: Compute Permutation Entropy at different time segments
num_segments = 10; % Divide signal into 5 parts to analyze PE evolution
segment_length = floor(length(t) / num_segments);
pe_values = zeros(num_segments, 1);

for seg = 1:num_segments
    % Define time segment
    seg_start = (seg-1) * segment_length + 1;
    seg_end = seg * segment_length;
    
    % Extract segment from the signal
    segment_signal = signal(seg_start:seg_end);
    
    % Calculate Permutation Entropy for the segment
    pe_values(seg) = permutation_entropy(segment_signal, m, tau);
end

% Step 2: Plot results
figure;

% Plot SNR in dB over time
subplot(3,1,1);
plot(t, snr_dB, 'b', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('SNR (dB)');
title('SNR Increase Over Time');
grid on;

% Plot original signal with SNR evolution
subplot(3,1,2);
plot(t, signal, 'k'); hold on;
xlabel('Time (s)');
ylabel('Signal Amplitude');
title('Original Signal with Increasing SNR');
grid on;

% Plot Permutation Entropy values for different time segments
subplot(3,1,3);
bar(1:num_segments, pe_values, 'FaceColor', 'b');
xlabel('Segment');
ylabel('Permutation Entropy');
title('Permutation Entropy for Different Time Segments');
grid on;



% Supporting function to compute permutation entropy
function H = permutation_entropy(signal, m, tau)
    N = length(signal);
    permutations = zeros(N - (m-1)*tau, 1);
    
    for i = 1:N - (m-1)*tau
        % Create embedding vector of length m
        embedding = signal(i:tau:i + (m-1)*tau);
        [~, permutation] = sort(embedding);  % Sort and get permutation
        
        % Create a unique number for each permutation
        perm_str = num2str(permutation');  % Convert to string representation
        permutations(i) = sum(2.^(0:m-1) .* (perm_str - '0')');  % Convert to number
    end
    
    % Calculate probability distribution of permutations
    unique_permutations = unique(permutations);
    p = histcounts(permutations, [unique_permutations; unique_permutations(end) + 1], 'Normalization', 'probability');
    
    % Permutation Entropy formula
    H = -sum(p .* log(p));
end
