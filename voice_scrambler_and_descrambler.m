
%  VOICE SCRAMBLER AND DESCRAMBLER SYSTEM
%  Signals and Systems Course Project
%  Technique: Frequency Inversion using Carrier Modulation

clear all; close all; clc;
%  STEP 1: LOAD AUDIO FILE

fprintf('Step 1: Loading Audio File...\n');
% File selection dialog
[filename, pathname] = uigetfile({'*.wav', 'Wave Files (*.wav)'}, ...
                                 'Select Audio File to Scramble');
if filename == 0
    error('No file selected. Exiting...');
end

% Read audio file
[original_audio, Fs] = audioread(fullfile(pathname, filename));

% Convert stereo to mono if necessary
if size(original_audio, 2) > 1
    original_audio = mean(original_audio, 2);
    fprintf('  - Converted stereo to mono\n');
end

% Audio information
N = length(original_audio);
duration = N / Fs;
t = (0:N-1) / Fs;  % Time vector

fprintf('  ✓ Audio loaded successfully\n');
fprintf('    Sample Rate (Fs): %d Hz\n', Fs);
fprintf('    Duration: %.2f seconds\n', duration);
fprintf('    Number of samples: %d\n\n', N);

% Play original audio
fprintf('Playing ORIGINAL audio...\n');
sound(original_audio, Fs);
pause(duration + 1);
%  STEP 2: DESIGN LOW-PASS FILTER

fprintf('Step 2: Designing Low-Pass Filter...\n');

% Filter specifications
fc = 3000;           % Cutoff frequency (Hz) - below carrier frequency
filter_order = 8;    % Filter order (higher = sharper cutoff)

% Design Butterworth Low-Pass Filter
% Butterworth provides maximally flat passband response
% Normalize cutoff frequency: fc/(Fs/2) where Fs/2 is Nyquist frequency
[b, a] = butter(filter_order, fc/(Fs/2), 'low');

fprintf('  ✓ Butterworth Low-Pass Filter Designed\n');
fprintf('    Filter Order: %d\n', filter_order);
fprintf('    Cutoff Frequency: %d Hz\n', fc);
fprintf('    Normalized Cutoff: %.4f\n\n', fc/(Fs/2));

%  STEP 3: SCRAMBLING PROCESS (Frequency Inversion)

fprintf('Step 3: SCRAMBLING - Frequency Inversion Process...\n');

% Carrier frequency selection fc_carrier should be > (maximum voice frequency) so we use 3500-4000 Hz
fc_carrier = 3700;   % Carrier frequency in Hz

fprintf('  - Carrier Frequency: %d Hz\n', fc_carrier);

% Generate carrier signal: cos(2xpi*fc*t)
carrier = cos(2 * pi * fc_carrier * t');

% MODULATION STEP
% Mathematical Operation: y(t) = x(t) × cos(2π*fc*t)
% In Frequency Domain (using Fourier Transform properties):
% If X(f) is spectrum of x(t), then:
% Y(f) = 0.5 * [X(f - fc) + X(f + fc)]
% This creates two sidebands:
% - Upper Sideband (USB): fc to fc + fmax
% - Lower Sideband (LSB): fc - fmax to fc
% The LSB is the INVERTED spectrum of original signal

modulated_signal = original_audio .* carrier;

fprintf('  - Modulation complete\n');

% LOW-PASS FILTERING to remove the upper sideband and keeps only the lower sideband
% The lower sideband contains the frequency-inverted speech

scrambled_audio = filter(b, a, modulated_signal);

% Normalize to prevent clipping
scrambled_audio = scrambled_audio / max(abs(scrambled_audio)) * 0.95;

fprintf('Low-pass filtering complete\n');
fprintf('SCRAMBLING COMPLETE\n\n');

% Save scrambled audio
audiowrite('scrambled_audio.wav', scrambled_audio, Fs);
fprintf('Scrambled audio saved as: scrambled_audio.wav\n\n');

% Play scrambled audio
fprintf('Playing SCRAMBLED audio (unintelligible)...\n');
sound(scrambled_audio, Fs);
pause(duration + 1);

%  STEP 4: DESCRAMBLING PROCESS (Reverse Frequency Inversion)

fprintf('Step 4: DESCRAMBLING - Reversing Frequency Inversion...\n');

% IMPORTANT: Descrambling uses THE SAME PROCESS as scrambling!
% Mathematical Proof:
% Let s(t) be scrambled signal = x(t) × cos(ωc*t) [filtered]
% Descrambling: d(t) = s(t) × cos(ωc*t) [filtered]
%             = [x(t) × cos(ωc*t)] × cos(ωc*t)
%             = x(t) × cos²(ωc*t)
%             = x(t) × [1 + cos(2ωc*t)] / 2
% After low-pass filtering (removes 2ωc term):
%             d(t) = x(t) / 2  (Original signal recovered!)

% Generate same carrier
carrier_descramble = cos(2 * pi * fc_carrier * t');

% Modulate scrambled signal with same carrier
demodulated_signal = scrambled_audio .* carrier_descramble;

fprintf('Demodulation complete\n');

% Low-pass filter to remove high-frequency components
descrambled_audio = filter(b, a, demodulated_signal);

% Normalize
descrambled_audio = descrambled_audio / max(abs(descrambled_audio)) * 0.95;

fprintf('Low-pass filtering complete\n');
fprintf('DESCRAMBLING COMPLETE\n\n');

% Save descrambled audio
audiowrite('descrambled_audio.wav', descrambled_audio, Fs);
fprintf('Descrambled audio saved as: descrambled_audio.wav\n\n');

% Play descrambled audio
fprintf('Playing DESCRAMBLED audio (recovered original)...\n');
sound(descrambled_audio, Fs);
pause(duration + 1);
%  STEP 5: VISUALIZATION - TIME DOMAIN

fprintf('Step 5: Generating Visualizations...\n');

% Plot duration (show first 0.1 seconds for clarity)
plot_duration = min(0.1, duration);
plot_samples = round(plot_duration * Fs);

figure('Position', [50, 50, 1400, 900], 'Name', 'Time Domain Analysis');
sgtitle('Voice Scrambler - Time Domain Comparison', 'FontSize', 14, 'FontWeight', 'bold');

% Original Signal
subplot(3,1,1);
plot(t(1:plot_samples), original_audio(1:plot_samples), 'b', 'LineWidth', 1.5);
grid on;
xlabel('Time (seconds)');
ylabel('Amplitude');
title('(a) Original Audio Signal');
xlim([0 plot_duration]);

% Scrambled Signal
subplot(3,1,2);
plot(t(1:plot_samples), scrambled_audio(1:plot_samples), 'r', 'LineWidth', 1.5);
grid on;
xlabel('Time (seconds)');
ylabel('Amplitude');
title('(b) Scrambled Audio Signal (Frequency Inverted)');
xlim([0 plot_duration]);

% Descrambled Signal
subplot(3,1,3);
plot(t(1:plot_samples), descrambled_audio(1:plot_samples), 'g', 'LineWidth', 1.5);
grid on;
xlabel('Time (seconds)');
ylabel('Amplitude');
title('(c) Descrambled Audio Signal (Recovered)');
xlim([0 plot_duration]);
%  STEP 6: VISUALIZATION - FREQUENCY DOMAIN

fprintf('  - Computing FFT...\n');

% Compute FFT for all three signals and use zero-padding for better frequency resolution
NFFT = 2^nextpow2(N);
freq_axis = (0:NFFT-1) * (Fs/NFFT);

% Only plot positive frequencies (0 to Fs/2)
freq_mask = freq_axis <= Fs/2;

% FFT computation
fft_original = fft(original_audio, NFFT);
fft_scrambled = fft(scrambled_audio, NFFT);
fft_descrambled = fft(descrambled_audio, NFFT);

% Magnitude spectrum
mag_original = abs(fft_original) / N;
mag_scrambled = abs(fft_scrambled) / N;
mag_descrambled = abs(fft_descrambled) / N;

% Frequency Domain Visualization
figure('Position', [100, 100, 1400, 900], 'Name', 'Frequency Domain Analysis');
sgtitle('Voice Scrambler - Frequency Domain (Magnitude Spectrum)', ...
        'FontSize', 14, 'FontWeight', 'bold');

% Original Spectrum
subplot(3,1,1);
plot(freq_axis(freq_mask), mag_original(freq_mask), 'b', 'LineWidth', 1.5);
grid on;
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('(a) Original Audio - Frequency Spectrum');
xlim([0 5000]);

% Scrambled Spectrum (INVERTED!)
subplot(3,1,2);
plot(freq_axis(freq_mask), mag_scrambled(freq_mask), 'r', 'LineWidth', 1.5);
grid on;
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('(b) Scrambled Audio - Frequency Spectrum (Note: Inverted!)');
xlim([0 5000]);
% Add marker for carrier frequency
line([fc_carrier fc_carrier], ylim, 'Color', 'k', 'LineStyle', '--', ...
     'LineWidth', 1.5, 'DisplayName', 'Carrier Frequency');
legend('Location', 'best');

% Descrambled Spectrum (Recovered)
subplot(3,1,3);
plot(freq_axis(freq_mask), mag_descrambled(freq_mask), 'g', 'LineWidth', 1.5);
grid on;
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('(c) Descrambled Audio - Frequency Spectrum (Recovered)');
xlim([0 5000]);
%  STEP 7: SPECTROGRAMS (Time-Frequency Analysis)

figure('Position', [150, 150, 1400, 900], 'Name', 'Spectrogram Analysis');
sgtitle('Voice Scrambler - Spectrogram Comparison (Time-Frequency)', ...
        'FontSize', 14, 'FontWeight', 'bold');

% Original Spectrogram
subplot(3,1,1);
spectrogram(original_audio, hamming(512), 256, 512, Fs, 'yaxis');
title('(a) Original Audio - Spectrogram');
colorbar;
ylim([0 5]);

% Scrambled Spectrogram
subplot(3,1,2);
spectrogram(scrambled_audio, hamming(512), 256, 512, Fs, 'yaxis');
title('(b) Scrambled Audio - Spectrogram (Inverted Frequencies)');
colorbar;
ylim([0 5]);

% Descrambled Spectrogram
subplot(3,1,3);
spectrogram(descrambled_audio, hamming(512), 256, 512, Fs, 'yaxis');
title('(c) Descrambled Audio - Spectrogram (Recovered)');
colorbar;
ylim([0 5]);
%  STEP 8: QUALITY METRICS

fprintf('  - Computing quality metrics...\n');

% Calculate Mean Squared Error (MSE) between original and descrambled
mse_value = mean((original_audio - descrambled_audio).^2);

% Calculate correlation coefficient
correlation = corrcoef(original_audio, descrambled_audio);
correlation_value = correlation(1,2);

% Signal-to-Noise Ratio (SNR)
signal_power = mean(original_audio.^2);
noise_power = mean((original_audio - descrambled_audio).^2);
snr_db = 10 * log10(signal_power / noise_power);

fprintf('  ✓ Visualizations complete\n\n');
fprintf('SYSTEM PARAMETERS:\n');
fprintf('  • Sampling Frequency: %d Hz\n', Fs);
fprintf('  • Carrier Frequency: %d Hz\n', fc_carrier);
fprintf('  • LPF Cutoff Frequency: %d Hz\n', fc);
fprintf('  • Filter Order: %d\n', filter_order);
fprintf('  • Audio Duration: %.2f seconds\n\n', duration);

fprintf('FREQUENCY INVERSION PRINCIPLE:\n');
fprintf('  • Original Speech: 0 to ~3000 Hz\n');
fprintf('  • After Modulation: Creates sidebands at fc±f\n');
fprintf('  • Lower Sideband: fc-3000 to fc (INVERTED spectrum)\n');
fprintf('  • Upper Sideband: fc to fc+3000 (NORMAL spectrum)\n');
fprintf('  • LPF removes upper sideband, keeps inverted spectrum\n\n');

fprintf('QUALITY METRICS:\n');
fprintf('  • Mean Squared Error: %.6f\n', mse_value);
fprintf('  • Correlation Coefficient: %.4f\n', correlation_value);
fprintf('  • Signal-to-Noise Ratio: %.2f dB\n\n', snr_db);

fprintf('OUTPUT FILES:\n');
fprintf('  • scrambled_audio.wav (unintelligible)\n');
fprintf('  • descrambled_audio.wav (recovered)\n\n');

fprintf('   Applying frequency inversion TWICE returns the original!\n');
fprintf('   This is why scrambling and descrambling use same process.\n\n');

fprintf('Would you like to hear all three versions in sequence? (y/n): ');
response = input('', 's');

if strcmpi(response, 'y')
    fprintf('\nPlaying comparative sequence...\n\n');
    
    fprintf('1/3: ORIGINAL AUDIO\n');
    sound(original_audio, Fs);
    pause(duration + 1);
    
    fprintf('2/3: SCRAMBLED AUDIO (unintelligible)\n');
    sound(scrambled_audio, Fs);
    pause(duration + 1);
    
    fprintf('3/3: DESCRAMBLED AUDIO (recovered)\n');
    sound(descrambled_audio, Fs);
    pause(duration + 1);
end

fprintf('\n✓ Voice Scrambler Project Complete!\n');
fprintf('All visualizations and audio files have been generated.\n\n');