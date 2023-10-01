close all
clear
clc

filename = "A2_guitar.wav";
[data, Fs] = audioread(filename);

y = data(:,1);

T = length(y) / Fs;
Dt = 1 / Fs;
t = 0 : Dt : T - Dt;

figure
plot(t, y)
title("A2 guitar note waveform")
xlabel("time [s]")
ylabel("amplitude")
grid on

% spectral estimation
NFFT = 5000;
window = hamming(NFFT);
n_overlap = NFFT / 2;

[Y, f] = pwelch(y, window, n_overlap, NFFT, Fs, 'centered');

figure
plot(f, 20*log10(abs(Y)));
title("A2 guitar note spectrum")
xlabel("frequency [Hz]")
ylabel("magnitude [dB]")
grid on
xlim([-500 500])

% normalized autocorrelation
[R, lags] = xcorr(y, 'normalized');

[~, maxima_idx] = findpeaks(R(length(y):end));
pitch_est = 1 / (Dt * (maxima_idx(2) - maxima_idx(1)));
disp("The estimated pitch is: " + num2str(pitch_est) + "Hz")

% my voice
% voice has been registered with the script record_voice.m
filename_voice = "my_voice.wav";
[y_voice, Fs_voice] = audioread(filename_voice);

T_voice = length(y_voice) / Fs_voice;
Dt_voice = 1 / Fs_voice;
t_voice = 0 : Dt_voice : T_voice - Dt_voice;

figure
plot(t_voice, y_voice)
title("My voice waveform")
xlabel("time [s]")
ylabel("amplitude")
grid on

NFFT_voice = 8000;
window_voice = hamming(NFFT_voice);
n_overlap_voice = NFFT / 2;

[Y_voice, f_voice] = pwelch(y_voice, window_voice, n_overlap_voice, NFFT_voice, Fs_voice, 'centered');

figure
plot(f_voice, 20*log10(abs(Y_voice)));
title("My voice spectrum")
xlabel("frequency [Hz]")
ylabel("magnitude [dB]")
grid on
xlim([-500 500])

[R_voice, ~] = xcorr(y_voice, 'normalized');

[~, maxima_idx_voice] = findpeaks(R_voice(length(y_voice):end), MinPeakProminence=0.9);
pitch_est_voice = 1 / (Dt_voice * (maxima_idx_voice(2) - maxima_idx_voice(1)));
disp("The estimated pitch of the recorded voice is: " + num2str(pitch_est_voice) + "Hz")

% cut vocalized part
y_vocal = y_voice(8835:17079);
T_vocal = length(y_vocal) / Fs_voice;
t_vocal = 0 : Dt_voice : T_vocal - Dt_voice;

[R_vocal, ~] = xcorr(y_vocal, 'normalized');

[~, maxima_idx_vocal] = findpeaks(R_vocal(length(y_vocal):end), MinPeakProminence=0.9);
pitch_est_vocal = 1 / (Dt_voice * (maxima_idx_vocal(2) - maxima_idx_vocal(1)));
disp("The estimated pitch of the recorded voice is: " + num2str(pitch_est_vocal) + "Hz")



