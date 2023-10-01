close all
clear
clc

filename = "A_major_scale.wav";
y = audioread(filename);
info = audioinfo(filename);

Fs = info.SampleRate;
N = info.TotalSamples;

y = y / max(abs(y));

y_obj = audioplayer(y, Fs);
playblocking(y_obj);

% spectrogram
figure
spectrogram(y,hamming(10001),8000,2^14,Fs,'yaxis')
title("A major scale spectrogram")
ylim([0 4])

% highest note pitch
y_highest = y(165943:173486);
[R_highest, ~] = xcorr(y_highest, 'normalized');

[~, maxima_idx_highest] = findpeaks(R_highest(length(y_highest):end));
pitch_est_highest = Fs / (maxima_idx_highest(2) - maxima_idx_highest(1));
disp("The estimated pitch is: " + num2str(pitch_est_highest) + "Hz")

% my tuner

% check if the tuner works with known signal
[note, dist] = my_tuner(y_highest, Fs, length(y_highest));
disp("The highest note on the A Major Scale is: " + note + " - difference: " + num2str(dist))

% try tuner on 3 notes on a guitar
note_1_file = "note_1.wav";
note_2_file = "note_2.wav";
note_3_file = "note_3.wav";
y_1 = audioread(note_1_file);
y_2 = audioread(note_2_file);
y_3 = audioread(note_3_file);
y_1 = y_1/max(abs(y_1));
y_2 = y_2/max(abs(y_2));
y_3 = y_3/max(abs(y_3));

[note_1, dist_1] = my_tuner(y_1, Fs, length(y_1));
[note_2, dist_2] = my_tuner(y_2, Fs, length(y_2));
[note_3, dist_3] = my_tuner(y_3, Fs, length(y_3));

disp("The first note recordered is: " + note_1 + " - difference: " + num2str(dist_1))
disp("The second note recordered is: " + note_2 + " - difference: " + num2str(dist_2))
disp("The third note recordered is: " + note_3 + " - difference: " + num2str(dist_3))


% low pass filtering
Wc = 5000;
Wc_rad = Wc * 2 * pi / Fs;
B = 500;
B_rad = B * 2 * pi / Fs;

N_filt = ceil((6.6 * pi) / B_rad);

n_filt = 0 : N_filt - 1;
window = 0.54 - 0.46 * cos(2 * pi * n_filt / (N_filt-1));

M = ceil((N_filt-1)/2);

h_id = 2 * (Wc / Fs) * sinc(2 * (Wc / Fs) * (n_filt - M));

h = window .* h_id;

y_filt = filtfilt(h, 1, y);


% under-sampling
Fs_under = Fs/4;
y_under = y_filt(1:4:end);

figure
spectrogram(y_under,hamming(10001),8000,2^14,Fs_under,'yaxis')
title("A major scale spectrogram filtered (\Omega_c=5000)", Interpreter="tex")
ylim([0 4])

y_obj_under = audioplayer(y_under, Fs_under);
playblocking(y_obj_under);

% filtering with half cuttoff frequency
Wc_2 = 2500;
Wc_rad_2 = Wc_2 * 2 * pi / Fs;

h_id_2 = 2 * (Wc_2 / Fs) * sinc(2 * (Wc_2 / Fs) * (n_filt - M));

h_2 = window .* h_id_2;

y_filt_2 = filtfilt(h_2, 1, y);

Fs_under_2 = Fs/8;
y_under_2 = y_filt_2(1:8:end);

figure
spectrogram(y_under_2,hamming(10001),8000,2^14,Fs_under_2,'yaxis')
title("A major scale spectrogram filtered (\Omega_c=2500)", Interpreter="tex")

y_obj_2 = audioplayer(y_under_2, Fs_under_2);
playblocking(y_obj_2);



