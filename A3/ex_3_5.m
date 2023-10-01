close all
clear
clc

Fs = 11025;
Fc = 1760;
Fm = 4.4;
I0 = 120;
A = 1;

T = 5;
Dt = 1/Fs;
t = 0 : Dt : T - Dt;

y_1 = A * cos(2 * pi * Fc * t + I0 * sin(2 * pi * Fm * t));

M = 256;
L = 128;
N = 512;
figure
spectrogram(y_1,hamming(M),L,N,Fs,'yaxis')
title("Alarm spectrogram")

Fc_2 = 200;
Fm_2 = 280;
I0_2 = 10;

T_2 = 10;
t_2 = 0 : Dt : T_2 - Dt;

a = exp(-t_2/2);

y_2 = A * a .* cos(2 * pi * Fc_2 * t_2 + I0_2 * a .* sin(2 * pi * Fm_2 * t_2));

figure
spectrogram(y_2,hamming(M),L,N,Fs,'yaxis')
title("Bell spectrogram")

y_1_obj = audioplayer(y_1, Fs);
y_2_obj = audioplayer(y_2, Fs);
playblocking(y_1_obj, Fs)
playblocking(y_2_obj, Fs)
audiowrite("alarm.wav", y_1, Fs)
audiowrite("bell.wav", y_2, Fs)


