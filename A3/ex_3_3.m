close all
clear
clc

filename = "A_major_scale.wav";
[x, Fs] = audioread(filename);

x = x' / max(abs(x));

M = 10001;
L = 8000;
N = 16384;
plot_flag = 0;

[S, f, t] = my_spectrogram(x, M, L, N, Fs, 0);

figure
imagesc(t, f, 20*log10(abs(S)))
title("Spectrogram custom function")
xlabel("time [s]")
ylabel("frequency [Hz]")
axis xy
ylim([0 4000])

[S0, f0, t0] = spectrogram(x, hamming(M), L, N, Fs, 'yaxis');

figure
imagesc(t0, f0, 20*log10(abs(S0)))
title("Spectrogram Matlab function")
xlabel("time [s]")
ylabel("frequency [Hz]")
axis xy
ylim([0 4000])


