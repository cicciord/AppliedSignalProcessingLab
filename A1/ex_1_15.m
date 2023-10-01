close all
clear
clc

N = 1000;
M = 50;
x = randn(1, N);

[S, f] = my_Bartlett(x, N, M);

NFFT = length(S);
window = rectwin(NFFT);
n_overlap = 0;
[S_welch, f_welch] = pwelch(x, window, n_overlap, NFFT, 1, 'centered');

figure
plot(f_welch, S_welch, LineWidth=0.8)
title("Comparison between custom bartlett and bartlett from pwelch")
xlabel("frequency [Hz]")
ylabel("magnitude")
grid on
hold on
plot(f, S)
legend("pwelch", "my\_Bartlett")
