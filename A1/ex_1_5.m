close all
clear
clc

% number of samples
N = 256;
% random vector on N samples
x = randn(1, N);

% compute DFT
X = my_DFT(x, N);
X_fft = fftshift(fft(x));

% frequency axis
f = -N/2 : N/2 - 1;

% comparison plot
figure;
subplot(2, 1, 1)
plot(f, abs(X), LineWidth=0.8)
title("Magnitude")
ylabel("magnitude")
grid on
hold on
plot(f, abs(X_fft))
legend("my\_DFT", "fft")
subplot(2, 1, 2)
plot(f, angle(X), LineWidth=0.8)
title("Phase")
xlabel("frequency [Hz]")
ylabel("phase")
grid on
hold on
plot(f, angle(X_fft))
legend("my\_DFT", "fft")
