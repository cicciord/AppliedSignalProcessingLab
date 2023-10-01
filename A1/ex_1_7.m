close all
clear
clc

% triangulat sequence length
N = 25;

% rectangular sequence length
M = (N + 1) / 2;

% rectangular pulse
r = ones(1, M);

% triangular signal
n = 0 : N-1;
x = my_conv(r, r);
x = x / M;

% plot triangular signal
figure
stem(n, x)
title("Triangular signal")
xlabel("sample number")
ylabel("value")
axis([0 N-1 0 1.1])
grid on

% compute DFT
f = -1/2 + 1/(2*N) : 1/N : 1/2 - 1/(2*N);
X = fftshift(fft(x));

% plot DFT
figure
stem(f, abs(X))
title("Triangular signal DFT magnitude")
xlabel("frequency [Hz]")
ylabel("magnitude")
grid on

% DTFT
N_DTFT = 1000;
f_DTFT = -1/2 : 1/N_DTFT : 1/2 - 1/N_DTFT;
X_DTFT = (sin(pi * f_DTFT * M).^2) ./ (M * sin(pi * f_DTFT).^2);

% zero padding
N_padd = [32, 64, 128];
for N_i = N_padd
    % add zeros to signal
    x_i = [x zeros(1, N_i - length(x))];
    f_i = -1/2 : 1/N_i : 1/2 - 1/N_i;
    X_i = fftshift(fft(x_i));
    
    % plot DFT and DTFT
    figure
    plot(f_DTFT, X_DTFT)
    title("Effect of zero padding (N=" + num2str(N_i) + ") and comparison with DTFT")
    xlabel("frequency [Hz]")
    ylabel("magniutude")
    grid on
    hold on
    stem(f_i, abs(X_i))
    legend("DTFT", "DFT")
end
