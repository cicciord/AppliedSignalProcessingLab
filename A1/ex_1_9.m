close all
clear
clc

% parameter definition
f0 = 20;
f1 = 16;
Fc = 64;
N = 128;
n = 0 : N-1;

% discrete signal definition
x_n = @(n) 3 * cos(2 * pi * f0 * n / Fc) + 2 * cos(2 * pi * f1 * n / Fc);

% find T0 and Df
T0 = N / Fc;
Df = Fc / N;

disp("The value of T0 is: " + num2str(T0))
disp("The value of Df is: " + num2str(Df))

% compute analog signal CTFT
syms t f
assume(t, 'real');
x_t = (3 * cos(2 * pi * f0 * t) + 2 * cos(2 * pi * f1 * t)) * rectangularPulse(0, T0, t);
X_f = fourier(x_t, f);

% compute a matlab function from the analytical formula
X = matlabFunction(X_f);
f = -Fc/2 : Df/100 : Fc/2 - Df/100;

% discrete signal DFT
X_n = fftshift(fft(x_n(n)));
f_n = -Fc/2 : Df : Fc/2 - Df;

figure
plot(f, abs(X(f)));
title("DFT of x(n) of " + num2str(N) + " number of samples")
xlabel("k , frequency [Hz]")
ylabel("magnitude")
hold on
stem(f_n, (T0 / N) * abs(X_n));
legend("CTFT", "DFT")

% change the value of N
N_diff = N + 2;
n_diff = 0 : N_diff-1;
Df_diff = Fc / N_diff;
f_n_diff = -Fc/2 : Df_diff : Fc/2 - Df_diff;
X_n_diff = fftshift(fft(x_n(n_diff)));

figure
plot(f, abs(X(f)));
title("DFT of x(n) of " + num2str(N_diff) + " number of samples")
xlabel("k , frequency [Hz]")
ylabel("magnitude")
hold on
stem(f_n_diff, (T0 / N_diff) * abs(X_n_diff));
legend("CTFT", "DFT")

% zero padding
N_pad = N + 128;
n_pad = 0 : N_pad-1;
Df_pad = Fc / N_pad;
T0_pad = N_pad * T0 / N;
f_n_pad = -Fc/2 : Df_pad : Fc/2 - Df_pad;
x_n_pad = [x_n(n) zeros(1, N_pad - N)];
X_n_pad = fftshift(fft(x_n_pad));

figure
plot(f, abs(X(f)));
title("DFT of x(n) zero-padded to " + num2str(N_pad) + " number of samples")
xlabel("k , frequency [Hz]")
ylabel("magnitude")
hold on
stem(f_n_pad, (T0_pad / N_pad) * abs(X_n_pad));
legend("CTFT", "DFT")

disp(' ')
disp("The new value of Df is: " + Df_pad)

