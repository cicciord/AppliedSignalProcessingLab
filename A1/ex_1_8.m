close all
clear
clc

% define symbolic variables
syms t f_n B;
sympref('FourierParameters', [1 -2*pi]);
assume(t, 'real');
assume(in(B, 'real') & B > 0);

% define time unlimited signal
x = 36 * t * exp(-12 * t) * heaviside(t);

% energy in time domain
E_t = int(abs(x)^2, t, 0, inf);
disp("Energy of the signal:")
disp(E_t)

% fourier transform of the signal
X = fourier(x, f_n);
disp("CTFT of the signal:")
disp(X)

% energy in frequency domain as function of the bandwidth
E_f = int(abs(X)^2, f_n, -B, B);
disp("Bandwidth integral:")
disp(E_f)

% find the bandwidth s.t. energy is 0.999 of total energy
B999 = vpasolve(E_f == 0.999 * E_t, B);
disp("The 99.9% energy bandwidth is:")
disp(num2str(eval(B999)))

Fc = 8 * eval(B999);
T0 = max(abs(eval(solve(x < 1e-3, t))));
N = T0 * Fc;

% find N as a power of 2
N = 2^nextpow2(N);

% recompute T0 for the new N
T0 = N / Fc;


% discrete time signal
Tc = 1 / Fc;
t_n = 0 : Tc : T0 - Tc;
x_n = 36 * t_n .* exp(-12 * t_n);

% evaluate conntinuous time signal
M = 100 * N;
Tc_t = T0 / M;
t_t = 0 : Tc_t : T0 - Tc_t;
x_t = matlabFunction(x);

% plot time domain signals
figure
plot(t_t, x_t(t_t), LineWidth=1.4)
title("Comparison between continuous-time and discrete-time signals")
xlabel("time [s]")
ylabel("x(t) , x(n)")
grid on
hold on
stem(t_n, x_n)

% DFT of discrete signal
Df_n = Fc / N;
f_n = -Fc/2 : Df_n : Fc/2 - Df_n;
X_n = fftshift(fft(x_n));
disp(" ")
disp("Df of discrete time is:")
disp(num2str(Df_n))

% evaluate Fourier Transform of continuous time signal
Df_t = Fc / M;
f_t = -Fc/2 : Df_t : Fc/2 - Df_t;
X_t = matlabFunction(X);

% plot frequency domain signals
figure
plot(f_t, abs(X_t(f_t)), LineWidth=1.4)
title("Comparison between CTFT and DFT")
xlabel("frequency [Hz]")
ylabel("X(f) , X(k)")
grid on
hold on
stem(f_n, Tc*abs(X_n))

figure
plot(f_t, abs(X_t(f_t)), LineWidth=1.4)
title("Comparison between CTFT and DFT (zoomed)")
xlabel("frequency [Hz]")
ylabel("X(f) , X(k)")
grid on
hold on
stem(f_n, Tc*abs(X_n))
xlim([Fc*0.25/2 Fc*1.025/2])


