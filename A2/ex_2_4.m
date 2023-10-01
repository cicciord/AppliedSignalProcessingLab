close all
clear
clc

Fs = 20;
B = 8;
T = 160;
M = 8000;

Df = Fs / M;
f = -Fs/2 : Df : Fs/2 - Df;

Dt = 1/Fs;
t = -T/2 : Dt : T/2 - Dt;

Si_l = (1/pi) * sinint(pi * T * (f + B/2));
Si_r = -(1/pi) * sinint(pi * T * (f - B/2));
Xt = Si_l + Si_r;

figure
plot(f, Si_l, "r")
title("Gibbs Phenomenon, B=" + num2str(B) + ", T=" + num2str(T))
xlabel("freqency [Hz]")
ylabel("Sine integral")
grid on
hold on
plot(f, Si_r, "g")
hold on
plot(f, Xt, "k")
legend("Si(\piT(f+B/2))/pi", "-Si(\piT(f-B/2))/pi", "[Si(\piT(f+B/2))-Si(\piT(f-B/2))]/pi")

N = T/Dt;
disp("The required number of samples is: " + num2str(N));

x = B * sinc(B*t);

figure
plot(t, x)
title("Sinc, B=" + num2str(B))
xlabel("time [s]")
ylabel("amplitude")
grid on

N_pad = 2^15;
L = N_pad - length(x);
x = [x zeros(1, L)];

X = fftshift(fft(x)) / Fs;
Df_pad = Fs/N_pad;
f_pad = -Fs/2 : Df_pad : Fs/2 - Df_pad;

figure
plot(f_pad, abs(X), LineWidth=0.8)
title("Comparison between DTF and CTFT")
xlabel("frequency [Hz]")
ylabel("magnitude")
grid on
hold on
plot(f, abs(Xt))
legend("DFT", "CTFT")

