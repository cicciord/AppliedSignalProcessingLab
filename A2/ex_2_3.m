close all
clear
clc

T = 20;
Fc = 500;
f1 = 5;
f2 = 100;
mu = 0;
sigma_sq = 0.1;

sigma = sqrt(sigma_sq);
N = T*Fc;

Dt = 1/Fc;
t = 0 : Dt : T - Dt;

w = mu + sigma * randn(1, N);

x = cos(2 * pi * f1 * t) + cos(2 * pi * f2 * t) + w;

NFFT = length(x)/25;
window = hamming(NFFT);
n_overlap = floor(NFFT/2);
[Sx, f] = pwelch(x, window, n_overlap, NFFT, Fc, 'centered');

figure
plot(f, 10*log10(abs(Sx)))
title("Estimated PSD")
xlabel("frequency [Hz]")
ylabel("magnitude [dB]")
grid on


Wp = [90 110];
Wp = Wp/(Fc/2);
Ws = [80 120];
Ws = Ws/(Fc/2);
Rp = 1;
Rs = 60;

[N, Wn] = buttord(Wp, Ws, Rp, Rs);
[B, A] = butter(N, Wn, "bandpass");

figure
freqz(B, A, 1024);
ylim([-100 10])

y1 = filter(B, A, x);
Sy1 = pwelch(y1, window, n_overlap, NFFT, Fc, 'centered');

figure
plot(f, 10*log10(abs(Sy1)))
title("PSD after filtering")
xlabel("frequency [Hz]")
ylabel("magnitude [dB]")
grid on
ylim([-20 -6])
xlim([-Fc/2 Fc/2])

[N2, Wp2] = ellipord(Wp, Ws, Rp, Rs);
[B2, A2] = ellip(N2, Rp, Rs, Wp2, "bandpass");
y2 = filter(B2, A2, x);
Sy2 = pwelch(y2, window, n_overlap, NFFT, Fc, 'centered');

figure
freqz(B2, A2, 1024);
ylim([-100 10])

figure
plot(f, 10*log10(abs(Sy2)))
title("PSD after filtering")
xlabel("frequency [Hz]")
ylabel("magnitude [dB]")
grid on
ylim([-20 -6])
xlim([-Fc/2 Fc/2])

[N3, F03, A03, W3] = firpmord([Ws(1) Wp Ws(2)], [0 1 0], [10^(-Rs/20) 1-10^(-Rp/20) 10^(-Rs/20)]);
B3 = firpm(N3, F03, A03, W3, "bandpass");
y3 = filter(B3, 1, x);
Sy3 = pwelch(y3, window, n_overlap, NFFT, Fc, 'centered');

figure
freqz(B3, 1, 1024);
ylim([-100 10])

figure
plot(f, 10*log10(abs(Sy3)))
title("PSD after filtering")
xlabel("frequency [Hz]")
ylabel("magnitude [dB]")
grid on
ylim([-20 -6])
xlim([-Fc/2 Fc/2])

Wp_cheb = 6;
Wp_cheb = Wp_cheb/(Fc/2);
Ws_cheb = 10;
Ws_cheb = Ws_cheb/(Fc/2);

[N4, Wp] = cheb1ord(Wp_cheb, Ws_cheb, Rp, Rs);
[B4, A4] = cheby1(N4, Rp, Wp);
y4 = filter(B4, A4, x);
Sy4 = pwelch(y4, window, n_overlap, NFFT, Fc, 'centered');

figure
freqz(B4, A4, 1024);
ylim([-100 10])

figure
plot(f, 10*log10(abs(Sy4)))
title("PSD after filtering")
xlabel("frequency [Hz]")
ylabel("magnitude [dB]")
grid on
ylim([-20 -6])
xlim([-Fc/2 Fc/2])


figure
plot(t(1:4*Fc), x(1:4*Fc))
title("Comparison of the signal after low pass filtering")
xlabel("time [s]")
ylabel("amplitude")
grid on
hold on
plot(t(1:4*Fc), y4(1:4*Fc), LineWidth=2)
legend("original signal", "filtered signal")





