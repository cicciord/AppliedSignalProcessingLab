close all
clear
clc

As1 = 40;
B1 = 400;
Fs1 = 8000;
Fc1 = 3200;
type1 = '-hp';

As2 = 60;
B2 = 200;
Fs2 = 8000;
Fc2 = [1600 2400];
type2 = '-bp';

[h1, N1, beta1] = my_Kaiser_filter(As1, B1, Fs1, Fc1, type1);
[h2, N2, beta2] = my_Kaiser_filter(As2, B2, Fs2, Fc2, type2);

b1 = fir1(N1-1, Fc1/(Fs1/2), kaiser(N1,beta1), 'high');
b2 = fir1(N2-1, Fc2/(Fs2/2), kaiser(N2,beta2), 'bandpass');

[H_cust1, f_cust1] = freqz(h1, 1, 1024, Fs1);
[H1, f1] = freqz(b1, 1, 1024, Fs1);

[H_cust2, f_cust2] = freqz(h2, 1, 1024, Fs2);
[H2, f2] = freqz(b2, 1, 1024, Fs2);

figure
plot(f_cust1, 20*log10(abs(H_cust1)), LineWidth=0.8)
title("High pass filter")
xlabel("frequency [Hz]")
ylabel("magnitude [dB]")
grid on
hold on
plot(f1, 20*log10(abs(H1)))
legend("Kaiser", "my\_Kaiser")

figure
plot(f_cust2, 20*log10(abs(H_cust2)), LineWidth=0.8)
title("Band pass filter")
xlabel("frequency [Hz]")
ylabel("magnitude [dB]")
grid on
hold on
plot(f2, 20*log10(abs(H2)))
legend("Kaiser", "my\_Kaiser")

