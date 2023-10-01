close all
clear
clc

S = 70;

[b, a] = butter(8, 0.1);
X = [1 zeros(1, S-1)];

y = my_filter(b, a, X);
t = 0 : S-1;

figure
impz(b, a, S)
title("Impulse Response Comparison")
xlabel("time [s]")
ylabel("amplitude")
grid on
hold on
stem(t, y)
legend("impz", "my\_filter")

B = my_DFT(b, 512);
A = my_DFT(a, 512);

H_cust = B ./ A;
[H, W] = freqz(b, a, 256, 1);

figure
plot(W, 20*log10(abs(H)), LineWidth=0.8)
title("Frequency Response Comparison")
xlabel("frequency [Hz]")
ylabel("magnitude [dB]")
ylim([-100 10])
grid on
hold on
plot(W, 20*log10(abs(H_cust(length(H_cust)/2+1:end))))
legend("freqz", "my\_filter")

S = 200;
t = 0 : S-1;
x = randn(1, S);
[b, a] = ellip(6, 3, 60, [0.2 0.3]);
y = filter(b, a, x);
y_cust = my_filter(b, a, x);

figure
plot(t, y)
title("Random Process Filtering Comparison")
xlabel("time [s]")
ylabel("amplitude")
grid on
hold on
stem(t, y_cust)
legend("filter", "my\_filter")

b = [1 1.3 0.49 -0.013 -0.029];
a = [1 -0.4326 -1.6656 0.1253 0.2877];

try
    y = my_filter(b, a, X);
catch e
    disp(e.message)
end

