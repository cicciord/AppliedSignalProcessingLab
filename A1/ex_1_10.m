close all
clear
clc

% parameters definition
B1 = 8;
B2 = 3;
B3 = 2;
f0 = 6;
Fc = 32;
N = 128;

% axis
n = -N/2 : N/2 - 1;

% signal definition
x = 2 * B1 * sinc(2 * B1 * n / Fc) + B2 * sinc(B2 * n /Fc).^2 / 2 + 4 * B3 * sinc(B3 * n /Fc).^2 .* cos(2 * pi * f0 * n / Fc);

% compute T0 and Df
T0 = N / Fc;
Df = Fc / N;
disp("The value of T0 is: " + num2str(T0))
disp("The value of Df is: " + num2str(Df))

% define frequency axis for analog signal
M = 100 * N;
Df_analog = Fc / M;
f_analog = -Fc/2 : Df_analog : Fc/2 - Df_analog;
Y = rectangularPulse(-B1, B1, f_analog) + triangularPulse(-B2, 0, B2, f_analog) / 2 + 2 * triangularPulse(-B3 + f0, f0, B3 + f0, f_analog) + 2 * triangularPulse(-B3 - f0, -f0, B3 - f0, f_analog);

% DFT
f = -Fc/2 : Df : Fc/2 - Df;
X = fftshift(fft(x)) / Fc;

% plot comparison between CTFT and DFT
figure
plot(f_analog, Y, LineWidth=1.2)
title("Comparison between CTFT and DFT")
xlabel("k , frequency [Hz]")
ylabel("magnitude")
grid on
hold on
stem(f, abs(X))
legend("CTFT", "DFT")

% zero padding
N_pad = N + 128;
n_pad = -N_pad/2 : N_pad/2 - 1;
Df_pad = Fc / N_pad;
f_pad = -Fc/2 : Df_pad : Fc/2 - Df_pad;

disp(" ")
disp("The new value of Df after zero-padding is: " + num2str(Df_pad))

x_pad = 2 * B1 * sinc(2 * B1 * n_pad / Fc) + B2 * sinc(B2 * n_pad /Fc).^2 / 2 + 4 * B3 * sinc(B3 * n_pad /Fc).^2 .* cos(2 * pi * f0 * n_pad / Fc);
X_pad = fftshift(fft(x_pad)) / Fc;

figure
plot(f_analog, Y, LineWidth=1.2)
title("Comparison between CTFT and DFT")
xlabel("k , frequency [Hz]")
ylabel("magnitude")
grid on
hold on
stem(f_pad, abs(X_pad))
legend("CTFT", "DFT")

figure
plot(f_analog, Y, LineWidth=1.2)
title("Comparison between CTFT and DFT (zoomed)")
xlabel("k , frequency [Hz]")
ylabel("magnitude")
grid on
hold on
stem(f_pad, abs(X_pad))
legend("CTFT", "DFT")
axis([7.5 9 0 1.5])

% convolution between Z and W
Y_first = find(Y, 1, 'first');
Y_last = find(Y, 1, 'last');
Bw = (Y_last - Y_first) * Df_analog;
disp("The signal frequency bandwidth is: " + num2str(Bw))

Y_tilde = Y(Y_first : Y_last-1);

Q = M * (1 + Bw/Fc);
Df_bw = (Fc + Bw) / Q;
f_bw = -(Fc+Bw)/2 + (Df_bw)/2 : (Fc+Bw)/Q : (Fc+Bw)/2 - Df_bw/2;
W = T0 * sinc(T0 * f_bw);

Z = zeros(1, length(Y));
for i = 1 : length(Y)
    Z(i) = trapz(Df_analog, Y_tilde .* W(i:i+length(Y_tilde)-1));
end

figure
plot(f_analog, abs(Z), LineWidth=1.2)
title("Comparison between CTFT after convolution and DFT")
xlabel("k , frequency [Hz]")
ylabel("magnitude")
grid on
hold on
stem(f, abs(X))
legend("CTFT", "DFT")

figure
plot(f_analog, abs(Z), LineWidth=1.2)
title("Comparison between CTFT after convolution and DFT (zoomed)")
xlabel("k , frequency [Hz]")
ylabel("magnitude")
grid on
hold on
stem(f, abs(X))
legend("CTFT", "DFT")
axis([7.5 9 0 1.5])

