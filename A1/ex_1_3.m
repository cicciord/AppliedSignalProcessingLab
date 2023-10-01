close all
clear
clc

% variables definiton
A = 1;
T = 1;

fi = [10, 20, 100];
phii = 2 * pi * rand(1,3);

% computing signal oversampling
[t, x, f, X_psd] = compute_signal(1000, T, fi, phii, A);

% signal plot
figure
plot(t, x)
title('Finite domain signal');
xlabel('time [t]')
ylabel('x(t)')
grid on

% plot psd sampled @1000Hz
figure
plot(f, X_psd)
title("PSD of signal sampled at 1000Hz")
xlabel("frequency [Hz]")
ylabel("magnitude")
grid on
ylim([0 1.1*max(X_psd)])

% psd @160Hz
[~, ~, f, X_psd] = compute_signal(160, T, fi, phii, A);

figure
plot(f, X_psd)
title("PSD of signal sampled at 160Hz")
xlabel("frequency [Hz]")
ylabel("magnitude")
grid on
ylim([0 1.1*max(X_psd)])

% psd @120Hz
[~, ~, f, X_psd] = compute_signal(120, T, fi, phii, A);

% repeating undersampling with a new signal
phii_2 = 2 * pi * rand(1,3);
[~, ~, ~, X_psd_2] = compute_signal(120, T, fi, phii_2, A);

phii_3 = 2 * pi * rand(1,3);
[~, ~, ~, X_psd_3] = compute_signal(120, T, fi, phii_3, A);

figure
subplot(3, 1, 1)
plot(f, X_psd)
title("PSD of signal sampled at 120Hz")
grid on
ylim([0 1])

subplot(3, 1, 2)
plot(f, X_psd_2)
ylabel("magnitude")
grid on
ylim([0 1])

subplot(3, 1, 3)
plot(f, X_psd_3)
xlabel("frequency [Hz]")
grid on
ylim([0 1])

