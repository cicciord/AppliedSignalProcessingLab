close all
clear
clc

% frequency axis
f = linspace(-1, 1, 1001);

% DTFT of sinc function
X = @(f, L) (sin(pi * f * L)) .* (exp(1j * pi * f * (L-1))) ./ (sin(pi * f));

% DTFT plots
DTFT_plot(f, X, 3, "Sinc DTFT magnitude");
DTFT_plot(f, X, 7, "Sinc DTFT magnitude");

Z = @(f, L) (sin(pi * f * L)) ./ (sin(pi * f));
Z_plot(f, Z, 3, "Periodic Sinc");
Z_plot(f, Z, 7, "Periodic Sinc");

L = 13;
W = L * sinc(f * L);

figure;
plot(f, W, LineWidth=0.8)
title("Z and W comparison (L = 13)");
xlabel("frequency [Hz]")
ylabel("Z, W")
grid on
hold on
plot(f, Z(f, L))
legend("W", "Z")
ylim([-5 15])

% MSE of Z and W
threshold = 0.16^2;
E = abs(Z(f, 13) - W).^2;
figure
plot(f, E)
title("MSE of Z and W")
xlabel("frequancy [Hz]")
ylabel("MSE")
grid on
xlim([0 0.5])
yline(threshold, "red")
legend("MSE", num2str(threshold))

% count the number of lobes
E = fftshift(E); % shift axis to only consider positive side
idx_threshold = find(E > threshold, 1); % index at which E is higher than threshold

n_lobes = sum(islocalmin(E(1:idx_threshold)));
disp("E(f) has " + num2str(n_lobes) + " lobes present in the range where E(f) < " + num2str(threshold) + " (onesided)")

