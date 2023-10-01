close all
clear
clc

% constants definition
T = 2;
Fc = 10e3;
mu = 0;
sigma_sq = 25;

sigma = sqrt(sigma_sq);

N = T * Fc;

Dt = 1 / Fc;
t = 0 : Dt : T - Dt;

% generate WSN signal
w = mu + sigma * randn(1, N);

% plot first 100ms of the sample
T1 = 100e-3;
N1 = T1 * Fc;

figure
plot(t(1:N1), w(1:N1));
title("WGN Signal")
xlabel("time [s]")
ylabel("w(t)")
grid on

% compute mean, standard deviation and variance
mu_w = mean(w);
sig_w = std(w);
sig_sq_w = var(w);

disp("W(n) mean value is:            " + num2str(mu_w))
disp("W(n) standard deviation is:    " + num2str(sig_w))
disp("W(n) variance is:              " + num2str(sig_sq_w))

% create vector of edges for histogram
Nbins = 120;
Dedges = 8 * sigma / (Nbins+1);
% to create 120 bins there are needed 121 edges non symmetrical for the 0
% not to be an edge
edges = -4*sigma : Dedges : 4*sigma - Dedges;
% difine points at ceter of each bin
binc = -4*sigma + Dedges/2 : Dedges : 4*sigma - 3*Dedges/2;

% plot histogram
figure
histogram(w, edges, Normalization="pdf")
title("WSN Signal PDF")
grid on
hold on
plot(binc, normpdf(binc, mu, sigma), LineWidth=1.2)
legend("experimental", "theoretical")

% compute autocorrelation
Rw = xcorr(w);

% plot autocorrelation
T_R = 10e-3;
N_R = T_R * Fc;

t_R = -T_R : Dt : T_R - Dt;
figure
plot(t_R, Rw(floor(length(Rw)/2) - N_R+1 : floor(length(Rw)/2) + N_R))
title("Auto-correlation of w")
xlabel("time [t]")
ylabel("R_w(t)")
grid on

% filter
N_filter = 80;
h = firpm(N_filter, [0 0.05 0.1 1], [1 1 0 0]);
z = filter(h,1,w);

% plot comparison with filtered signal
figure
plot(t(N_filter/2+1:N1), w(N_filter/2+1:N1));
title("Comparison of w(t) and z(t)")
xlabel("time [s]")
ylabel("w(t),z(t)")
grid on
hold on
plot(t(N_filter/2+1:N1), z(N_filter/2+1:N1), LineWidth=2);
legend("w(t)", "z(t)")

% plot pdf after filtering
sigma_sq_z = sigma_sq * sum(abs(h).^2);
sigma_z = sqrt(sigma_sq_z);
Dedges_z = 8 * sigma_z / (Nbins+1);
edges_z = -4*sigma_z : Dedges_z : 4*sigma_z - Dedges_z;
binc_z = -4*sigma_z + Dedges_z/2 : Dedges_z : 4*sigma_z - 3*Dedges_z/2;

figure
histogram(z, edges_z, Normalization="pdf")
title("WSN Signal PDF after filtering")
grid on
hold on
plot(binc_z, normpdf(binc_z, mu, sigma_z), LineWidth=1.2)
legend("experimental", "theoretical")

% auto-correlation comparison
Rz = xcorr(z);

figure
plot(t_R, Rw(floor(length(Rw)/2) - N_R+1 : floor(length(Rw)/2) + N_R))
title("Comparison between Rw and Rz")
xlabel("time [t]")
ylabel("R_w(t),R_z(t)")
grid on
hold on
plot(t_R, Rz(floor(length(Rz)/2) - N_R+1 : floor(length(Rz)/2) + N_R), LineWidth=1.2)
legend("R_w(t)", "R_z(t)")

