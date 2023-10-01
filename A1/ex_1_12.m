close all
clear
clc

% parameters definition
T_tot = 10;
T_C = 4;
T_B = 2e-3;
Fc = 2e3;
mu = 0;
sigma_sq = 25;

sigma = sqrt(sigma_sq);

N = T_tot * Fc;

Tc = 1 / Fc;
B = T_C / T_B;

t = 0 : Tc : T_tot - Tc;

% generate Bernoulli process
b = binornd(1, 0.5, 1, B);

% map 0 to -1
b = -(b == 0) + (b == 1);

% generate c from b
N_B = T_B * Fc;
c = repmat(b, N_B, 1);
c = c(:)';

% number of samples of c
M = length(c);

% generate T_D
T_D = round((T_tot - T_C) * rand(1, 1), 2);

% generate delayed signal
N_D = T_D * Fc;
s = [zeros(1, N_D) c zeros(1, N - N_D - M)];

% generate w and x
w = mu + sigma * randn(1, N);
x = s + w;

% plot signal with and without noise
figure
subplot(2, 1, 1)
plot(t, s)
title("Useful signal s(t)")
xlabel("time [t]")
ylabel("s(t)")
grid on
subplot(2, 1, 2)
plot(t, x)
title("Signal with noise")
xlabel("time [t]")
ylabel("x(t)")
grid on

% compute cross-correlation
r = my_xcorr(x, c);
lags = -M+1 : N-1;

% plot auto-correlation
figure
plot(lags, r)
title("Cross-correlation between x(t) and c(t)")
xlabel("lags")
ylabel("r(t)")
grid on

% find maximum value
[~, d_hat] = max(r);
disp("The position of the max in the cross-correlation sequence d_hat is: " + num2str(d_hat))

T_D_hat = Tc * lags(d_hat);
disp("The extimated delay is :    " + num2str(T_D_hat))
disp("The delay set was      :    " + num2str(T_D))
