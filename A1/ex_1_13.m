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
T_D = round((T_tot - T_C) * rand(1), 2);

% generate delayed signal
N_D = round(T_D * Fc); % round used for safety (error for value of the order of 1e3)

P = 20;
F0 = randi(P) * 50;

p = c .* cos(2 * pi * F0 * t(1:length(c)));

s = [zeros(1, N_D) p zeros(1, N - N_D - M)];

% generate w and x
w = mu + sigma * randn(1, N);
x = s + w;

% compute P cross-correlation
A = zeros(P, N + M - 1);
for i = 1 : P
    F0_i = i * 50;
    p_i = c .* cos(2 * pi * F0_i * t(1:length(c)));
    A(i,:) = my_xcorr(x, p_i);
end

% plot cross-correlation mesh
lags = -M+1 : N-1;
figure
mesh(lags, 50:50:P*50, A)
title("Cross-correlation matrix for all F0 values")
xlabel("lags")
ylabel("frequency [Hz]")
grid on

% find max position
[~, lin_idx] = max(A, [], 'all');

% get x and y indices from linear index
f_hat = mod(lin_idx - 1, P) + 1;
d_hat = floor((lin_idx - 1) / P) + 1;

T_D_hat = Tc * lags(d_hat);
F0_hat = f_hat * 50;
disp("The extimated delay is                :    " + num2str(T_D_hat))
disp("The delay set was                     :    " + num2str(T_D))
disp("The extimated modulation frequency is :    " + num2str(F0_hat))
disp("The modulation frequency set was      :    " + num2str(F0))
