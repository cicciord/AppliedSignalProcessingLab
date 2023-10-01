close all
clear
clc

% parameters definition
N = 8;
Fc = 900e6;
c = 3e8;

lambda = c / Fc;

d = lambda / 2;
theta_1 = 0;

% compute and plot the array factor
[pattern, theta] = compute_pattern(d, lambda, theta_1, N);

figure
plot(theta, abs(pattern))
title("Array factor")
xlabel("theta [rad]")
ylabel("magnitude")
grid on

Psi = 2 * pi * d * (sin(theta) - sin(deg2rad(theta_1))) / lambda;
figure
plot(theta, abs(diric(Psi, N)))
title("Array factor (diric)")
xlabel("theta [rad]")
ylabel("magnitude")
grid on

figure
polarplot(theta, abs(pattern), 'r', 'LineWidth', 2);
ax = gca;
ax.ThetaZeroLocation = 'top';
ax.ThetaDir = 'clockwise';
ax.ThetaLim = [-90 90];
title(['Conventional beamforming with N=',num2str(N),', DoA \theta_1=', num2str(theta_1), '^o'], 'FontSize' ,14);

% change values for N
N_1 = 2;
N_2 = 12;
N_3 = 20;
N_4 = 50;
[pattern_N_1, ~] = compute_pattern(d, lambda, theta_1, N_1);
[pattern_N_2, ~] = compute_pattern(d, lambda, theta_1, N_2);
[pattern_N_3, ~] = compute_pattern(d, lambda, theta_1, N_3);
[pattern_N_4, ~] = compute_pattern(d, lambda, theta_1, N_4);

disp("The maximum of the array factor for N=" + num2str(N_1) + ...
    " is " + num2str(max(abs(pattern_N_1))) + " and the sqrt(N)=" + num2str(sqrt(N_1)))
disp("The maximum of the array factor for N=" + num2str(N_2) + ...
    " is " + num2str(max(abs(pattern_N_2))) + " and the sqrt(N)=" + num2str(sqrt(N_2)))
disp("The maximum of the array factor for N=" + num2str(N_3) + ...
    " is " + num2str(max(abs(pattern_N_3))) + " and the sqrt(N)=" + num2str(sqrt(N_3)))
disp("The maximum of the array factor for N=" + num2str(N_4) + ...
    " is " + num2str(max(abs(pattern_N_4))) + " and the sqrt(N)=" + num2str(sqrt(N_4)))

figure
subplot(2, 2, 1)
plot(theta, abs(pattern_N_1))
title("Array factor (N=" + num2str(N_1) + ")")
xlabel("theta [rad]")
ylabel("magnitude")
grid on
subplot(2, 2, 2)
plot(theta, abs(pattern_N_2))
title("Array factor (N=" + num2str(N_2) + ")")
xlabel("theta [rad]")
ylabel("magnitude")
grid on
subplot(2, 2, 3)
plot(theta, abs(pattern_N_3))
title("Array factor (N=" + num2str(N_3) + ")")
xlabel("theta [rad]")
ylabel("magnitude")
grid on
subplot(2, 2, 4)
plot(theta, abs(pattern_N_4))
title("Array factor (N=" + num2str(N_4) + ")")
xlabel("theta [rad]")
ylabel("magnitude")
grid on

% change values for d/lambda
d_1 = 0.25 * lambda;
d_2 = 0.5 * lambda;
d_3 = 1 * lambda;
d_4 = 2 * lambda;
[pattern_d_1, ~] = compute_pattern(d_1, lambda, theta_1, N);
[pattern_d_2, ~] = compute_pattern(d_2, lambda, theta_1, N);
[pattern_d_3, ~] = compute_pattern(d_3, lambda, theta_1, N);
[pattern_d_4, ~] = compute_pattern(d_4, lambda, theta_1, N);

figure
subplot(2, 2, 1)
plot(theta, abs(pattern_d_1))
title("Array factor (d/\lambda=" + num2str(d_1/lambda) + ")")
xlabel("theta [rad]")
ylabel("magnitude")
grid on
subplot(2, 2, 2)
plot(theta, abs(pattern_d_2))
title("Array factor (d/\lambda=" + num2str(d_2/lambda) + ")")
xlabel("theta [rad]")
ylabel("magnitude")
grid on
subplot(2, 2, 3)
plot(theta, abs(pattern_d_3))
title("Array factor (d/\lambda=" + num2str(d_3/lambda) + ")")
xlabel("theta [rad]")
ylabel("magnitude")
grid on
subplot(2, 2, 4)
plot(theta, abs(pattern_d_4))
title("Array factor (d/\lambda=" + num2str(d_4/lambda) + ")")
xlabel("theta [rad]")
ylabel("magnitude")
grid on

theta_max_1 = rad2deg(find_maxima(theta_1, lambda, d_1));
theta_max_2 = rad2deg(find_maxima(theta_1, lambda, d_2));
theta_max_3 = rad2deg(find_maxima(theta_1, lambda, d_3));
theta_max_4 = rad2deg(find_maxima(theta_1, lambda, d_4));
disp(" ")
disp("The maxima for d\lambda=" + num2str(d_1/lambda) + " are in rad:")
for i = 1 : length(theta_max_1)
    disp(num2str(theta_max_1(i)))
end
disp("The maxima for d\lambda=" + num2str(d_2/lambda) + " are in rad:")
for i = 1 : length(theta_max_2)
    disp(num2str(theta_max_2(i)))
end
disp("The maxima for d\lambda=" + num2str(d_3/lambda) + " are in rad:")
for i = 1 : length(theta_max_3)
    disp(num2str(theta_max_3(i)))
end
disp("The maxima for d\lambda=" + num2str(d_4/lambda) + " are in rad:")
for i = 1 : length(theta_max_4)
    disp(num2str(theta_max_4(i)))
end


% nulls
x = sym("x");
last_null = asin(sin(deg2rad(theta_1)) + 0.5 * x);
sol = solve(last_null == pi/2, x);
sol = 1/eval(sol);

disp(" ")
disp("The value of d/lambda which gives the best directivity w/o grating lobes is: " + num2str(sol))


% FNBW
theta_1 = 30;
first_null = asin(sin(deg2rad(theta_1)) + 1/N * lambda/d);
first_neg_null = asin(sin(deg2rad(theta_1)) - 1/N * lambda/d);
FNBW = abs(first_null-first_neg_null);
disp(" ")
disp("The FNBW is: " + num2str(rad2deg(FNBW)) + char(176))

N_min = N;
while isreal(asin(sin(deg2rad(theta_1)) + 1/N_min * lambda/d))
    N_min = N_min - 1;
end
N_min = N_min + 1;
disp("The minimum N to have theta_null_1 less than pi/2 is: " + num2str(N_min))

% INTERNAL FUNCTIONS
function maxima = find_maxima(theta_deg, lambda, d)
    maxima = theta_deg;
    m = 1;
    while 1
        maxima(2*m) = asin(sin(deg2rad(theta_deg)) + m * lambda/d);
        maxima(2*m+1) = asin(sin(deg2rad(theta_deg)) - m * lambda/d);
        if ~isreal(maxima(2*m+1))
            break
        end
        m = m + 1;
    end
    maxima = maxima(1 : end - 2);
end


function [pattern, theta] = compute_pattern(d, lambda, theta_deg, N)
    S = 3601;
    theta = linspace(-pi/2, pi/2, S);
    
    Psi = 2 * pi * d * sin(deg2rad(theta_deg)) / lambda;

    m = (0 : N-1)';
    V = exp(2j * pi * d * m * sin(theta) / lambda);
    
    w = exp(1j * m * Psi) / sqrt(N);

    pattern = w' * V;
end


