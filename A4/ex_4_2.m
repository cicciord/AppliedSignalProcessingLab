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

theta_int = [20, -40, 60, -75, 80];

[pattern, theta] = compute_pattern(d, lambda, theta_1, N, theta_int);

figure
polarplot(theta, abs(pattern), 'r', 'LineWidth', 2);
ax = gca;
ax.ThetaZeroLocation = 'top';
ax.ThetaDir = 'clockwise';
ax.ThetaLim = [-90 90];
[~, min_int] = min(abs(theta_1 - theta_int));
min_int = theta_int(min_int);
title(['Capon beamforming with N=',num2str(N),', DoA \theta_1=', num2str(theta_1), '^o', ', DoA closest interferer = ', num2str(min_int), '^o'], 'FontSize' ,14);
hold on
polarplot([deg2rad(theta_int); deg2rad(theta_int)], [zeros(1, length(theta_int)); ones(1, length(theta_int))]*max(abs(pattern)), '--k');


theta_int(1) = 10;

[pattern, theta] = compute_pattern(d, lambda, theta_1, N, theta_int);

figure
polarplot(theta, abs(pattern), 'r', 'LineWidth', 2);
ax = gca;
ax.ThetaZeroLocation = 'top';
ax.ThetaDir = 'clockwise';
ax.ThetaLim = [-90 90];
[~, min_int] = min(abs(theta_1 - theta_int));
min_int = theta_int(min_int);
title(['Capon beamforming with N=',num2str(N),', DoA \theta_1=', num2str(theta_1), '^o', ', DoA closest interferer = ', num2str(min_int), '^o'], 'FontSize' ,14);
hold on
polarplot([deg2rad(theta_int); deg2rad(theta_int)], [zeros(1, length(theta_int)); ones(1, length(theta_int))]*max(abs(pattern)), '--k');




theta_int(1) = 5;

[pattern, theta] = compute_pattern(d, lambda, theta_1, N, theta_int);

figure
polarplot(theta, abs(pattern), 'r', 'LineWidth', 2);
ax = gca;
ax.ThetaZeroLocation = 'top';
ax.ThetaDir = 'clockwise';
ax.ThetaLim = [-90 90];
[~, min_int] = min(abs(theta_1 - theta_int));
min_int = theta_int(min_int);
title(['Capon beamforming with N=',num2str(N),', DoA \theta_1=', num2str(theta_1), '^o', ', DoA closest interferer = ', num2str(min_int), '^o'], 'FontSize' ,14);
hold on
polarplot([deg2rad(theta_int); deg2rad(theta_int)], [zeros(1, length(theta_int)); ones(1, length(theta_int))]*max(abs(pattern)), '--k');


theta_int = linspace(-90, 90, 31);

[pattern, theta] = compute_pattern(d, lambda, theta_1, N, theta_int);

figure
polarplot(theta, abs(pattern), 'r', 'LineWidth', 2);
ax = gca;
ax.ThetaZeroLocation = 'top';
ax.ThetaDir = 'clockwise';
ax.ThetaLim = [-90 90];
[~, min_int] = min(abs(theta_1 - theta_int));
min_int = theta_int(min_int);
title(['Capon beamforming with N=',num2str(N),', DoA \theta_1=', num2str(theta_1), '^o', ', DoA closest interferer = ', num2str(min_int), '^o'], 'FontSize' ,14);
hold on
polarplot([deg2rad(theta_int); deg2rad(theta_int)], [zeros(1, length(theta_int)); ones(1, length(theta_int))]*max(abs(pattern)), '--', 'Color', [0, 0, 0, 0.2]);




% INTERNAL FUNCTIONS

function [pattern, theta] = compute_pattern(d, lambda, theta_deg, N, theta_int_deg)
    m = (0 : N-1)';
    A = exp(2j * pi * d * m * sin(deg2rad([theta_deg theta_int_deg])) / lambda);

    sigma_n = 1e-5;
    R = A * A' + sigma_n * eye(N);
    
    w = R^-1 * A(:,1) / (A(:,1)' * R^-1 * A(:,1));

    S = 3601;
    theta = linspace(-pi/2, pi/2, S);
    
    Psi = 2 * pi * d * (sin(theta) - sin(deg2rad(theta_deg))) / lambda;
    V = exp(1j * m * Psi);

    pattern = w' * V;
end

