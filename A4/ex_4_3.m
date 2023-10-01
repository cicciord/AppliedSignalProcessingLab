close all
clear
clc

% Parameters
c = 3e8;
Fc = 6e9;

lambda = c/Fc;

Nz = 8;
Ny = 8;

d_z = lambda/2;
d_y = lambda/2;

theta_1_deg = 90;
theta_1 = deg2rad(theta_1_deg);

phi_1_deg = 0;
phi_1 = deg2rad(phi_1_deg);

directive = 0;

% array in zy plane
Tz = -Nz/2+1 : Nz/2;
Tz = (2 * Tz - 1) / 3;

Ty = -Ny/2+1 : Ny/2;
Ty = 2 * Ty - 1;

T = Tz' * ones(1, Ny) ./ Nz + 1j * ones(Nz, 1) * Ty / Ny;
T = T(:);

% antenna directivity function
if directive == 1
    dir = @(th, ph) 0.25 * (1 - cos(2 * th)) .* (1 + cos(ph));
end

% beamforming filter
epsilon_z = 2 * pi * d_z * cos(theta_1) / lambda;
m = (0 : Nz-1)';
az = exp(1j * epsilon_z * m);

epsilon_y = 2 * pi * d_y * sin(theta_1) * sin(phi_1) / lambda;
n = (0 : Ny-1)';
ay = exp(1j * epsilon_y * n);

w = kron(az, ay) / sqrt(Nz*Ny);

% UPA pattern
angles_1 = 361;
theta = linspace(0, pi, angles_1);

angles_2 = 721;
phi = linspace(-pi, pi, angles_2);

V = zeros(Nz*Ny, angles_1, angles_2);
for i = 1:angles_1
    epsilon_z = 2 * pi * d_z * cos(theta(i)) / lambda;
    az = exp(1j * epsilon_z .* m);
    for j = 1:angles_2
        epsilon_y = 2 * pi * d_y * sin(theta(i)) .* sin(phi(j)) / lambda;
        ay = exp(1j * epsilon_y .* n);

        V(:, i, j) = kron(az, ay);
    end
end

pattern = zeros(angles_1, angles_2);
for i = 1 : angles_1
    pattern(i, :) = w' * squeeze(V(:, i, :));
end

[Phi, Theta] = meshgrid(phi, theta);
% directive case
if directive == 1
    pattern = dir(Theta, Phi) .* pattern;
end

% Plot the pattern
r = abs(pattern);
X = r .* sin(Theta) .* cos(Phi);
Y = r .* sin(Theta) .* sin(Phi);
Z = r .* cos(Theta);
r_C = sqrt(X.^2 + Y.^2 + Z.^2);

r_broad = max(max(abs(pattern)));
x_broad = r_broad .* sin(pi/2) .* cos(0);
y_broad = r_broad .* sin(pi/2) .* sin(0);
z_broad = r_broad .* cos(pi/2);

figure
mesh(X, Y, Z, r_C, FaceAlpha=0.5, EdgeAlpha=0.5)
axis equal
view(140, 15)
hold on
hidden on
plot3(zeros(Nz*Ny,1), imag(T), real(T), 'ks', 'MarkerFaceColor', 'r')
hold on
plot3([0 x_broad], [0 y_broad], [0 z_broad], '--r', LineWidth=2)
text1 = "Broadside";
text(r_broad, 0, 0, text1, "FontSize", 12)
if theta_1 ~= pi/2 || phi_1 ~= 0
    r_steer = max(max(abs(pattern)));
    x_steer = r_steer .* sin(theta_1) .* cos(phi_1);
    y_steer = r_steer .* sin(theta_1) .* sin(phi_1);
    z_steer = r_steer .* cos(theta_1);
    hold on
    plot3([0 x_steer], [0 y_steer], [0 z_steer], '--r', LineWidth=2)
    text2 = "UE el. (" + num2str(theta_1_deg) + char(176) + ", " + num2str(phi_1_deg) + char(176) + ")";
    text(x_steer, y_steer, z_steer, text2, "FontSize", 12)
end
title("Conventional beamforming, UPA, N_z=" + num2str(Nz) + ", N_y=" + num2str(Ny) + ", DoA (\theta=" + num2str(theta_1_deg) + ", \phi=" + num2str(phi_1_deg) + ")", FontSize=12)
xlabel("x")
ylabel("y")
zlabel("z")












