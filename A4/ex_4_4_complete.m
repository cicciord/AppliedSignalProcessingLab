close all
clear
clc

upa_pattern_capon(16, 16, 90, 0, [86 85 80 100 105], [4 20 5 -15 15], 1e-5, 0)
upa_pattern_capon(16, 16, 90, 0, [86 85 80 100 105], [4 20 5 -15 15], 1e-5, 1)

upa_pattern_capon(16, 16, 90, 0, [88 85 80 100 105], [2 20 5 -15 15], 1e-5, 0)
upa_pattern_capon(16, 16, 90, 0, [88 85 80 100 105], [2 20 5 -15 15], 1e-5, 1)

upa_pattern_capon(16, 16, 90, 0, [88 85 80 100 105], [2 20 5 -15 15], 1e3, 0)
upa_pattern_capon(16, 16, 90, 0, [88 85 80 100 105], [2 20 5 -15 15], 1e3, 1)


function [] = upa_pattern_capon(Nz, Ny, theta_1deg, phi_1deg, theta_intfdeg, phi_intfdeg, sigman2, directive)

% Parameters with spacing d=lambda\2
N_tot = Nz*Ny;

theta_1 = deg2rad(theta_1deg);

phi_1 = deg2rad(phi_1deg);

theta_intf = deg2rad(theta_intfdeg);
phi_intf = deg2rad(phi_intfdeg);
n_intf = length(theta_intf);


% array in zy plane
Tz = -Nz/2+1 : Nz/2;
Tz = 2 * Tz -1;

Ty = -Ny/2+1 : Ny/2;
Ty = 2 * Ty -1;

T = Tz'*ones(1,Ny)./Nz+1i*ones(Nz,1)*Ty/Ny;
T = T(:)./sqrt(N_tot);

% antenna directivity function
if directive == 1
    dir = @(th, ph) 0.25 * (1 - cos(2 * th)) .* (1 + cos(ph));
else
    dir = @(th, ph) 1;
end

% beamforming filter
m = (0 : Nz-1)';
n = (0 : Ny-1)';

epsilon_z = pi * cos([theta_1 theta_intf]);
Az = exp(1j * m * epsilon_z);

epsilon_y = pi * sin([theta_1 theta_intf]) .* sin([phi_1 phi_intf]);
Ay = exp(1j * n * epsilon_y);

A = zeros(N_tot, n_intf+1);
for i = 1 : n_intf+1
    A(:, i) = kron(Az(:, i), Ay(:, i));
end
Dir = dir([theta_1 theta_intf], [phi_1 phi_intf]);
A = Dir .* A;

R = A * A' + sigman2 * eye(N_tot);

w = R^-1 * A(:,1) / (A(:,1)' * R^-1 * A(:,1));

% UPA pattern
angles_1 = 361;
theta = linspace(0, pi, angles_1);

angles_2 = 721;
phi = linspace(-pi, pi, angles_2);

V = zeros(N_tot, angles_1, angles_2);
for i = 1:angles_1
    epsilon_z = pi * cos(theta(i));
    az = exp(1j * epsilon_z .* m);
    for j = 1:angles_2
        epsilon_y = pi * sin(theta(i)) .* sin(phi(j));
        ay = exp(1j * epsilon_y .* n);

        V(:, i, j) = kron(az, ay);
    end
end

pattern = zeros(angles_1, angles_2);
for i = 1 : angles_1
    pattern(i, :) = w' * squeeze(V(:, i, :));
end

% directive case
[Phi, Theta] = meshgrid(phi, theta);
if directive == 1
    pattern = dir(Theta, Phi) .* pattern;
end

% Plot the pattern
r = abs(pattern);
X = r .* sin(Theta) .* cos(Phi);
Y = r .* sin(Theta) .* sin(Phi);
Z = r .* cos(Theta);
r_C = sqrt(X.^2 + Y.^2 + Z.^2);

r_steer = 1;
x_steer = r_steer .* sin(theta_1) .* cos(phi_1);
y_steer = r_steer .* sin(theta_1) .* sin(phi_1);
z_steer = r_steer .* cos(theta_1);

r_intf = max(max(abs(pattern)));
x_intf = r_intf .* sin(theta_intf) .* cos(phi_intf);
y_intf = r_intf .* sin(theta_intf) .* sin(phi_intf);
z_intf = r_intf .* cos(theta_intf);

figure
mesh(X, Y, Z, r_C, FaceAlpha=0.5, EdgeAlpha=0.5)
axis equal
view(80, 5)
hold on
hidden on
plot3(zeros(N_tot,1), imag(T), real(T), 'ks', 'MarkerFaceColor', [0.75 0.75 0.75]);
hold on
plot3([0 x_steer], [0 y_steer], [0 z_steer], '--r', LineWidth=2)
text2 = "UE el. (" + num2str(theta_1deg) + char(176) + ", " + num2str(phi_1deg) + char(176) + ")";
text(x_steer, y_steer, z_steer, text2, "FontSize", 12)
hold on
plot3(x_steer, y_steer, z_steer, 'r.', 'MarkerSize', 8)
if theta_1 ~= pi/2 || phi_1 ~= 0
    r_broad = max(max(abs(pattern)));
    x_broad = r_broad .* sin(pi/2) .* cos(0);
    y_broad = r_broad .* sin(pi/2) .* sin(0);
    z_broad = r_broad .* cos(pi/2);
    plot3([0 x_broad], [0 y_broad], [0 z_broad], '--r', LineWidth=2)
    text1 = "Broadside";
    text(r_broad, 0, 0, text1, "FontSize", 12)
end
hold on
plot3([zeros(1, n_intf); x_intf], [zeros(1, n_intf); y_intf], [zeros(1, n_intf); z_intf], '--k')
text3 = strings(1, n_intf);
for i = 1:n_intf
    text3(i) = "(" + num2str(theta_intfdeg(i)) + char(176) + ", " + num2str(phi_intfdeg(i)) + char(176) + ")";
end
text(x_intf, y_intf, z_intf, text3, "FontSize", 12)
title("Capon beamforming, UPA, N_z=" + num2str(Nz) + ", N_y=" + num2str(Ny) + ", DoA (\theta=" + num2str(theta_1deg) + ", \phi=" + num2str(phi_1deg) + ")", FontSize=12)
xlabel("x")
ylabel("y")
zlabel("z")

end



