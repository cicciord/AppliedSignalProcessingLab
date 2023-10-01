close all
clear
clc

% variables definitions
T = 1;
f = 1 : 1 : 100;
phi = rand(1, 100) .* (2 * pi);

% axis definition
Fs = 20 * f(100);
Ts = 1 / Fs;
t = 0 : Ts : 1 - Ts;

s = zeros(100, Fs);
for i = 1 : 100
    s(i, :) = sin(2 * pi * f(i) .* t + phi(i));
end

s_a = sum(s(1:6, :));
s_b = sum(s(1:24, :));
s_c = sum(s(1:96, :));

% plots
figure

subplot(2, 3, 1)
plot(t, s(6,:));
title("Sine signal @6Hz");
ylabel("s_6(t)");
grid on

subplot(2, 3, 2)
plot(t, s(24, :));
title("Sine signal @24Hz");
xlabel("time [s]");
ylabel("s_{24}(t)");
grid on

subplot(2, 3, 3)
plot(t, s(96, :));
title("Sine signal @96Hz");
ylabel("s_{96}(t)");
grid on

subplot(2, 3, 4)
plot(t, s_a);
title("Sum from 1Hz to 6Hz");
ylabel("s_a(t)");
grid on

subplot(2, 3, 5)
plot(t, s_b);
title("Sum from 1Hz to 24Hz");
xlabel("time [s]");
ylabel("s_b(t)");
grid on

subplot(2, 3, 6)
plot(t, s_c);
title("Sum from 1Hz to 96Hz");
ylabel("s_c(t)");
grid on
