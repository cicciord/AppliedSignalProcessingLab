close all
clear
clc

b1 = [1 1.3 0.49 -0.013 -0.029];
a1 = [1 -0.4326 -1.6656 0.1253 0.2877];

b2 = [1 -1.4 0.24 0.3340 -0.1305];
a2 = [1 0.5913 -0.6436 0.3803 -1.0091];

b3 = [0.0725 0.22 0.4085 0.4883 0.4085 0.22 0.0725];
a3 = [1 -0.5835 1.7021 -0.8477 0.8401 -0.2823 0.0924];

z1 = roots(b1);
p1 = roots(a1);

z2 = roots(b2);
p2 = roots(a2);

z3 = roots(b3);
p3 = roots(a3);

if sum(abs(p1) >= 1) ~= 0
    disp("Filter 1 is unstable")
end

if sum(abs(p2) >= 1) ~= 0
    disp("Filter 2 is unstable")
end

if sum(abs(p3) >= 1) ~= 0
    disp("Filter 3 is unstable")
end

figure
zplane(z1, p1)
title("Zeros and Poles of filter 1")
xlabel("Re")
ylabel("Im")
grid on

figure
zplane(z2, p2)
title("Zeros and Poles of filter 2")
xlabel("Re")
ylabel("Im")
grid on

figure
zplane(z3, p3)
title("Zeros and Poles of filter 3")
xlabel("Re")
ylabel("Im")
grid on

N = 50;

X = [1 zeros(1, N-1)];
Y1 = filter(b1, a1, X);
Y2 = filter(b2, a2, X);
Y3 = filter(b3, a3, X);

figure
stem(Y1)
title("Impulse response of filter 1")
xlabel("samples")
ylabel("amplitude")
grid on

figure
stem(Y2)
title("Impulse response of filter 2")
xlabel("samples")
ylabel("amplitude")
grid on

figure
stem(Y3)
title("Impulse response of filter 3")
xlabel("samples")
ylabel("amplitude")
grid on
