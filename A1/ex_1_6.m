close all
clear
clc

% rectangular signal definition
x = rand(1, 100);
y = rand(1, 150);

% linear convolution
z_custom_l = my_conv(x, y, 'l');
z_conv = conv(x, y);
axis = 0 : length(z_conv)-1;

figure
plot(axis, z_custom_l, LineWidth=0.8)
title("Linear convolution")
xlabel("index")
ylabel("value")
grid on
hold on
plot(axis, z_conv)
legend("my\_conv", "conv")

% circular convolution
z_custom_c = my_conv(x, y, 'c');
z_cconv = cconv(x, y, max(length(x), length(y)));
axis = 0 : length(z_cconv)-1;

figure
plot(axis, z_custom_c, LineWidth=0.8)
title("Circular convolution")
xlabel("index")
ylabel("value")
grid on
hold on
plot(axis, z_cconv)
legend("my\_conv", "cconv")
