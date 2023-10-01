close all
clear
clc

filename = "liberty.jpeg";
img = imread(filename);
img = rgb2gray(img);

figure
imshow(img)

[M, N] = size(img, [1 2]);

P = 2^nextpow2(2*max(M, N));

IMG = fftshift(fft2(img, P, P));
f = -P/2 : P/2 - 1/2;

D0sq = (P/16)^2;
[U, V] = meshgrid(f);
Dsq = U.^2 + V.^2;

H = exp(-(Dsq)./(2*D0sq));

figure
mesh(f, f, abs(H))

IMG_filt = IMG .* H;

img_filt = ifft2(IMG_filt);
img_filt = uint8(abs(img_filt(1:M, 1:N)));

figure
imshow(img_filt)

figure
mesh(f, f, log10(1+abs(IMG)))
figure
mesh(f, f, log10(1+abs(IMG_filt)))

figure
subplot(1, 2, 1)
imshow(log10(1+abs(IMG)), [])
subplot(1, 2, 2)
imshow(log10(1+abs(IMG_filt)), [])

imwrite(img_filt, "liberty_filtered.jpeg")



