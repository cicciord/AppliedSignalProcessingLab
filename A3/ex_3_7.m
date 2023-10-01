close all
clear
clc

filename = "lena.bmp";
img = imread(filename);

figure
imshow(img)

[M, N] = size(img);

D4 = [0  8  2 10; 12  4 14  6; 3 11  1  9; 15  7 13  5];

D = repmat(D4, M, N);

L = repelem(img, 4, 4);
L = L/16;

img_dith = L > D;

figure
imshow(img_dith)

D8 = [4*D4 4*D4+2; 4*D4+3 4*D4+1];

D_2 = repmat(D8, M, N);

L_2 = repelem(img, 8, 8);
L_2 = L_2/4;

img_dith_2 = L_2 > D_2;

figure
imshow(img_dith_2)

