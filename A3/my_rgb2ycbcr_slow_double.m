function YCCimg = my_rgb2ycbcr_slow_double(RGBimg)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
[M, N] = size(RGBimg, [1 2]);

YCCimg = zeros(M, N, 3);
T = [0.299 0.587 0.114; -0.169 -0.331 0.5; 0.5 -0.419 -0.081];

for i = 1:M
    for j = 1:N
        YCCimg(i, j, :) = T * squeeze(RGBimg(i, j, :));
    end
end

