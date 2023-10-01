function YCCimg = my_rgb2ycbcr_slow(RGBimg)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
[M, N] = size(RGBimg, [1 2]);

YCCimg = zeros(M, N, 3);
% T = [0.299 0.587 0.114; -0.169 -0.331 0.5; 0.5 -0.419 -0.081];
T = [ ...
    65.481 128.553 24.966;...
    -37.797 -74.203 112; ...
    112 -93.786 -18.214];
offset = [16; 128; 128];
for i = 1:M
    for j = 1:N
        YCCimg(i, j, :) = T * squeeze(im2double(RGBimg(i, j, :))) + offset;
    end
end

YCCimg = uint8(255 * mat2gray(YCCimg));

