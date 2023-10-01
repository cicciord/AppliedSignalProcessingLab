function YCCimg = my_rgb2ycbcr_fast(RGBimg)
%MY_RGB2YCBCR_FAST Summary of this function goes here
%   Detailed explanation goes here
[M, N] = size(RGBimg, [1 2]);

% T = [0.299 0.587 0.114; -0.169 -0.331 0.5; 0.5 -0.419 -0.081];
T = [ ...
    65.481 128.553 24.966;...
    -37.797 -74.203 112; ...
    112 -93.786 -18.214];
offset = [16; 128; 128];

offset = repmat(offset, 1, M*N);

YCCimg_reshaped = T * reshape(im2double(RGBimg), M*N, 3)' + offset;

YCCimg = reshape(uint8(255 * mat2gray(YCCimg_reshaped))', M, N, 3);

% RGBimg = im2double(RGBimg);
% R = RGBimg(:,:,1);
% G = RGBimg(:,:,2);
% B = RGBimg(:,:,3);
% 
% Y = T(1,1)*R + T(1,2)*G + T(1,3)*B + offset(1);
% Cr = T(2,1)*R + T(2,2)*G + T(2,3)*B + offset(2);
% Cb = T(3,1)*R + T(3,2)*G + T(3,3)*B + offset(3);
% 
% YCCimg = cat(3, Y, Cr, Cb);
% YCCimg = uint8(255 * mat2gray(YCCimg));

end

