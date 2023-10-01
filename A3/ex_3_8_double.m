close all
clear
clc

filename = "toucan.jpeg";
img = imread(filename);

figure
subplot(2, 2, 1)
imshow(img)
title("Full image")
img = double(img);
subplot(2, 2, 2)
imshow(img(:,:,1), [])
title("Red")
subplot(2, 2, 3)
imshow(img(:,:,2), [])
title("Green")
subplot(2, 2, 4)
imshow(img(:,:,3), [])
title("Blue")

tic
img_ycc = my_rgb2ycbcr_fast_double(img);
toc

img_ycc_matlab = rgb2ycbcr(img);

tic
img_ycc_slow = my_rgb2ycbcr_slow_double(img);
toc

figure
subplot(3, 3, 1)
imshow(img_ycc(:,:,1), [])
title("Y (fast)")
subplot(3, 3, 2)
imshow(img_ycc(:,:,2), [])
title("Cr (fast)")
subplot(3, 3, 3)
imshow(img_ycc(:,:,3), [])
title("Cb (fast)")
subplot(3, 3, 4)
imshow(img_ycc_slow(:,:,1), [])
title("Y (slow)")
subplot(3, 3, 5)
imshow(img_ycc_slow(:,:,2), [])
title("Cr (slow)")
subplot(3, 3, 6)
imshow(img_ycc_slow(:,:,3), [])
title("Cb (slow)")
subplot(3, 3, 7)
imshow(img_ycc_matlab(:,:,1), [])
title("Y (matlab)")
subplot(3, 3, 8)
imshow(img_ycc_matlab(:,:,2), [])
title("Cr (matlab)")
subplot(3, 3, 9)
imshow(img_ycc_matlab(:,:,3), [])
title("Cb (matlab)")





