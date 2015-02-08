% %% Load image data
imdir = dir('./ESE650 P1');
num_imgs = length(imdir)-2;
% load('pixel_data.mat');
% barrel_pixels = barrel_pixels(1:10000, :);
% other_pixels = other_pixels(1:100000, :);
%% Color Classification
% classify each pixel in the test image as barrel color or not barrel color

%% thresholding in YCbCr
test_im_num = 6;
red_coordinates = [];
test_im_rgb = imread(['ESE650 P1/', imdir(test_im_num).name]);
test_im = rgb2ycbcr(test_im_rgb);

[red_coordinates_x, red_coordinates_y] = ind2sub(size(test_im), find(test_im(:,:,2)<140 & test_im(:,:,2)>110  & test_im(:,:,3)<250 & test_im(:,:,3)>150));

red_mask = zeros(size(test_im(:,:,1)));
red_mask(find(test_im(:,:,2)<140 & test_im(:,:,2)>110  & test_im(:,:,3)<250 & test_im(:,:,3)>150)) = 1;
figure();
imshow(red_mask);
title(imdir(test_im_num).name);
