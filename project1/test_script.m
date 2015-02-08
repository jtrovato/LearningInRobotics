
clear all; close all;  
option = 3;
test_dir_name = './ESE650_PROJ1_TEST/';

%% Option 1: If you really want to start from the begining:
if option == 1
    label_images;
    label_barrel_widths;
    EMGMMfit; 
    [barrel_dist_func, gof] = fit_barrel_dist_func(barrel_width_lut, 1:10);
end
%% Option 2: If you want to skip the labeling:
if option == 2

    load('pixel_data.mat');
    load('barrel_width_lut');
    EMGMMfit;
    [barrel_dist_func, gof] = fit_barrel_dist_func(barrel_width_lut, 1:10);
end
%% Option 3: or if you just want to load everything:
if option == 3
    %load trained color models
    load('color_model.mat');
    load('barrel_width_lut.mat'); %load the barrel widths
    load('fitted_func');
end

%% After everything is ready to go

% Load test image data 

imdir = dir(test_dir_name);
num_imgs = length(imdir)-2;  
for cur_im=3:num_imgs+2
    classify_colors; 
    find_barrel; %centroid and distance are store in variables x, y, d
    disp(['img name:', imdir(cur_im).name, '     x= ',num2str(x),'   y=',num2str(y),'   d=',num2str(d)]);
    pause;
end