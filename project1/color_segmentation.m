% ESE650 Project 1 - Color Segementation
close all; clear all;

%% Preprocess data (toy data)

%% generate toy dataset:
% test_mu = [5;5;5];
% test_mu2 = [1,1,1];
% test_sigma = 0.1*eye(3);
% data = mvnrnd(test_mu', test_sigma, 1005);
% data = [data; mvnrnd(test_mu2', test_sigma, 1005)];
% figure();
% plot3(data(:,1), data(:,2), data(:,3), '+')
% grid on;
% axis equal;
% y = data;

%% load labeled data

% Load image data
imdir = dir('./ESE650 P1');
num_imgs = length(imdir)-2;

%shuffle data, split train and test
split = 0.2;
img_inds = randperm(num_imgs)+2;
train_ims = img_inds(1:round(num_imgs*split));
test_ims = img_inds(round(num_imgs*split)+1:end);
% All image data
load('pixel_data.mat');
y_b = double(barrel_pixels);
y_o = double(other_pixels(1:8:end,:));

%% Preprocess data (Images)

% Load image data
imdir = dir('./ESE650 P1');
num_imgs = length(imdir)-2;

%shuffle data, split train and test
split = 0.05;
img_inds = randperm(num_imgs)+2;
train_ims = img_inds(1:round(num_imgs*split));
test_ims = img_inds(round(num_imgs*split)+1:end);

barrel_pixels = [];
other_pixels = [];
%label training images
for i=train_ims
    imrgb = imread(['ESE650 P1/', imdir(i).name]); %read in image
    im_els = numel(imrgb(:,:,1));
    imycbcr = rgb2ycbcr(imrgb);% convert image to YCbCr color space
    barrel_mask = roipoly(imrgb); close;
    barrel_inds = find(barrel_mask);
    barrel_pixels = [barrel_pixels; imycbcr(barrel_inds), imycbcr(barrel_inds+im_els), imycbcr(barrel_inds+2*im_els)];
    other_inds = find(barrel_mask==0);
    other_pixels = [other_pixels; imycbcr(other_inds), imycbcr(other_inds+im_els), imycbcr(other_inds+2*im_els)];
end
[n ,d] = size(barrel_pixels);
y_b = double(barrel_pixels);


%% Plot the Data (just for me)
figure();
scatter(other_pixels(:,2), other_pixels(:,3), 'b.');
hold on;
scatter(barrel_pixels(:,2), barrel_pixels(:,3), 'r.');

%% EM for GMM fitting
% this algorithm is used to find the optimal parameters of the gaussians
% that make up out GMMs. 

max_iters = 10000;
[n_b,d] = size(y_b);
[n_o,~] = size(y_o);
k = 3; %number of gaussians in GMM
delta = 0.0001; %threshold to detemine convergence

%initialization
rng(13);
w_b = 1/k*ones(1,k); % [1 x k]
w_o = 1/k*ones(1,k); % [1 x k]
mu_b = [[60, 120, 170]', [50, 120, 200]', [100, 140, 150]']; % [d x k]
mu_o = [[60,120,120]', [100,100,100]', [120,120,130]'];
%mu = ones(3);
SIGMA_b = reshape([255*rand*eye(d), 255*rand*eye(d), 255*rand*eye(d)], d,d,k); %stacked covariance matrices [d x d x k] matrices
SIGMA_o = reshape([255*rand*eye(d), 255*rand*eye(d), 255*rand*eye(d)], d,d,k); %stacked covariance matrices [d x d x k] matrices
%SIGMA = SIGMA/255
L_b = log_likelihood(y_b, w_b, mu_b, SIGMA_b); %scalar
L_o = log_likelihood(y_o, w_o, mu_o, SIGMA_o); %scalar


figure();
mu_image = reshape(mu_b', 1,k,d);
image(ycbcr2rgb(uint8(mu_image)));
drawnow;
for p=1:max_iters
    disp(['iteration:' num2str(p)]);
    %E-Step:
    gamma_b = assign_probs_to_data(y_b, w_b, mu_b, SIGMA_b); % predict the data labels (with probabilities)
    prob_sums_b = sum(gamma_b);
    gamma_o = assign_probs_to_data(y_o, w_o, mu_o, SIGMA_o); % predict the data labels (with probabilities)
    prob_sums_o = sum(gamma_o);
    %M-Step (update the model parameters):
    w_b = prob_sums_b./n_b; %new weights are accroding to probability mass assigned to each gaussian. Could do raw counts.
    w_o = prob_sums_o./n_o; %new weights are accroding to probability mass assigned to each gaussian. Could do raw counts.
    mu_b = (repmat((1./prob_sums_b), d, 1).*(gamma_b'*y_b)');
    mu_o = (repmat((1./prob_sums_o), d, 1).*(gamma_o'*y_o)');

    for j=1:k
        mean_centered_b = bsxfun(@minus, y_b, mu_b(:,j)');
        SIGMA_b(:,:,j) = 1/prob_sums_b(j)*bsxfun(@times, gamma_b(:,j), mean_centered_b)'*mean_centered_b;
        mean_centered_o = bsxfun(@minus, y_o, mu_o(:,j)');
        SIGMA_o(:,:,j) = 1/prob_sums_o(j)*bsxfun(@times, gamma_o(:,j), mean_centered_o)'*mean_centered_o;
    end
    mu_image = reshape(mu_b', 1,k,d);
    image(ycbcr2rgb(uint8(mu_image))); drawnow;
    
    %Check for convergence:
    Lnew_b = log_likelihood(y_b, w_b, mu_b, SIGMA_b);
    disp(['barrel log likelihood diff =', num2str(abs(Lnew_b-L_b))]);
    Lnew_o = log_likelihood(y_o, w_o, mu_o, SIGMA_o);
    disp(['other log likelihood diff =', num2str(abs(Lnew_o-L_o))]);
    if abs(Lnew_b-L_b) < delta && abs(Lnew_o-L_o) < delta
        disp('EM converged for other pixels')
        w_b
        mu_b
        SIGMA_b
        w_o
        mu_o
        SIGMA_o
        break;
    end
    L_b = Lnew_b;
    L_o = Lnew_o;
end

%% Color Segmentation using GGM model found above
% Estimate probabilites of test image pixels in the test image
for cur_im=3:num_imgs+2
    test_im_rgb = imread(['ESE650 P1/', imdir(cur_im).name]);
    num_pix = size(test_im_rgb,1)*size(test_im_rgb,2);
    test_im = rgb2ycbcr(test_im_rgb);
    y_test = reshape(test_im, [num_pix,3]);

    [n,d] = size(y_test);
    probs_r = sum(repmat(w_b, n, 1).*gaussian_model(double(y_test), mu_b, SIGMA_b), 2);%find pixels with a high probability of being generated by a gaussian
    probs_o = sum(repmat(w_o, n, 1).*gaussian_model(double(y_test), mu_o, SIGMA_o), 2);%find pixels with a high probability of being generated by a gaussian
    red_inds = find(probs_r > probs_o);

    red_mask = zeros(size(test_im(:,:,1)));
    red_mask(red_inds)=1;
    
    %Some morphology to clean up the mask
    SE = strel('square',15);
    red_mask2 = imerode(red_mask, SE);
    SE2 = strel('square',20);
    red_mask2 = imdilate(red_mask2, SE2);
    

%     figure(99);
%     subplot(2,1,1);  
%     imshow(red_mask);
%     subplot(2,1,2);
%     imshow(red_mask2);
%     title(imdir(i).name);

    red_mask = red_mask2;
    find_barrel
    pause
end
        