%% collect data
imdir = dir('./ESE650 P1');
num_imgs = length(imdir)-2;

barrel_pixels = [];
other_pixels = [];

for i=3:6%num_imgs+2
    imrgb = imread(['ESE650 P1/', imdir(i).name]); %read in image
    im_els = numel(imrgb(:,:,1));
    imycbcr = rgb2ycbcr(imrgb);% convert image to YCbCr color space
    imhsv = rgb2hsv(imrgb);
    barrel_mask = roipoly(imrgb); close;
    
    
    barrel_inds = find(barrel_mask);
    barrel_pixels = [barrel_pixels; imycbcr(barrel_inds), imycbcr(barrel_inds+im_els), imycbcr(barrel_inds+2*im_els)];
    %[other_i, other_j] = ind2sub(size(imrgb(:,:,1)), find(barrel_mask==0));
    other_inds = find(barrel_mask==0);
    other_pixels = [other_pixels; imycbcr(other_inds), imycbcr(other_inds+im_els), imycbcr(other_inds+2*im_els)];
%     %need to use linear indexing with images to avoid running out of memory
%     %[barrel_i, barrel_j] = ind2sub(size(imrgb(:,:,1)), find(barrel_mask));
%     barrel_inds = find(barrel_mask);
%     barrel_pixels = [barrel_pixels; imycbcr(barrel_inds), imycbcr(barrel_inds+im_els), imycbcr(barrel_inds+2*im_els)];
%     %[other_i, other_j] = ind2sub(size(imrgb(:,:,1)), find(barrel_mask==0));
%     other_inds = find(barrel_mask==0);
%     other_pixels = [other_pixels; imycbcr(other_inds), imycbcr(other_inds+im_els), imycbcr(other_inds+2*im_els)];

%% Plot the Data (just for me)
figure();
scatter(other_pixels(1:100000,2), other_pixels(1:100000,3), 'b.');
hold on;
scatter(barrel_pixels(1:10000,2), barrel_pixels(1:10000,3), 'r.');
%%
figure();
scatter3(barrel_pixels(1:100000,1), barrel_pixels(1:100000,2), barrel_pixels(1:100000,3), 'r.');
hold on;
xlabel('Y');
ylabel('Cb');
zlabel('Cr');
scatter3(mean(barrel_pixels(1:100000,1)), mean(barrel_pixels(1:100000,2)), mean(barrel_pixels(1:100000,3)), 'g*');
scatter3(other_pixels(1:100000,1), barrel_pixels(1:100000,2), barrel_pixels(1:100000,3), 'b.');
scatter3(mean(other_pixels(1:100000,1)), mean(other_pixels(1:100000,2)), mean(other_pixels(1:100000,3)), 'y*');
