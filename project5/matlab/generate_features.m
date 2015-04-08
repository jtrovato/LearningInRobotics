function F = generate_features(I)
%a script to generate features from the aerial image

Ig = rgb2gray(I);
[m,n] = size(Ig);


%HSV
%Ihsv = single(rgb2hsv(I));
%YCbCr
Iycbcr = rgb2ycbcr(I);
%Lab
Ilab = rgb2lab(I);


% %edge detection
% canny = single(edge(Ig, 'canny'));
% G = fspecial('gaussian',[5 5],2);
% canny_blur = imfilter(canny, G, 'same');
% 
% %textures - convolution filters
% %Irange_rgb = single(rangefilt(I));
% Irange_gr = single(rangefilt(Ig));
% Ientrop = single(entropyfilt(Ig));
% Istd = single(stdfilt(Ig, ones(5)));


%Hough transform

% %Kmeans
% k1=8;
% cim = reshape(double(Iycbcr),[], 3); 
% [cluster_idx1, cluster_center] = kmeans(cim,k1,'distance','sqEuclidean');
% pixel_labels1 = reshape(cluster_idx1, size(I(:,:,1)));
% imagesc(pixel_labels1);
% axis equal
% %binarize clusters
% cluster_bins1 = zeros(m,n,k1);
% for i=1:k1
%     cluster_bins1(:,:,i) = (pixel_labels1 == i);
% end

k2=5;
cim = reshape(Ilab(:,:,2:3),[],2); 
[cluster_idx2, cluster_center] = kmeans(cim,k2,'distance','sqEuclidean');
pixel_labels2 = reshape(cluster_idx2, size(I(:,:,1)));
imagesc(pixel_labels2);
axis equal
%binarize clusters
cluster_bins2 = zeros(m,n,k2);
for i=1:k2
    cluster_bins2(:,:,i) = (pixel_labels2 == i);
end

% d = k1+k2+5;
 F = zeros(m*n, k2+k2^2);
% F(:,1:k1) = reshape(cluster_bins1, [], k1);
% F(:,k1+1:k1+k2) = reshape(cluster_bins2, [], k2);
% F(:,end-4) = Irange_gr(:);
% F(:,end-3) = Ientrop(:);
% F(:,end-2) = Istd(:);
% F(:,end-1) = canny(:);
% F(:,end) = canny_blur(:);

F(:,1:k2) = reshape(cluster_bins2, [], k2);
for j=1:k2
    for i=1:k2
        F(:,k2+j+i) = F(:,j).*F(:,i);
    end
end

%maybe use a kernel here?


%normalize the features:
Fmean = mean(F,1);
Fstd = std(F,0,1);
F = bsxfun(@minus, F, Fmean);
F = bsxfun(@rdivide, F, Fstd);

save('feats.mat', 'F');

