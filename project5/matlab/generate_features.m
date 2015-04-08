function F = generate_features(I)
%a script to generate features from the aerial image

Ig = rgb2gray(I);
[m,n] = size(Ig);


%HSV Thesholding
Ihsv = single(rgb2hsv(I));
roads = Ihsv(:,:,1) <.35 & Ihsv(:,:,1) >.29 & Ihsv(:,:,2) <.11 & Ihsv(:,:,2) >.08 & Ihsv(:,:,3) < 0.45 & Ihsv(:,:,3) > 0.15;
roads = imdilate(roads, ones(9));
roads = imerode(roads, ones(5));
roads = imclose(roads, ones(7));
sidewalks = Ihsv(:,:,1) <.3 & Ihsv(:,:,1) >.14 & Ihsv(:,:,2) <.2 & Ihsv(:,:,2) >.1 & Ihsv(:,:,3) < 0.9 & Ihsv(:,:,3) > 0.6;
trees = Ihsv(:,:,1) <.35 & Ihsv(:,:,1) >.29 & Ihsv(:,:,2) <0.45 & Ihsv(:,:,2) >.25 & Ihsv(:,:,3) < 0.45 & Ihsv(:,:,3) > 0.25;
trees = imclose(trees, ones(15));

%YCbCr
Iycbcr = rgb2ycbcr(I);
%Lab
Ilab = rgb2lab(I);

%edge detection
canny = single(edge(Ig, 'canny'));
G = fspecial('gaussian',[5 5],2);
canny_blur = imfilter(canny, G, 'same');
canny_dilate = imdilate(canny, ones(9));

% %textures - convolution filters
%Irange_rgb = single(rangefilt(I));
%Irange_gr = single(rangefilt(Ig));
%Ientrop = single(entropyfilt(Ig));
%Istd = single(stdfilt(Ig, ones(5)));




%region props

%Kmeans
k2=6;
cim = reshape(Ihsv(:,:,:),[],3); 
[cluster_idx2, cluster_center] = kmeans(cim,k2,'distance','sqEuclidean');
pixel_labels2 = reshape(cluster_idx2, size(I(:,:,1)));
pixel_labels2 = imerode(pixel_labels2, ones(5));
imagesc(pixel_labels2);
axis equal
%binarize clusters
cluster_bins2 = zeros(m,n,k2);
for i=1:k2
    cluster_bins2(:,:,i) = (pixel_labels2 == i);
end

d = k2+15+6;
F = zeros(m*n, d);
F(:,1:k1) = reshape(cluster_bins1, [], k1);
F(:,k1+16:k1+18) = reshape(Ihsv, [], 3);
F(:,end-5) = canny_dilate;
F(:,end-4) = canny(:);
F(:,end-3) = canny_blur(:) > 0.5;
F(:,end-2) = trees(:);
F(:,end-1) = sidewalks(:);
F(:,end) = roads(:);


 for j=1:k2
     for i=j:k2
         F(:,k1+j+i-1) = F(:,j).*F(:,i);
     end
 end


%maybe use a kernel here?
F = F-0.5;

% %normalize the features:
% Fmean = mean(F,1);
% Fstd = std(F,0,1);
% F = bsxfun(@minus, F, Fmean);
% F = bsxfun(@rdivide, F, Fstd);

save('feats.mat', 'F');

