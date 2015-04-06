function F = generate_features(I)
%a script to generate features from the aerial image

Ig = single(rgb2gray(I));
[m,n] = size(Ig);
d = 3+1+3+3+4+1;

%HSV
%Ihsv = single(rgb2hsv(I));
%YCbCr
Iycbcr = single(rgb2ycbcr(I));

%gradient map
[grad_mag, grad_dir] = imgradient(Ig);
grad_mag = single(grad_mag);
grad_dir = single(grad_dir);

%edge detection
canny = single(edge(Ig, 'canny'));

%connected components
thres = 70;
Ibin = single(zeros(size(Ig)));
Ibin(Ig > thres) = 1;
%CC = bwconncomp(Ig);

%textures - convolution filters
Irange_rgb = single(rangefilt(I));
Irange_gr = single(rangefilt(Ig));

%Hough transform

%Kmeans

F = single(zeros(m*n,d));
F(:,1:3) = reshape(I, [], 3);
F(:,4) = Ig(:);
F(:,5) = grad_mag(:);
F(:,6) = grad_dir(:);
F(:,7) = canny(:);
F(:,8) = Ibin(:);
F(:,9:11) = reshape(Iycbcr, [], 3);
F(:,12:14) = reshape(Irange_rgb, [], 3);
F(:,15) = Irange_gr(:);

%maybe use a kernel here?

Fmean = mean(F,1);
Fstd = std(F,0,1);
F = bsxfun(@minus, F, Fmean);
F = bsxfun(@rdivide, F, Fstd);

save('feats.mat', 'F');

