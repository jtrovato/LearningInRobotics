%function h = stich_images(cam, rots, h)
% stitch images together given corresponding homographies and images
close all;

pano = uint8(zeros(800, 1600, 3)); 
imw = 320;
imh = 240;
sx = 4.8e-3;
sy = 3.6e-3;
f = 3.67e-3;
FOV = 60; %degreees
ppd = imw/FOV;
hand = imshow(pano);
xscale = size(pano,2)/(2*pi);%px/rad
yscale = size(pano,1)/pi;%px/rad
[h,w] = size(pano(:,:,1));
imsz = h*w;

%create the vector image
pixx = reshape(repmat(1:imw, imh, 1), 1, imw*imh);
pixy = repmat(1:imh, 1, imw);
w_x = (pixx-imw/2)/imw*(sx/f);
w_y = (pixy-imh/2)/imh*(sy/f);

world = [ones(1,imh*imw); -w_x;-w_y];
R_0 = eye(3);
for i=1:length(cam_ts)-20
    t=cam_ts(i);
    [min_val, ind] = min(abs(imu_ts - t));
    I = cam(:,:,:,i);
    %R = rots(:,:,ind);
    R = kf_R(:,:,ind);
    
    %transform the vector image
    % R_bias = [-1 0 0; 0 0 1;0 1 0];
    t_r =R*world;
    azimuth = atan2(t_r(2,:), t_r(1,:)); %radians
    elevation = atan2(t_r(3,:), sqrt(t_r(1,:).^2 + t_r(2,:).^2));%radians
    px = ceil(xscale*azimuth + w/2);
    py = ceil(yscale*elevation + h/2);
    ind = sub2ind([h,w], py, px);
    ind_source = sub2ind([imh,imw],pixy,pixx);
    pano(ind) = I(ind_source);
    pano(ind+imsz) = I(ind_source+imh*imw);
    pano(ind+2*imsz) = I(ind_source+2*imh*imw);
    
    set(hand, 'CData', rot90(pano,2));
    drawnow;
    %pause(0.05);
end
    

%end
