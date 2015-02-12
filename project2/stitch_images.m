%function h = stich_images(cam, rots, h)
% stitch images together given corresponding homographies and images
for i=1:length(rots)
    I = cam(:,:,:,i);
    R = rots(:,:,i);
    imw = 320;
    imh = 240;
    FOV = 60; %degreees
    ppd = imw/FOV;
    pano = zeros(180*ppd, 360*ppd, 3); 
    
    vector_image = ones(imh, imw, 3);
    vector_image_cart = ones(imh, imw, 3);
    vector_image_trans = ones(imh, imw, 3);
    vector_image_trans_cart = ones(imh, imw, 3);
    thetas = repmat(linspace(-pi/6,pi/6,imw), imh,1) ;
    phis = repmat(linspace(-pi/4/2,pi/4/2,imh)',1,imw);
    vector_image(:,:,1) = thetas;
    vector_image(:,:,2) = phis;
    for u=1:imh
        for v=1:imw
            [x,y,z] = sph2cart(vector_image(u,v,1), vector_image(u,v,2), vector_image(u,v,3));
            vector_image_cart(u,v,:) = reshape([x,y,z],1,1,3);
            vector_image_trans_cart(u,v,:) = R*[vector_image_cart(u,v,1); vector_image_cart(u,v,2); vector_image_cart(u,v,3)];
            [t,p,r] = cart2sph(vector_image_trans_cart(u,v,1), vector_image_trans_cart(u,v,2), vector_image_trans_cart(u,v,3));
            vector_image_trans(u,v,:) = reshape([t,p,r],1,1,3);
            pano(theta2pix(t), phi2pix(p),:) = I(u,v,:);
        end
    end
    
    figure();
    imshow(pano);
end
    

%end
