%label minor axis lengths for lookup table 
figure();
dists = zeros(num_imgs,1);
barrel_width_sum = zeros(10,1);
barrel_width_num = zeros(10,1);
for i=3:num_imgs+2
    imrgb = imread(['ESE650 P1/', imdir(i).name]); %read in image
    imshow(imrgb);
    title(imdir(i).name);
    h = imdistline(gca);
    api = iptgetapi(h);
    fcn = makeConstrainToRectFcn('imline',...
                                  get(gca,'XLim'),get(gca,'YLim'));
    api.setDragConstraintFcn(fcn);  
    pause
    dist = api.getDistance();
    dists(i-2) = dist;
    im_label = imdir(i).name(1:2);
    if im_label(2) == '.'
        im_label = im_label(1);
    end
    im_label = str2num(im_label);
    barrel_width_sum(im_label) = barrel_width_sum(im_label) + dist
    barrel_width_num(im_label) = barrel_width_num(im_label) + 1
    
end
barrel_width_lut = barrel_width_sum./barrel_width_num;
