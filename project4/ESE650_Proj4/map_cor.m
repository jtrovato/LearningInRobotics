function val = map_cor(ogrid, vects_w, s_off, xmin, ymin, xmax, ymax, res)
    %NOTE: vects_w must b ein homogeneous coordinates
    
    %rotate and shift lidar points according to offsets
    H = getH(s_off(1)*res,s_off(2)*res,0,0,0,s_off(3));
    vects_w_r = H*vects_w; 
    %filter values with small z or large x and y ou
    vects_w_r(:, vects_w_r(1,:) < xmin | vects_w_r(1,:) > xmax) = [];
    vects_w_r(:, vects_w_r(2,:) < ymin | vects_w_r(2,:) > ymax) = [];
    %map the points to o-grid coordinates
    %max([round((vects_w_r(1,:)-xmin)/res), round((vects_w_r(2,:)-ymin)/res)])
    %min([round((vects_w_r(1,:)-xmin)/res), round((vects_w_r(2,:)-ymin)/res)])

    xinds = round((vects_w_r(1,:)-xmin)/res);
    %xinds(xinds < xmin | xinds > xmax) = [];
    yinds = round((vects_w_r(2,:)-ymin)/res);
    %yinds(yinds < ymin | yinds > ymax) = [];

    
    inds = sub2ind(size(ogrid), xinds, yinds);
    val = sum(ogrid(inds));
end

