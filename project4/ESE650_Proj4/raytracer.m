function cells = raytracer(source, dest)
    pts1 = repmat([source; 1], 1, size(dest, 2));
    pts2 = [dest, ones(1, size(dest, 2))];
    ls = cross(pts1, pts2); %lines
    
    %TODO finish
    
end

