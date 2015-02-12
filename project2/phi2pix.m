function u = phi2pix( phi )
px_angle = linspace(-pi/2, pi/2, 960);
[~,u] = min(abs(px_angle-phi));
end

