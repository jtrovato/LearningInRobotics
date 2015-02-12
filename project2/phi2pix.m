function u = phi2pix( phi )
px_angle = linspace(-90, 90, 960);
[~,u] = min(abs(px_angle-phi));
end

