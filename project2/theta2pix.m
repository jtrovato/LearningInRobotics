function v = theta2pix( theta )
px_angle = linspace(-180, 180, 1920);
[~,v] = min(abs(px_angle-theta));
end

