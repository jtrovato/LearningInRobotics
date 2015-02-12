function v = theta2pix( theta )
px_angle = linspace(-pi, pi, 1920);
[~,v] = min(abs(px_angle-theta));
end

