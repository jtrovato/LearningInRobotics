function rot = quat2rot( q, deltat)
alpha_w = 2*arccos(q(0));
e_w = 1/(sin(arccos(q(0))));
rot = (alpha_w*e_w/deltat)*q(2:4);
end

