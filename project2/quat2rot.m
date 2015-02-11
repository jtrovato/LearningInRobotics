function rot = quat2rot(q)
alpha_w = 2*acos(q(1));
e_w = 1/(sin(acos(q(1))));
rot = (alpha_w*e_w)*q(2:4);
end

