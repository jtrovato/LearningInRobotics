function H = getH(x, y,z, roll, pitch, yaw)
R = getR(roll, pitch, yaw);
H = [[R, [x y z]'];[0 0 0 1]];

end

