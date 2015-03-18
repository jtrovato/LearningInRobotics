function R = getR(roll, pitch, yaw)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    Ryaw = [cos(yaw)  -sin(yaw) 0; sin(yaw)  cos(yaw) 0; 0 0 1];
    Rroll = [1 0 0; 0 cos(roll), -sin(roll); 0 sin(roll) cos(roll)];
    Rpitch = [cos(pitch) 0 sin(pitch);0 1 0;-sin(pitch) 0 cos(pitch)];
    R = Ryaw*Rpitch*Rroll; %rotation matix

end

