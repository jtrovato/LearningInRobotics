function [x,y,z] = rot2euler(r)
%remove the thrid dimension for 2D conversion
x = atan2(r(2,3,:),r(3,3,:));
y = atan2(-r(1,3,:), sqrt(r(2,3,:).^2 + r(3,3,:).^2));
z = atan2(r(1,2,:),r(1,1,:));

end