%function [p3, p2, s1, s2] = rotplot(R)
% This is a simple function to plot the orientation
% of a 3x3 rotation matrix R in 3-D
% You should modify it as you wish for the project.


lx = 3.0;
ly = 1.5;
lz = 1.0;

x = .5*[+lx -lx +lx -lx +lx -lx +lx -lx;
        +ly +ly -ly -ly +ly +ly -ly -ly;
        +lz +lz +lz +lz -lz -lz -lz -lz];

xp = R*x;
ifront = [1 3 7 5 1];
iback = [2 4 8 6 2];
itop = [1 2 4 3 1];
ibottom = [5 6 8 7 5];

x = plot3([0 1],[0 0],[0 0], 'r');
hold on;
y = plot3([0 0],[0 1],[0 0], 'r');
z = plot3([0 0],[0 0],[0 1], 'r');
hold off
% p3 = plot3(xp(1,itop), xp(2,itop), xp(3,itop), 'k-');
% hold on;
% p2 = plot3(xp(1,ibottom), xp(2,ibottom), xp(3,ibottom), 'k-');
% s1 = patch(xp(1,ifront), xp(2,ifront), xp(3,ifront), 'b');
% s2 = patch(xp(1,iback), xp(2,iback), xp(3,iback), 'r');
% hold off;

axis equal;
axis([-2 2 -2 2 -2 2]);

drawnow
