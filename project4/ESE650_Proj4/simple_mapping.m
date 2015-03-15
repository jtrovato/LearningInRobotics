%Simple Mapping Program
load lidar0.mat
load depth0.mat

figure();
pose = lidar{1}.pose; %[x,y,theta]
x=pose(1); y=pose(2); theta = pose(3);
T = [x+0.7, y; x-0.3, y-0.3; x-0.3, y+0.3;x+0.7, y];
h = plot(T(:,1),T(:,2));
axis([-10,10,-10,10]);
grid on;

for i=10:10:numel(lidar)
    % remove noisy data out of valid range
    lidar{i}.scan(find(lidar{i}.scan > 30)) = 0;
    
    
    pose = lidar{i}.pose; %[x,y,theta]
    x=pose(1); y=pose(2); theta = pose(3);
    T = [0.7, 0; -0.3, -0.3; -0.3, 0.3;0.7, 0];
    rot = [cos(theta)  -sin(theta) ; sin(theta)  cos(theta)];
    T=(rot*T')';
    T(:,1) = T(:,1)+x;
    T(:,2) = T(:,2)+y;
    set(h, 'XDATA', T(:,1), 'YDATA', T(:,2));
    
    
    pause(0.025);
    i
end