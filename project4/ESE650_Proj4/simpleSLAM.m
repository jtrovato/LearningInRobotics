%Simple Mapping Program
close all;
verbose =1;
if verbose
    load lidar0.mat
    load depth0.mat
    load joints0.mat
end

thetas = 0:0.25:270;
thetas = thetas*pi/180;
num_angs = length(thetas);
rots = zeros(2,2,num_angs);

xmin = -15;
xmax = 15;
ymin = -15;
ymax = 15;

% the initial plot of the map
figure();
cmap = colormap(autumn);
cmap(1,:) = [0 0 0];

%get data
pose = lidar{1}.pose; %[x,y,theta]
rpy = lidar{1}.rpy';
x = pose(1); y=pose(2);
iNeck = get_joint_index('Neck'); % head yaw
iHead = get_joint_index('Head'); % head pitch
t = ts(1);
head_angles = [pos(1,iNeck), pos(1,iHead)]; %[yaw, pitch]

%localization
T = [x+0.7, y, 0; x-0.3, y-0.3, 0; x-0.3, y+0.3, 0;x+0.7, y, 0];
%rotation
Rh2b = getR(0,head_angles(2), head_angles(1));
Rb2w = getR(rpy(1),rpy(2),pose(3));
R = Rb2w*Rh2b;
Tr=(R*T');
%tranlation
Tr(1,:) = Tr(1,:)+x;
Tr(2,:) = Tr(2,:)+y;
h = fill(Tr(1,:),Tr(2,:), 'g');
axis([xmin,xmax,ymin,ymax]);
grid on;
hold on;
%mapping
lidar{1}.scan(find(lidar{1}.scan > 20)) = 0;
dists = lidar{1}.scan;
[vects_h_x, vects_h_y] = pol2cart(thetas-135, dists);
vects_h = [dists.*cos(thetas-2.3562); dists.*sin(thetas-2.3562); zeros(1, num_angs)];
vects_w = R*vects_h;
vects_w(:, vects_w(3,:) < 0) = [];
vects_w(1,:) = vects_w(1,:)+x;
vects_w(2,:) = vects_w(2,:)+y;
%g = plot(vects_w(1,:), vects_w(2,:), 'b.', 'MarkerSize', 1);

%2D occupancy grid:
res = 0.05;
ogrid = zeros((xmax-xmin)*(1/res), (ymax-ymin)*(1/res));
ogrid(sub2ind(size(ogrid), round((vects_w(1,:)-xmin)/res), round((vects_w(2,:)-ymin)/res))) = 1;
[occupied_i, occupied_j] = find(ogrid == 1);
g = plot((occupied_i*res)+xmin, (occupied_j*res)+ymin, 'b.', 'MarkerSize', 1);

for i=10:10:numel(lidar)

    %get data
    pose = lidar{i}.pose; %[x,y,theta]
    rpy = lidar{i}.rpy';
    iNeck = get_joint_index('Neck'); % head yaw
    iHead = get_joint_index('Head'); % head pitch
    t = ts(i);
    head_angles = [pos(i,iNeck), pos(i,iHead)]; %[yaw, pitch]

    %localization
    %rotation
    Rh2b = getR(0,head_angles(2), head_angles(1));
    Rb2w = getR(rpy(1),rpy(2),pose(3));
    R = Rb2w*Rh2b;
    Tr=(R*T');
    %translation
    Tr(1,:) = Tr(1,:)+pose(1);
    Tr(2,:) = Tr(2,:)+pose(2);
    set(h, 'XDATA', Tr(1,:), 'YDATA', Tr(2,:));
    
    %mapping
    % remove noisy data out of valid range
    lidar{i}.scan(lidar{i}.scan > 20) = 0;
    dists = lidar{i}.scan;    
    vects_h = [dists.*cos(thetas-2.3562); dists.*sin(thetas-2.3562); zeros(1, num_angs)];
    vects_w = R*vects_h;
    vects_w(:, vects_w(3,:) < -.5) = [];
    vects_w(1,:) = vects_w(1,:)+pose(1);
    vects_w(2,:) = vects_w(2,:)+pose(2);
    %set(g, 'XDATA', vects_w(1,:), 'YDATA', vects_w(2,:));
    %g = plot(vects_w(1,:), vects_w(2,:), 'b.', 'MarkerSize', 1);

    %2d ogrid
    ogrid(sub2ind(size(ogrid), round((vects_w(1,:)-xmin)/res), round((vects_w(2,:)-ymin)/res))) = 1;
    [x_between, y_between] = getMapCellsFromRay(round(x-xmin)/res,round(y-ymin)/res, vects_w(1,:), vects_w(2,:));
    between = unique([x_between, y_between], 'rows');
    
    ogrid(sub2ind(between(:,1), between(:,2))) = 0;
    [occupied_i, occupied_j] = find(ogrid == 1);
    set(g, 'XDATA', (occupied_i*res)+xmin, 'YDATA', (occupied_j*res)+ymin);
    
    pause(0.025);
    i
end