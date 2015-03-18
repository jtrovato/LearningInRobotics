%Simple Mapping Program
close all;

verbose = 1;
if verbose
    load joints0.mat
    load lidar0.mat
    ts_l = zeros(numel(lidar), 1);
    for i=1:numel(lidar)
        ts_l(i) = (lidar{i}.t -t0); %make the lidar time on the same scale as the joint time
    end
    load depth0.mat
end

thetas = 0:0.25:270;
thetas = thetas*pi/180;
num_angs = length(thetas);
rots = zeros(2,2,num_angs);

xmin = -15;
xmax = 15;
ymin = -15;
ymax = 15;
zmin = -2;
zmax = 3;

% the initial plot of the map
figure();
z = 1.34;
iNeck = get_joint_index('Neck'); % head yaw
iHead = get_joint_index('Head'); % head pitch
%localization object
T = [0.7, 0, 0
    -0.3, -0.3, 0
    -0.3, 0.3, 0
    0.7, 0, 0]';

%rotation
% ax = eye(3);
% h1 = plot3([0 ax(1,1)],[0 ax(2,1)],[0 ax(3,1)],'r');
% hold on;
% h2 = plot3([0 ax(1,2)],[0 ax(2,2)],[0 ax(3,2)],'b');
% h3 = plot3([0 ax(1,3)],[0 ax(2,3)],[0 ax(3,3)],'g');
%tranlation
h = fill(T(1,:),T(2,:), 'g');
hold on;
axis([xmin,xmax,ymin,ymax,zmin,zmax]);
axis equal;
grid on;
%mapping object
%g = plot3(0, 0, 0, 'b.', 'MarkerSize', 1);


%2D occupancy grid:
res = 0.05;
inc = 1;
occ_thres = 0.8;
ogrid = 0.5*ones((xmax-xmin)*(1/res), (ymax-ymin)*(1/res));
o = plot(0, 0, 'b.', 'MarkerSize', 1);

for i=20000:10:numel(ts)

    %get data and match up times
    t = ts(i); %taken from joints which is relative time
    fprintf('time: %4.2f \r', t);
    %find correct lidar (these use a different time than the joints)
    [~,idx] = min(abs(ts_l - t));
    pose = lidar{idx}.pose; %[x,y,theta]
    rpy = lidar{idx}.rpy';
    %get the correct joint angles (these use the joint time)
    head_angles = [pos(i,iNeck), pos(i,iHead)]; %[yaw, pitch]

    
    %localization
    %rotation
    Rh2b = getR(0,head_angles(2), head_angles(1));
    Rb2w = getR(rpy(1),rpy(2),pose(3));
    R = Rb2w*Rh2b;
%    axr = R*ax;
%     set(h1, 'XDATA', [0 axr(1,1)/abs(axr(:,1))]+pose(1), 'YDATA', [0 axr(2,1)/abs(axr(:,1))]+pose(2), 'ZDATA', [0 axr(3,1)/abs(axr(:,1))]+z);
%     set(h2, 'XDATA', [0 axr(1,2)/abs(axr(:,2))]+pose(1), 'YDATA', [0 axr(2,2)/abs(axr(:,2))]+pose(2), 'ZDATA', [0 axr(3,2)/abs(axr(:,2))]+z);
%     set(h3, 'XDATA', [0 axr(1,3)/abs(axr(:,3))]+pose(1), 'YDATA', [0 axr(2,3)/abs(axr(:,3))]+pose(2), 'ZDATA', [0 axr(3,3)/abs(axr(:,3))]+z);

    Tr=(R*T);
    %translation
    Tr(1,:) = Tr(1,:)+pose(1);
    Tr(2,:) = Tr(2,:)+pose(2);
    Tr(3,:) = Tr(3,:)+z;
    set(h, 'XDATA', Tr(1,:), 'YDATA', Tr(2,:));
    
    %MAPPING --- OCCUPANCY GRID IMPLEMENTATION
    
    % remove noisy data out of valid range
    lidar{idx}.scan(lidar{idx}.scan > 20) = 0;
    
    dists = lidar{idx}.scan;    
    vects_h = [dists.*cos(thetas-2.3562); dists.*sin(thetas-2.3562); zeros(1, num_angs)];
    vects_w = R*vects_h;
    vects_w(1,:) = vects_w(1,:)+pose(1);
    vects_w(2,:) = vects_w(2,:)+pose(2);
    vects_w(3,:) = vects_w(3,:)+1.4;
    %set(g, 'XDATA', vects_w(1,:), 'YDATA', vects_w(2,:), 'ZDATA', vects_w(3,:));
    %g = plot3(vects_w(1,:), vects_w(2,:), vects_w(3,:), 'b.', 'MarkerSize', 1);

    %2d ogrid
    %filter values with small z or large x and y ou
    vects_w(:, vects_w(3,:) < .5) = [];
    vects_w(:, vects_w(1,:) < xmin | vects_w(1,:) > xmax) = [];
    vects_w(:, vects_w(2,:) < ymin | vects_w(2,:) > ymax) = [];



    
    ogrid(sub2ind(size(ogrid), round((vects_w(1,:)-xmin)/res), round((vects_w(2,:)-ymin)/res))) = ...
        ogrid(sub2ind(size(ogrid), round((vects_w(1,:)-xmin)/res), round((vects_w(2,:)-ymin)/res))) + 5*inc;
    between = [];
    xs = round((pose(1)-xmin)/res);
    ys = round((pose(2)-ymin)/res);
    for j=1:100:numel(vects_w(3,:))-100;

        xf = (vects_w(1,j:j+99)-xmin)/res;
        yf = (vects_w(2,j:j+99)-ymin)/res;
        
        [x_between, y_between] = getMapCellsFromRay(double(xs),double(ys),double(xf),double(yf));
        x_between(x_between <= 0) = [];
        y_between(y_between <= 0) = [];
        %between = unique([x_between, y_between], 'rows');
        between = [x_between, y_between];
        ogrid(sub2ind(size(ogrid),between(:,1), between(:,2))) = ...
            ogrid(sub2ind(size(ogrid),between(:,1), between(:,2))) - inc;
    end
    [occupied_i, occupied_j] = find(siglim(ogrid) > occ_thres);
    set(o, 'XDATA', (occupied_i*res)+xmin, 'YDATA', (occupied_j*res)+ymin);
    if(mod(i, 50) == 0)
        pause(0.025);
    end
end