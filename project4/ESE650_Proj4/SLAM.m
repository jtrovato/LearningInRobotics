%Simultaneous localization and mapping
close all;

%% load the data
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

%% initialization
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

z = 1.34;
iNeck = get_joint_index('Neck'); % head yaw
iHead = get_joint_index('Head'); % head pitch
T = [0.7, 0, 0
    -0.3, -0.3, 0
    -0.3, 0.3, 0
    0.7, 0, 0]';

% rotation
% ax = eye(3);
% h1 = plot3([0 ax(1,1)],[0 ax(2,1)],[0 ax(3,1)],'r');
% hold on;
% h2 = plot3([0 ax(1,2)],[0 ax(2,2)],[0 ax(3,2)],'b');
% h3 = plot3([0 ax(1,3)],[0 ax(2,3)],[0 ax(3,3)],'g');
%tranlation

%2D occupancy grid:
res = 0.05;
inc = 0.05;
occ_thres = 0.8;
ogrid = zeros((xmax-xmin)*(1/res), (ymax-ymin)*(1/res));

%particle filter (kinda)
%particle_priors_offsets = [];

gridx_offsets = -1:1;
gridy_offsets = -1:1;
H = getH(0, 0, 0, 0, 0, 0);
previous_pose = [0 0 0];
theta0 = lidar{1}.rpy(3);
s = [lidar{1}.pose(1:2), theta0]';

%plot all the things:
figure();
im_hand = imshow(siglim(ogrid));
hold on;
%axis([xmin,xmax,ymin,ymax,zmin,zmax]);
%axis equal;
p = plot(s(1),s(2),'r.', 'MarkerSize', 1);
%o = plot(0, 0, 'b.', 'MarkerSize', 1);
h = fill(T(1,:),T(2,:), 'g');
g = plot(0, 0, 'g.', 'MarkerSize', 1);

s_hist = s;
%% run in time
for i=20000:10:numel(ts)

    %get data and match up times
    t = ts(i); %taken from joints which is relative time
    fprintf('time: %4.2f \r', t);
    %find correct lidar (these use a different time than the joints)
    [~,idx] = min(abs(ts_l - t));
    pose = lidar{idx}.pose; %[x,y,theta]
    rpy = lidar{idx}.rpy';
    theta = rpy(3)-theta0;
    %get the correct joint angles (these use the joint time)
    head_angles = [pos(i,iNeck), pos(i,iHead)]; %[yaw, pitch]
    % remove noisy data out of valid range
    lidar{idx}.scan(lidar{idx}.scan > 20) = 0;
    dists = lidar{idx}.scan; 
    vects_l = [dists.*cos(thetas-2.3562); dists.*sin(thetas-2.3562); zeros(1, num_angs)];
    vects_w = H*[vects_l; ones(1, length(vects_l))];

    
    % LOCALIZATION -- PARTICLEish FILTER IMPLEMENTATION
    %generate "school of fish" particles
    ds = [pose(1:2) - previous_pose(1:2), theta - previous_pose(3)];
    
    theta_offsets = normrnd(0, abs(ds(3)), [1, 5]);
    particle_weights = zeros(length(gridy_offsets), length(gridx_offsets), length(theta_offsets));
    for v = 1:length(gridx_offsets)
        for u=1:length(gridy_offsets)
            for w =1:length(theta_offsets)
                particle_weights(v, u, w) = map_cor(ogrid, vects_w, [ds(1)+gridx_offsets(v),ds(2)+gridy_offsets(u),ds(3)+theta_offsets(w)], xmin, ymin, xmax, ymax, res);
            end
        end
    end
    
    [max_cor, max_ind] = max(reshape(particle_weights, numel(particle_weights), 1));
    [idx, idy, idz] = ind2sub(size(particle_weights), max_ind);
    fprintf('%f     %i ,%i, %i    ', range(theta_offsets), idx, idy, idz);
    master_particle = [s(1)+gridx_offsets(idx)*res, s(2)+gridy_offsets(idy)*res, s(3)+theta_offsets(idz)]';
    s = master_particle;
    
    %rotation and translation homographies (after the master particle)
    Rb2n = getR(0,head_angles(2), head_angles(1));
    Rw2b = getR(rpy(1),rpy(2),rpy(3));
    tw2b = [s(1), s(2), .93]'; %from world to cog
    tb2n = [0, 0, .40]'; %from imu to neck
    tn2l = [0, 0, .10]'; %from neck to lidar
    Hw2b = [[Rw2b, tw2b];[0 0 0 1]];
    Hb2n = [[Rb2n, tb2n];[0 0 0 1]];
    Hn2l = [[eye(3), tn2l];[0 0 0 1]];
    H = Hw2b*Hb2n*Hn2l; %build the world to lidar transformation
    
    Tr=(H*[T;ones(1, 4)]); %project the triangle robot to the proper location
    
    %MAPPING --- OCCUPANCY GRID IMPLEMENTATION  
    vects_w(:, vects_w(3,:) < .2) = [];%trim in z

    %g = plot3(vects_w(1,:), vects_w(2,:), vects_w(3,:), 'b.', 'MarkerSize', 1);

    %2d ogrid
    %filter values with small z or large x and y ou
    vects_w(:, vects_w(1,:) < xmin | vects_w(1,:) > xmax) = [];
    vects_w(:, vects_w(2,:) < ymin | vects_w(2,:) > ymax) = [];
    
    ogrid(sub2ind(size(ogrid), round((vects_w(1,:)-xmin)/res), round((vects_w(2,:)-ymin)/res))) = ...
        ogrid(sub2ind(size(ogrid), round((vects_w(1,:)-xmin)/res), round((vects_w(2,:)-ymin)/res))) + 10*inc;

    between = [];
    xs = round((s(1)-xmin)/res);
    ys = round((s(2)-ymin)/res);
    for j=1:100:numel(vects_w(3,:))-100;

        xf = (vects_w(1,j:j+99)-xmin)/res;
        yf = (vects_w(2,j:j+99)-ymin)/res;

        
        [x_between, y_between] = getMapCellsFromRay(repmat(xs, size(xf)),repmat(ys, size(yf)),double(xf),double(yf));
        x_between(x_between <= 0) = [];
        y_between(y_between <= 0) = [];
        %between = unique([x_between, y_between], 'rows');
        between = [x_between, y_between];
        ogrid(sub2ind(size(ogrid),between(:,1), between(:,2))) = ...
            ogrid(sub2ind(size(ogrid),between(:,1), between(:,2))) - inc;
    end
    %limit the occupancy grid
    ogrid(ogrid < -2) = -2;
    ogrid(ogrid > 2) = 2;
    [occupied_i, occupied_j] = find(siglim(ogrid) > occ_thres);
    
    %plotting
    set(im_hand, 'CDATA', 1-siglim(ogrid));
    %set(o, 'XDATA', (occupied_i*res)+xmin, 'YDATA', (occupied_j*res)+ymin);
    set(g, 'XDATA', round((vects_w(2,:)-ymin)/res), 'YDATA', round((vects_w(1,:)-xmin)/res));
    set(h, 'XDATA', round((Tr(2,:)-ymin)/res), 'YDATA', round((Tr(1,:)-xmin)/res));
    if(mod(i, 25) == 0)
        pause(0.025);
    end
    previous_pose = [pose(1:2), theta];
    s_hist = [s_hist, s];
end