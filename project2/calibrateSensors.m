% ESE650 Project 2
% author: Joe Trovato
% 
% This script loads training data and calibrates the data to get meaningful
% values.
verbose = 0;


%% Load the data



%     imu_dir = dir('./ESE650 P2/imu');
%     vicon_dir = dir('./ESE650 P2/vicon');
%     cam_dir = dir('./ESE650 P2/cam');
    vicon_filename = 'viconRot2.mat';
    cam_filename = 'cam2';
    imu_filename = 'imuRaw2';
    
if ~test
    vicon = load(['./ESE650 P2/vicon/',vicon_filename]);
    vicon_ts = vicon.ts;
    rots = vicon.rots;
    imu = load(['./ESE650 P2/imu/', imu_filename]);
    imu_ts = imu.ts;
    vals = imu.vals;
    cam = load(['./ESE650 P2/cam/', cam_filename]);
    cam_ts = cam.ts;
    cam = cam.cam;
else
    imu = load('./P2_TEST/imu_test.mat');
    imu_ts = imu.ts;
    vals = imu.vals;
    cam = load('./P2_TEST/cam_test.mat');
    cam_ts = cam.ts;
    cam = cam.cam;
end


%% plot data
if verbose
    plot(imu_ts, vals(1,:), imu_ts, vals(2,:), imu_ts, vals(3,:));
    legend('x', 'y', 'z');
    xlabel('time');
    ylabel('acceleration');
    figure()
    plot(imu_ts, vals(5,:), imu_ts, vals(6,:), imu_ts, vals(4,:));
    legend('z', 'x', 'y');
    xlabel('time');
    ylabel('gyro value');

    %% plot the rotations
    R = rots(:,:,1); 
    rotplot;
    for i = 2:length(rots)
        tic;
        deltat = vicon_ts(i+1) - vicon_ts(i);
        xp = rots(:,:,i)*x;
        %set(p3, 'Xdata',xp(1,itop), 'YData',xp(2,itop),'ZData',xp(3,itop));
        set(p2, 'Xdata',xp(1,ibottom), 'YData',xp(2,ibottom),'ZData',xp(3,ibottom));
        set(s1, 'Xdata',xp(1,ifront), 'YData',xp(2,ifront),'ZData',xp(3,ifront));
        set(s2, 'Xdata',xp(1,iback), 'YData',xp(2,iback),'ZData',xp(3,iback));
        drawnow
        %while(toc < deltat)
        %end
    end
end
if ~test
    %% plot the vicon angular velocity
    Rd = zeros(3,3,length(rots)-1);
    for i = 2:length(rots)
        Rd = rots(:,:,i-1)'*rots(:,:,i);
        r = vrrotmat2vec(Rd);
    end
end

%% Calibrate the IMU data

% % this does not work well
% vref = 3300; %mV
% bit_res = 1023;
% sensitivitya = 2048;
% sensitivityg = 16.4;
% vals(1:3,:) = vals(1:3,:)*(vref/bit_res)/sensitivitya;
% vals(4:6,:) = vals(4:6,:)*(vref/bit_res)*(pi/180)/sensitivityg;


% This does:
ascale = [-0.097;-0.097;0.097]; %/9.8; % g/bit %make x and y negative!!!
abias = [510;500;505];
acc = bsxfun(@times,bsxfun(@minus,vals(1:3,:),abias), ascale);

wscale =  0.0158; %
wbias = [373.567;375.3606;369.7041];
w = [vals(5:6,:);vals(4,:)];
w = bsxfun(@minus, w, wbias)*wscale;

calibrated_vals = [acc; w];
verbose =0;
if verbose
    figure()
    plot(imu_ts, acc(1,:), imu_ts, acc(2,:), imu_ts, acc(3,:));
    legend('x', 'y', 'z');
    xlabel('time');
    ylabel('g');

    figure()
    plot(imu_ts, w(1,:), imu_ts, w(2,:), imu_ts, w(3,:));
    legend('x', 'y', 'z');
    xlabel('time');
    ylabel('rad/s');
end


