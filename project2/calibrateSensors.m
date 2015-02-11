% ESE650 Project 2
% author: Joe Trovato
% 
% This script loads training data nad calibrates the data to get meaningful
% values.
verbose = 0;


%% Load the data

imu_dir = dir('./ESE650 P2/imu');
total_imus = length(imu_dir)-2;
vicon_dir = dir('./ESE650 P2/vicon');
total_vicons = length(vicon_dir)-2;
percent_train = .01;
num_vicons = ceil(percent_train*10);
num_imus = ceil(percent_train*13);

for i=3:num_imus+2
    if i-2 <= num_vicons
        vicon = load(['./ESE650 P2/vicon/', vicon_dir(i).name]);
        vicon_ts = vicon.ts;
        rots = vicon.rots;
    end
    imu = load(['./ESE650 P2/imu/', imu_dir(i).name]);
    imu_ts = imu.ts;
    vals = imu.vals;
end

%% plot data
if verbose
    plot(imu_ts, vals(1,:), imu_ts, vals(2,:), imu_ts, vals(3,:));
    legend('x', 'y', 'z');
    xlabel('time');
    ylabel('acceleration');
    figure()
    plot(imu_ts, vals(4,:), imu_ts, vals(5,:), imu_ts, vals(6,:));
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
%% Calibrate the IMU data

% % this does not work well
% vref = 3300; %mV
% bit_res = 1023;
% sensitivitya = 2048;
% sensitivityg = 16.4;
% vals(1:3,:) = vals(1:3,:)*(vref/bit_res)/sensitivitya;
% vals(4:6,:) = vals(4:6,:)*(vref/bit_res)*(pi/180)/sensitivityg;


% This does:
ascale = 0.016; % g/bit
abias = [510;500;505];
acc = bsxfun(@minus,vals(1:3,:),abias)*ascale;

wscale = 1.5708e-4;
wbias = [370;373.6;375.5];
w = [vals(2:3,:);vals(1,:)];
w = bsxfun(@minus, -vals(4:6,:), wbias)*wscale;

calibrated_vals = [acc; w];
verbose =0;
if verbose
    figure()
    plot(imu_ts, acc(1,:), imu_ts, acc(2,:), imu_ts, acc(3,:));
    legend('x', 'y', 'z');
    xlabel('time');
    ylabel('g');

    figure()
    angles = cumsum(w')';
    plot(imu_ts, angles(1,:), imu_ts, angles(2,:), imu_ts, angles(3,:));
    legend('x', 'y', 'z');
    xlabel('time');
    ylabel('rad');
end


