% simple integration filter for ESE650
clear all; close all;

calibrateSensors

% Parameters:
x = zeros(4, length(calibrated_vals));
kf_R = zeros(3,3,length(calibrated_vals));
h = [];

% initialization:
q = [1 0 0 0]';
x(:,1) = q;

for i = 1:length(imu_ts)-1
    tic;
    disp(['iteration', num2str(i)]);
    %%%%%%%%%%%%%%%  start of prediction step  %%%%%%%%%%%%%%%%%%%%%
    %just apply motion model to orientation to get new orientation
    omega = calibrated_vals(4:6,i);
    delta_t = imu_ts(i+1)-imu_ts(i);
    q_delta = omega2quatdelta(omega, delta_t);
    q_new = quat_mult(x(:,i),q_delta); % transpose everything because the function takes in row vectors
    x(:,i+1) = q_new;
    kf_R(:,:,i) = quat_to_rot(x(:,i+1));

    
    verbose = 1;
    if verbose   
        x(:,i+1);
        t = imu_ts(i);
        [min_val, ind] = min(abs(vicon_ts - t));
        vicon_R = rots(:,:,ind);
        h = newrotplot2(kf_R(:,:,i), vicon_R, h);
        title(imu_ts(i)-imu_ts(1));
    end
    pause(delta_t - toc);
end