% implementation of a Kalman Filter for ESE650
clear all; close all;

calibrateSensors

% Parameters: 
n = 6; %number of dimensions
h = []; % for plotting
x = zeros(n+1,length(imu_ts)+1);
P = zeros(n,n,length(imu_ts)+1);

% initialization:
% initialize to first orientation with very large covariance
w1 = [0, 0, 0]'; %calibrated_vals(4:6,1);
q = [1 0 0 0]';
x(:,1) = [q;w1];
P(:,:,1) = eye(6);
Q = eye(6); %process noise covariance (mean = 0)
R = eye(6); % measurement noise covariance (mean = 0)

for i = 2:length(imu_ts)
    tic;
    disp(['iteration', num2str(i)]);
    %%%%%%%%%%%%%%%  start of prediction step  %%%%%%%%%%%%%%%%%%%%%
    % calculate sigma points:
    %S = chol(P(:,:,i)+Q,'upper'); % Cholesky Decomposition
    [U,S,V] = svd(P(:,:,i-1)+Q);
    S = V*sqrt(S)*V';
    W = [-sqrt(2*n)*S, sqrt(2*n)*S]; %[6x12]
    X = zeros(n+1, 2*n);
    for j = 1:2*n
        Wq = rot2quat(W(1:3,j));
        x_q = quat_mult(x(1:4),Wq);
        x_omega = x(5:7,i)+ W(4:6,j); %[7x12] % rot to quat in W to add to state vector which is 7x1
        X(:,j) = [x_q;x_omega];
    end
    %Transform Sigma Points Xi -> Yi -> Zi (and calculate means and
    %variance alonf the way)
    
    % project sigma points ahead in time:
    % apply process model to X (sigma points) to get Y
    delta_t = imu_ts(i+1)-imu_ts(i);
    for j = 1:2*n
        omega = X(5:7,j);q = X(1:4,j);
        q_delta = omega2quatdelta(omega, delta_t);
        q_new = quat_mult(q,q_delta); % transpose everything because the function takes in row vectors
        omega_new = omega; %stays the same according to model
        Y(:,j) = [q_new;omega_new];
    end
    %Compute the mean of Y into W' (quaternions and angular vel are treated
    %differently) 
    % taking the principle eigenvector of the mean of the outter product of
    % the quaternions. (explained in class)
    Quats = Y(1:4,:);% just the quaternions
    q_mean = quatmean(Quats); %prior estimate
    omega_mean = mean(Y(5:7,:), 2);
    x_hat_bar = [q_mean; omega_mean]; % mean of sigma points Y
    
    % Transform sigma points Y [7x2n] to W_prime [6x2n]
    for j = 1:2*n
        r_W_prime = quat2rot(Y(1:4, j)); %convert quaternion to rotation vector
        omega_W_prime = (Y(5:7, j) - omega_mean);
        W_prime(:,j) = [r_W_prime; omega_W_prime]; %[6 x 2n]
    end
    P_bar = cov(W_prime');%(1/2*n)*(W_prime*W_prime'); % prior estimate
    
    %%%%%%%%%%%%%%%%%%%%%%  start of update step   %%%%%%%%%%%%%%%%%%%%%
    % Apply measurement model to Y (sigma points) to get Z (Z is for acc
    % and rot)
    g = [0 0 0 1]'; %vector quaternion of gravity
    for j=1:2*n
        q = Y(1:4,j);
        q_inv = [q(1); -q(2:4)]; % inv quat is in reverse direction, scaler stays the same.
        g_prime = quat_mult(quat_mult(q, g), q_inv);
        %z_acc = quat2rot(g_prime); 
        z_acc = g_prime(2:4);
        z_quat = Y(1:4, j); % Hrot = identity
        z_rot = quat2rot(z_quat); % should be our estimate of the rotational velocity
        Z(:,j) = [z_acc; Y(5:7, j)]; %estimate of our measurement vector [6x1] for each sigma point
    end
    
    %calculate mean  and covariance of Z:
    Z_bar = mean(Z,2);
    Pzz = cov(Z'); %[6x6]
    
    %calculate innovation and covariance
    measurement = calibrated_vals(:,i); %actual measurement
    measurement(4:6) = measurement(4:6)*(180/pi);
    v = measurement - Z_bar; %innovation
    Pvv = Pzz + R; %[6x6]
    
    %compute cross-correlation matrix
    Pxz = zeros(n);
    for j=1:2*n
        Pxz = Pxz + (1/2*n)*W_prime(:,j)*(Z(:,j)-Z_bar)'; %[6x6]
    end
    %compute Kalman Gain
    K = Pxz/Pvv; %[6x6]
    update = K*v;
    K_q = rot2quat(update(1:3));
    x(1:4,i+1) = quat_mult(x_hat_bar(1:4),K_q); %quaternion product
    x(5:7,i+1) = x_hat_bar(5:7) + update(4:6);
    P(:,:,i+1) = P_bar - K*Pvv*K';
    pause(delta_t - toc);
    
    verbose = 1;
    
    if verbose
        
        only_pred =0;
        if(only_pred)
            x(:,i+1) = x_hat_bar;
            P(:,:,i+1) = P_bar;
        end
        
        x(:,i+1)
        P(:,:,i+1)
        t = imu_ts(i);
        [min_val, ind] = min(vicon_ts - t);
        vicon_R = rots(:,:,ind);
        kf_R = quat_to_rot(x(1:4,i+1));
        h = newrotplot2(kf_R, vicon_R, h);
    end
end