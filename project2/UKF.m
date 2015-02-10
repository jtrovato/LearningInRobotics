% implementation of a Kalman Filter for ESE650
% Parameters: 
n = 6; %number of dimensions
x = zeroes(n,length(imu_ts));

% initialization:
% initialize to first orientation with very large covariance
w1 = vals(4:6,1);
q = [1 0 0 0]';
x(1) = [q;w1];
P = 10*eye(6);
Q = eye(6); %process noise covariance (mean = 0)
R = eye(6); % measurement noise covariance (mean = 0)

for i = 1:length(imu_ts)
    %%%%%%%%%%%%%%%  start of prediction step  %%%%%%%%%%%%%%%%%%%%%
    % calculate sigma points:
    S = chol(P+Q); % Cholesky Decomposition
    W = [-sqrt(2*n)*S, sqrt(2*n)*S]; %[6x12]
    X = bsxfun(@plus, x, [zeros(1,size(W,2));W]); %[7x12] % augmented the W vector so it could be added to steate vector which is 7x1
    
    %Transform Sigma Points Xi -> Yi -> Zi (and calculate means and
    %variance alonf the way)
    
    % project sigma points ahead in time:
    % apply process model to X (sigma points) to get Y
    delta_t = ts(i+1)-ts(i);
    for j = 1:2*n
        omega = X(5:7,j);q = X(1:4,j);
        q_delta = omega2quatdelta(omega, delta_t);
        q_new = quatmultiply(q',qdelta')'; % transpose everything because the function takes in row vectors
        omega_new = omega; %stays the same according to model
        Y(:,j) = [q_new;omega_new];
    end
    %Compute the mean of Y into W' (quaternions and angular vel are treated
    %differently) 
    % taking the principle eigenvector of the mean fo the outter product of
    % the quaternions. (explained in class)
    q_mean = quatmean(Y); %prior estimate
    omega_mean = mean(Y, 2);
    x_hat_bar = [q_mean; omega_mean]; % mean of sigma points Y
    
    % Transform sigma points Y [7x2n] to W_prime [6x2n]
    for j = 1:2*n
        r_W_prime = quat2rot(Y(1:4)); %convert quaternion to rotation vector
        omega_W_prime = Y(5:7) - omega_mean;
        W_prime = [r_W_prime; w_W_prime]; %[6 x 2n]
    end
    %calculate covariance in the process model prediction
    P_bar = (1/2*n)*sum(W_prime*W_prime', 2); % prior estimate
    
    %%%%%%%%%%%%%%%%%%%%%%  start of update step   %%%%%%%%%%%%%%%%%%%%%
    % Apply measurement model to Y (sigma points) to get Z (Z is for acc
    % and rot)
    g = [0 0 0 1]'; %vector quaternion of gravity
    for j=1:2*n
        q = Y(:,j);
        q_inv = [q(1); -q(2:4)]; % inv quat is in reverse direction, scaler stays the same.
        g_prime = quatmultiply(quatmultiply(q, g), q_inv);
        z_acc = g_prime; 
        z_quat = Y(1:4, :); % Hrot = identity
        z_rot = quat2rot(z_quat); % should be our estimate of the rotational velocity
        Z(:,j) = [z_acc, z_rot]; %estiamte of our measurement vector [6x1] for each sigma point
    end
    
    %calculate mean  and covariance of Z:
    Z_bar = mean(Z,2);
    Pzz = cov(Z'); %[6x6]
    
    %calculate innovation and covariance
    measurement = vals(:,i);
    v = measurement - Z_bar; %innovation
    Pvv = Pzz + R; %[6x6]
    
    %compute cross-correlation matrix
    Pxz = (1/2*n)*sum(W_prime*(bsxfun(@minus,Z, Z_bar)'), 2); %[6x6]
    %compute Kalman Gain
    K = Pxz/Pvv;
    x(i) = x_hat_bar + K*v;
    P = P_bar - K*Pvv*K';
    
end