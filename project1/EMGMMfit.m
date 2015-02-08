%% EM for GMM fitting
% this algorithm is used to find the optimal parameters of the gaussians
% that make up out GMMs. 

max_iters = 10000;
[n_b,d] = size(y_b);
[n_o,~] = size(y_o);
k = 3; %number of gaussians in GMM
delta = 0.0001; %threshold to detemine convergence

%data
y_b = double(barrel_pixels);
y_o = double(other_pixels(1:8:end,:));

%initialization
rng(13);
w_b = 1/k*ones(1,k); % [1 x k]
w_o = 1/k*ones(1,k); % [1 x k]
mu_b = [[60, 120, 170]', [50, 120, 200]', [100, 140, 150]']; % [d x k]
mu_o = [[60,120,120]', [100,100,100]', [120,120,130]'];
%mu = ones(3);
SIGMA_b = reshape([255*rand*eye(d), 255*rand*eye(d), 255*rand*eye(d)], d,d,k); %stacked covariance matrices [d x d x k] matrices
SIGMA_o = reshape([255*rand*eye(d), 255*rand*eye(d), 255*rand*eye(d)], d,d,k); %stacked covariance matrices [d x d x k] matrices
%SIGMA = SIGMA/255
L_b = log_likelihood(y_b, w_b, mu_b, SIGMA_b); %scalar
L_o = log_likelihood(y_o, w_o, mu_o, SIGMA_o); %scalar


figure();
mu_image = reshape(mu_b', 1,k,d);
image(ycbcr2rgb(uint8(mu_image)));
drawnow;
for p=1:max_iters
    disp(['iteration:' num2str(p)]);
    %E-Step:
    gamma_b = assign_probs_to_data(y_b, w_b, mu_b, SIGMA_b); % predict the data labels (with probabilities)
    prob_sums_b = sum(gamma_b);
    gamma_o = assign_probs_to_data(y_o, w_o, mu_o, SIGMA_o); % predict the data labels (with probabilities)
    prob_sums_o = sum(gamma_o);
    %M-Step (update the model parameters):
    w_b = prob_sums_b./n_b; %new weights are accroding to probability mass assigned to each gaussian. Could do raw counts.
    w_o = prob_sums_o./n_o; %new weights are accroding to probability mass assigned to each gaussian. Could do raw counts.
    mu_b = (repmat((1./prob_sums_b), d, 1).*(gamma_b'*y_b)');
    mu_o = (repmat((1./prob_sums_o), d, 1).*(gamma_o'*y_o)');

    for j=1:k
        mean_centered_b = bsxfun(@minus, y_b, mu_b(:,j)');
        SIGMA_b(:,:,j) = 1/prob_sums_b(j)*bsxfun(@times, gamma_b(:,j), mean_centered_b)'*mean_centered_b;
        mean_centered_o = bsxfun(@minus, y_o, mu_o(:,j)');
        SIGMA_o(:,:,j) = 1/prob_sums_o(j)*bsxfun(@times, gamma_o(:,j), mean_centered_o)'*mean_centered_o;
    end
    mu_image = reshape(mu_b', 1,k,d);
    image(ycbcr2rgb(uint8(mu_image))); drawnow;
    
    %Check for convergence:
    Lnew_b = log_likelihood(y_b, w_b, mu_b, SIGMA_b);
    disp(['barrel log likelihood diff =', num2str(abs(Lnew_b-L_b))]);
    Lnew_o = log_likelihood(y_o, w_o, mu_o, SIGMA_o);
    disp(['other log likelihood diff =', num2str(abs(Lnew_o-L_o))]);
    if abs(Lnew_b-L_b) < delta && abs(Lnew_o-L_o) < delta
        disp('EM converged for other pixels')
        w_b
        mu_b
        SIGMA_b
        w_o
        mu_o
        SIGMA_o
        break;
    end
    L_b = Lnew_b;
    L_o = Lnew_o;
end