function [Pi_new, A_new, B_new]  = baumwelch(Pi, A, B, O)

N = size(A,1); %number of states
T = size(O, 1); %number of observations in time
K = size(B, 1); %number of discrete observations
eps = 1e-15;
%calculate forward probabilities

%initialization:
alpha = zeros(T, N);
za = zeros(T, 1);
a1 = Pi.*B(1,:); %first row of alpha (one row for each time)
norm_factor = sum(a1);
alpha(1,:) = a1./max(norm_factor, eps);
za(1) = log(norm_factor);

%induction:
for t=1:T-1
    a = (alpha(t,:)*A).*B(O(t+1),:);
    norm_factor = sum(a);
    alpha(t+1,:) = a./max(norm_factor, eps); % B(k,:) is the emission prob of seeing the obsevation k at all N states
    za(t+1) = za(t) + log(norm_factor);
end
alpha_real = bsxfun(@times, exp(za), alpha);
alpha;
%termination
L = sum(alpha(end,1), 2); 


%calclate backward probabilities

%intialize:
beta = zeros(T,N);
beta(T,:) = 1/N;
beta2 = zeros(T,N);
beta2(T,:) = 1/N;
zb = zeros(T,1);
zb(T) = log(N);
zb2 = zeros(T,1);
zb2(T) = log(N);
%induction:
for t=T-1:-1:1
    b2 = zeros(1,N);
    for i=1:N
       summ = 0;
       for j =1:N
           summ = summ + A(i,j)*B(O(t+1),j)*beta(t+1,j);
       end
       b2(i) = summ;
    end
    
    b = (B(O(t+1),:).*beta(t+1,:))*A'; % B(k,:) is the emission prob of seeing the obsevation k at all N states
    norm_factor = sum(b);
    norm_factor2 = sum(b2);

    beta(t, :) = b./max(norm_factor, eps);
    beta2(t, :) = b2./max(norm_factor2, eps);
    zb(t) = zb(t+1) + log(norm_factor);
    zb2(t) = zb2(t+1) + log(norm_factor2);

end
beta_real = bsxfun(@times, exp(zb), beta);
beta;
beta2;

gamma_2 = zeros(T,N);
for t=1:T
    gamma_2(t,:) = bsxfun(@rdivide, alpha(t,:).*beta(t,:), sum(alpha(t,:).*beta(t,:), 2));
end
gamma_2;
% [row_max, SPred] = max( gamma_2, [], 2 );
% SPred = SPred';


% calcualte xi
xi = zeros(N,N,T-1);
for t=1:T-1
    num = repmat(alpha(t,:)',1,N).*A.*repmat(B(O(t+1),:).*beta(t+1,:),N,1);
    xi(:,:,t) = num/max(sum(sum(num)), eps);
end
xi;
% calculate gamma
gamma = zeros(T,N);
for t=1:T-1
    gamma(t,:) = reshape(sum(xi(:,:,t),2),1,N);
    
end
gamma = gamma_2; %to use the other calculation for gamma

%M-Step: calcualte new model parameters
% statistics over gamma and xi
exp_num_trans_from_i = sum(gamma(1:end-1,:), 1); %Nx1
exp_num_trans_from_i_to_j = sum(xi, 3); %NxN

Pi_new = gamma(1,:);
A_new = exp_num_trans_from_i_to_j./repmat(max(exp_num_trans_from_i, eps)', 1, N);
B_new = zeros(K,N);
for k=1:K
    inds_where_o_is_k = (O == k);
    B_new(k, :) = sum(gamma(inds_where_o_is_k, :))./max(sum(gamma), eps);
end
B_new = bsxfun(@rdivide, (B_new + eps), sum(B_new));