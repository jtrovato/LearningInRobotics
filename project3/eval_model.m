function likelihood = eval_model(Pi, A, B, O)
% compute the likelihood of a sequence of observations (O) given a HMM
% model (Pi,A,B) using the forward backward algorithm (really jsut the
% forward part)

N = size(A,1); %number of states
T = size(O, 1); %number of observations in time
%initialization:
alpha = zeros(T, N);
za = zeros(T, 1);
alpha(1,:) = Pi.*B(1,:); %first row of alpha (one row for each time)

%induction:
for t=1:T-1
    a = (alpha(t,:)*A).*B(O(t),:);
    norm_factor = sum(a);
    alpha(t+1,:) = a./norm_factor; % B(k,:) is the emission prob of seeing the obsevation k at all N states
    za(t+1) = za(t) + log(norm_factor);
end
alpha_real = bsxfun(@times, exp(za), alpha);

%termination
likelihood = sum(alpha(end,1), 2); 


end