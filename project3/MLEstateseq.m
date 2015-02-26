function [state_seq, logprob] = MLEstateseq(Pi, A, B, O)
%computes most probable state sequence and logprob
N = length(A);
T = length(O);
state_seq = zeros(T,1);
phi = zeros(T, N);
psi = zeros(T, N);
phi(1, :) = log(Pi) + log(B(O(1), :));
for t = 2:T
    [max_val, max_inds] = max(bsxfun(@plus, phi(t-1, :)', log(A)));
    phi(t, :) =  max_val + log(B(O(t), :));
    psi(t, :) = max_inds;
end
[~, state_seq(T)] = max(phi(T, :));
for t = T-1:-1:1
    state_seq(t) = psi(t+1, state_seq(t+1));
end

logprob = max(phi(T,:));

