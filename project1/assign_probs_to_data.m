function [ gamma ] = assign_probs_to_data(y, w, mu, SIGMA)
%This function guesses the probability that data point i was generated by
%gaussian j
%   gamma should be a [n x k] vector
%   the sum of a row is 1 because it weights each probability so the sum is
%   1
[n,d] = size(y);
k = size(mu, 2);
% the w may or may not be necessary
gamma = repmat(w, n, 1) .* gaussian_model(y, mu, SIGMA);
norm_factor = sum(gamma, 2);
gamma = gamma ./ repmat(norm_factor, 1, k);
end

