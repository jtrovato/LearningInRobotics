function [ L ] = log_likelihood(y, w, mu, SIGMA)
%calcualtes the log likelihood that the data fits the model
%   L is a scalar
n = size(y,1);
L = (1/n).*sum(log(gaussian_model(y, mu, SIGMA)*w')); 

end

