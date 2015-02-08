function [ prob ] = gaussian_model(y, mu, sigma)
%estiamte of the probability of x given the gaussian specified by mu and
%sigma
%   
d = size(y, 2);
prob = 1/((2*pi)^(d/2) * sigma^.5) * exp(-0.5*(y-mu)'/(sigma)*(y-mu));

end

