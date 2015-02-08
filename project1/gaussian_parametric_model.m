function [ prob ] = phi(y, mu, sigma, d)
%estiamte of the probability of x given the gaussian specified by mu and
%sigma
%   Detailed explanation goes here

prob = 1/((2*pi)^(d/2) * sigma^.5) * exp(-0.5*(y-mu)'/(sigma)*(y-mu));

end

