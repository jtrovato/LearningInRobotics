function [ phi ] = gaussian_model(y, mu, SIGMA)
%estiamte of the probability of x given the gaussian specified by mu and
%sigma
%   phi is a [n x k] vector
[n,d] = size(y);
k = size(mu, 2);
phi = zeros(n,k);

for j = 1:k
    sigma = SIGMA(:,:,j);
    phi(:,j) = 1/((2.*pi)^(d/2) .* sqrt(det(sigma)))...
        .* exp(-0.5.*sum(((y-repmat(mu(:,j)', n,1))*...
        inv(sigma)).*(y-repmat(mu(:,j)', n,1)), 2));
end

end

