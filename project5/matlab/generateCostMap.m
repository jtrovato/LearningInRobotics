function C = generateCostMap(a, F)

cmin = 0.01; %no idea what to intialize this

% %first way
% d = length(a);
% C = reshape(F, [], d)';
% C = a*C;
% C = reshape(C, size(F(:,:,1)));

%second way
 C = bsxfun(@times, F, a);
 C = sum(C, 2);
 Cmean = mean(C,1);
 Cstd = std(C,0,1);
 C = bsxfun(@minus, C, Cmean);
 C = bsxfun(@rdivide, C,Cstd);

C = exp(C);
end