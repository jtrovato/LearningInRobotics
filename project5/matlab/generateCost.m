function C = generateCost(a,F,x,y)

cmin = 0.001; %no idea what to intialize this

% %first way
% F = F(x,y,:);
% C = reshape(F, [], d)';
% C = a*C;
% C = reshape(C, size(F(:,:,1)));

%second way
C = bsxfun(@times, a, F);
C = sum(C, 3);

C = cmin + exp(C);

end

