function X = LinearTriangulation(K, C1, R1, C2, R2, points1, points2)
%triangulates points from two images in 3D space
P1 = K*[R1 -C1];
P2 = K*[R2, -C2];
vec2skew = @(v) [0, -v(3), v(2); v(3), 0, -v(1); -v(2), v(1), 0];
numpts = length(points1);

X = zeros(numpts, 3);
for i=1:numpts
    skew1 = vec2skew(points1(i,:)); 
    skew2 = vec2skew(points2(i,:));

    A = [skew1*P1; skew2*P2];
    %A = normx(A')';
    [~,~,V] = svd(A, 0);
    xrow = V(:,end)/V(end,end);
    X(i, :) = xrow(1:3)';
    %X(i,:) = V(1:3,end);
end

% x = normx(x)  Normalize MxN matrix so that norm of each its column is 1.
function x = normx(x)
if ~isempty(x)
  x = x./(ones(size(x,1),1)*sqrt(sum(x.*x)));
end
end

end

