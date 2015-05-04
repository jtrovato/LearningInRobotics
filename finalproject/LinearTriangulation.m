function X = LinearTriangulation(K, C1, R1, C2, R2, x1, x2)
%triangulates points from two images in 3D space
if size(x1,2) == 2 || size(x2,2)==2
    x1 = [x1, ones(size(x1,1),1)];
    x2 = [x2, ones(size(x2,1),1)];
end

P1 = K*R1*[eye(3) -C1];
P2 = K*R2*[eye(3), -C2];
vec2skew = @(v) [0, -v(3), v(2); v(3), 0, -v(1); -v(2), v(1), 0];
numpts = length(x1);

X = zeros(numpts, 3);
for i=1:numpts
    skew1 = vec2skew(x1(i,:)); 
    skew2 = vec2skew(x2(i,:));

    A = [skew1*P1; skew2*P2];
    [~,~,V] = svd(A);
    X(i, :) = V(1:3,end)'/V(end,end); %normalize points
 
end
