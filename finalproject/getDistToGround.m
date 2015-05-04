function h = getDistToGround(X, pitch, median)
%use a RANSAC like procedure to guess normal vectors to the image
%plane and subsequently determein distance of camera from the plane. 
sigma = median/50.0;
weight = .5*(sigma^2);


% project to 2D, because the ground plane only varies in y and z
X = X(:, 2:3);

% d is guesses of normal vector to plane.
n = [cos(-pitch), sin(-pitch)];
d = n*X';
N = length(d);

sum = 0;
best_sum = 0;
for i=1:N
    if d(i) > median/100
        for j=1:N
            dist = d(j)-d(i);
            sum = sum + exp(dist);
        end
        if sum > best_sum
            best_sum = sum;
            best_idx = i;
        end
    end
end

h = d(best_idx);

end

