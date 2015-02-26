function [O] = rawtrain2obs(data, centroids)
% transform raw continuous sensor input to discrete observations
O = knnclassify(data, centroids, 1);
%maybe convert cluster to gaussain distributions

%maybe use GMM instead of kmeans

end

