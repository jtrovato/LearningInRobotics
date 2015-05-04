function [Tr, inliers2, x2] = process(I, replace, K)
%processes a new image in the VO pipeline
    
    %parameters
    maxfeats = 6;
    verbose = 0;
    bucket_width = 50;
    bucket_height = 50;
    world_height = 1.6;
    pitch = -0.08;
    motion_threshold = 100;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % perform feature detection, matching, and tracking
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tic; 
    % if we didnt use the last frame replace it with the new one and try again
    if replace
        matcherMex('replace', I);
    else
        matcherMex('push',I);
    end
    
    disp(['Feature detection: ' num2str(toc) ' seconds']);
    tic; 
    matcherMex('match',0); 
    matcherMex('bucketing', maxfeats, bucket_width, bucket_height); %bucket features
    p_matched = matcherMex('get_matches',0);
    i_matched = matcherMex('get_indices',0);
    disp(['Feature matching:  ' num2str(toc) ' seconds']);

    x1 = p_matched(3:4,:)';
    x2 = p_matched(1:2,:)';
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % use matches to determine camera pose
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [inliers1, inliers2, inds] = GetInliersRANSAC(x1, x2); 
    if sum(inds) < 10
        Tr = [];
        disp('not enough inliers');
        return;
    end
    F = EstimateFundamentalMatrix(inliers1 , inliers2);
    E = EssentialMatrixFromFundamentalMatrix(F, K);
    [tset, Rset] = ExtractCameraPose(E);
    
    Xset = cell(4,1);
    for i=1:4
        %start the transformation with the identity and build from there
        Xset{i} = LinearTriangulation(K, [0;0;0], eye(3), tset{i}, Rset{i}, inliers1, inliers2);
        verbose = 0;
        if verbose
            figure(3);
            points = Xset{i};
            Rplot = [0 0 1; -1 0 0; 0 -1 0];
            points_r = (Rplot*points')';
            subplot(2,2,i);
            showPointCloud(points_r(:,1), points_r(:,2), points_r(:,3));
        end
    end
    
    [t,R,X] = DisambiguateCameraPose(tset, Rset, Xset);
    % remove feature behind image plane (or rather keep only those in front
    X = X(X(:,3) > 0, :);
    if verbose
        figure(4)
        Rplot = [0 0 1; -1 0 0; 0 -1 0];
        points_r = (Rplot*points')';
        showPointCloud(points_r(:,1), points_r(:,2), points_r(:,3));
        ind
        pause
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Detect ground plane and get the scale factor
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % given 3D points find the distances from the origin
    l1_dist = sum(abs(X), 2);
    l2_dist = sum(X.^2, 2);
    median_dist = median(l1_dist);
    if median_dist > motion_threshold
        Tr = [];
        disp('does not pass motion threshold');
        return;
    end
    
    gnd_X = X(l1_dist<median_dist, :);
    if length(gnd_X) < 10
        Tr = [];
        disp('not enough ground points');
        return;
    end
    
    % find distance of camera from the ground plane
    h = getDistToGround(gnd_X, pitch, median_dist);
    t = t*world_height/h;
    
    Tr = [R , t; 0 0 0 1];

end

