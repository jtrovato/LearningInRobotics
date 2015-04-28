% demonstrates monocular feature tracking (via feature indices)
disp('===========================');
clear all; dbstop error; close all;
addpath(genpath('labviso2'));

%video init
vid = VideoWriter('monovo_progress');
open(vid);

K = [9.842439e+02 0.000000e+00 6.900000e+02;
    0.000000e+00 9.808141e+02 2.331966e+02;
    0.000000e+00 0.000000e+00 1.000000e+00];

% matching parameters
param.nms_n                  = 2;   % non-max-suppression: min. distance between maxima (in pixels)
param.nms_tau                = 50;  % non-max-suppression: interest point peakiness threshold
param.match_binsize          = 50;  % matching bin width/height (affects efficiency only)
param.match_radius           = 200; % matching radius (du/dv in pixels)
param.match_disp_tolerance   = 1;   % du tolerance for stereo matches (in pixels)
param.outlier_disp_tolerance = 5;   % outlier removal: disparity tolerance (in pixels)
param.outlier_flow_tolerance = 5;   % outlier removal: flow tolerance (in pixels)
param.multi_stage            = 1;   % 0=disabled,1=multistage matching (denser and faster)
param.half_resolution        = 1;   % 0=disabled,1=match at half resolution, refine at full resolution
param.refinement             = 2;   % refinement (0=none,1=pixel,2=subpixel)

% init matcher
matcherMex('init',param);

% push back first image
I = imread('./08/mono_gray/000000.png');
matcherMex('push',I);

pos = [0 0 0 1]';
R = eye(3);

%intial plot to set up handles
fig = figure(1);
subplot(2,1,1);
hi = imshow(I); hold on;
hp = plot(1, 1, 'xr');
hf = plot(1, 1, 'xg');
subplot(2,1,2);
%hpos = plot3(pos(1), pos(2), pos(3), 'b-o');
hpos = plot(pos(1,:), pos(3,:), 'b-o');
grid on;


% feature tracks
tracks = {};
tracked_feats = [];
previous_feats = [];

% start matching
for im=1:100
    I = imread(['./08/mono_gray/' num2str(im*2,'%06d') '.png']);
    tic; matcherMex('push',I);
    disp(['Feature detection: ' num2str(toc) ' seconds']);
    tic; matcherMex('match',0);
    p_matched{im} = matcherMex('get_matches',0);
    i_matched{im} = matcherMex('get_indices',0);
    disp(['Feature matching:  ' num2str(toc) ' seconds']);

    x1 = p_matched{end}(3:4,:)';
    x2 = p_matched{end}(1:2,:)';
    [inliers1, inliers2, inds] = GetInliersRANSAC(x1, x2);
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
    
    [t,R,X, ind] = DisambiguateCameraPose(tset, Rset, Xset);
    if verbose
        figure(4)
        Rplot = [0 0 1; -1 0 0; 0 -1 0];
        points_r = (Rplot*points')';
        showPointCloud(points_r(:,1), points_r(:,2), points_r(:,3));
        ind
        pause
    end
    Tk = [eye(3), t; 0 0 0 1]; %[4x4]
    
    pos(:,end+1) = Tk*pos(:, end);
    pos(:,end) = pos(:,end)/pos(end,end);
    
    
    % update the plot with new matches
    set(hi, 'CDATA', I);
    %set(hp, 'XDATA', p_matched{end}(3, :), 'YDATA', p_matched{end}(4, :));
    set(hp, 'XDATA', inliers1(:, 1), 'YDATA', inliers1(:, 2));
    if ~isempty(tracked_feats)
        %set(hf, 'XDATA', tracked_feats(1, :), 'YDATA', tracked_feats(2, :));
    end
    %set(hpos, 'XDATA', pos(1,:), 'YDATA', pos(3,:), 'ZDATA', pos(3,:));
    set(hpos, 'XDATA', pos(1,:), 'YDATA', pos(3,:));
    pause(0.025);
    frame = getframe(fig);
    writeVideo(vid, frame);
    tracked_feats = [];
    
end

% close matcher
matcherMex('close');
%close video
close(vid);

