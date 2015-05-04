% demonstrates monocular feature tracking (via feature indices)
disp('===========================');
clear all; dbstop error; close all;
addpath(genpath('labviso2'));

%video init
%vid = VideoWriter('monovo_progress');
%open(vid);

K = [7.070912000000e+02 0.000000000000e+00 6.018873000000e+02
    0.000000000000e+00 7.070912000000e+02 1.831104000000e+02
    0.000000000000e+00 0.000000000000e+00 1.000000000000e+00];


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
numims = 100;

% init transformation matrix array
Tr_total = eye(4);

%intial plot to set up handles
fig = figure(1);
subplot(2,1,1);
hi = imshow(I); hold on;
hp = plot(1, 1, 'xr');
hf = plot(1, 1, 'xg');
ha = plot(1, 1, '+g');
subplot(2,1,2);
%hpos = plot(Tr_total(1,4), Tr_total(3,4), 'b-o');
hpos = plot3(Tr_total(1,4), Tr_total(3,4), -Tr_total(2,4), 'b-o');
grid on;

replace = 0;
% start matching
for k=1:numims
    I = imread(['./08/mono_gray/' num2str(k,'%06d') '.png']);
    [Tr, inliers2, x2] = process(I, replace, K);
    % accumulate egomotion, starting with second frame
    if k>1

        % if motion estimate failed: set replace "current frame" to "yes"
        % this will cause the "last frame" in the ring buffer unchanged
        if isempty(Tr)
            replace = 1;
            Tr_total(:,:,k) = Tr_total(:,:,k-1);

        % on success: update total motion (=pose)
        else
            replace = 0;
            Tr_total(:,:,k) = Tr_total(:,:,k-1)*Tr; %inv becuase this is the motions from the second frame to the first?
        end
    end
    
    % update the plot with new matches
    set(hi, 'CDATA', I);
    %set(hp, 'XDATA', p_matched{end}(3, :), 'YDATA', p_matched{end}(4, :));
    set(hp, 'XDATA', inliers2(:, 1), 'YDATA', inliers2(:, 2));
    set(ha, 'XDATA', x2(:,1), 'YDATA', x2(:, 2));
    %if ~isempty(tracked_feats)
        %set(hf, 'XDATA', tracked_feats(1, :), 'YDATA', tracked_feats(2, :));
    %end
    set(hpos, 'XDATA', Tr_total(1,4,:), 'YDATA', Tr_total(3,4,:), 'ZDATA', -Tr_total(2,4,:));
    %set(hpos, 'XDATA', Tr_total(1,4,:), 'YDATA', Tr_total(3,4,:));
    pause(0.025);
    %frame = getframe(fig);
    %writeVideo(vid, frame);
    
end

% close matcher
matcherMex('close');
%close video
%close(vid);

