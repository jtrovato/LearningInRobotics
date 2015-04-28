% demonstrates monocular feature tracking (via feature indices)
disp('===========================');
clear all; dbstop error; close all;
addpath(genpath('labviso2'));

% matching parameters
param.nms_n                  = 6;   % non-max-suppression: min. distance between maxima (in pixels)
param.nms_tau                = 80;  % non-max-suppression: interest point peakiness threshold
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

%intial plot to set up handles
hi = imshow(I); hold on;
hp = plot(1, 1, 'xr');
hf = plot(1, 1, 'xg');

% feature tracks
tracks = {};
tracked_feats = [];

% start matching
for im=1:200
    I = imread(['./08/mono_gray/' num2str(im,'%06d') '.png']);
    tic; matcherMex('push',I);
    disp(['Feature detection: ' num2str(toc) ' seconds']);
    tic; matcherMex('match',0);
    p_matched{im} = matcherMex('get_matches',0);
    i_matched{im} = matcherMex('get_indices',0);
    disp(['Feature matching:  ' num2str(toc) ' seconds']);

    % for all matches in the last frame
    for i=1:size(i_matched{end},2)

        ind = i;

        % init track with latest matches
        p1 = p_matched{end}(3:4,ind);
        p2 = p_matched{end}(1:2,ind);
        p  = [p1 p2];

        % augment track into the past
        for j=length(p_matched)-1:-1:1

            % find backwards
            ind = find(i_matched{j}(2,:)==i_matched{j+1}(1,ind));
            if isempty(ind)
              break;
            end

            p3 = p_matched{j}(1:2,ind);
            p  = [p p3];
        end

        track_length = length(p_matched)-j;
        if track_length>8
            tracked_feats = [tracked_feats, p(:, 1)];
        end
    end
    
    % update the plot with new matches
    set(hi, 'CDATA', I);
    set(hp, 'XDATA', p_matched{end}(3, :), 'YDATA', p_matched{end}(4, :));
    if ~isempty(tracked_feats)
        set(hf, 'XDATA', tracked_feats(1, :), 'YDATA', tracked_feats(2, :));
    end
    pause(0.025);
    tracked_feats = [];
    
end

% close matcher
matcherMex('close');


