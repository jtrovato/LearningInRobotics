% from binary images output from color segmentation find barrels in image
close all;

%% Find barrel

CC = bwconncomp(red_mask);
STATS1 = regionprops(CC, 'BoundingBox', 'ConvexArea', 'ConvexHull', 'Centroid',...
    'MajorAxisLength', 'MinorAxisLength', 'Extent', 'Orientation');
[STATS1, RectMetric1, AxisMetric1] = rejectCCs(STATS1);

if isempty(STATS1)
    SE = strel('square',25);
    red_mask = imdilate(red_mask, SE);
    CC = bwconncomp(red_mask);
    STATS2 = regionprops(CC, 'BoundingBox', 'ConvexArea', 'ConvexHull', 'Centroid',...
    'MajorAxisLength', 'MinorAxisLength', 'Extent', 'Orientation');
    [STATS2, RectMetric2, AxisMetric2] = rejectCCs(STATS2);

    SE = strel('square',20);
    red_mask = imdilate(red_mask, SE);
    CC = bwconncomp(red_mask);
    STATS3 = regionprops(CC, 'BoundingBox', 'ConvexArea', 'ConvexHull', 'Centroid',...
    'MajorAxisLength', 'MinorAxisLength', 'Extent', 'Orientation');
    [STATS3, RectMetric3, AxisMetric3] = rejectCCs(STATS3);
    STATS = [STATS2; STATS3];
    AxisMetric = [AxisMetric2, AxisMetric3];
    RectMetric = [RectMetric2, RectMetric3];
else
    STATS = STATS1;
    AxisMetric = AxisMetric1;
    RectMetric = RectMetric1;
    
end

if ~isempty(STATS)
    [m, ind] = max(RectMetric);

    figure();
    verbose = 0;
    if(verbose)
        subplot(2,1,1);
        imshow(red_mask);
        subplot(2,1,2);
        imshow(test_im_rgb);
    else
        imshow(test_im_rgb)
    end
    hold on;
    for i=1:length(STATS)
        plot(STATS(i).ConvexHull(:,1), STATS(i).ConvexHull(:,2), 'w', 'LineWidth', 1);
    end
    plot(STATS(ind).ConvexHull(:,1), STATS(ind).ConvexHull(:,2), 'g', 'LineWidth', 1);
    

    plot(y, x, 'g+');
    title(imdir(cur_im).name);

    %% Projection Geometery to find barrel distance
    barrelh = 0.57; %meters
    barrelw = 0.40; %meters

    sensorh = 0.0036;
    sensorw = 0.0048;
    f = 0.00367; %meters

    imageh = 1200;
    imagew = 1600;

    p2x = sensorh/imageh;
    p2y = sensorw/imagew;

    x_pix = 0.85*STATS(ind).MinorAxisLength; %scaling the axis because the convex hullis often wider than the actual barrel
    y_pix = 0.85*STATS(ind).MajorAxisLength;
    x = x_pix*p2x;
    y = y_pix*p2y;

    zx = (barrelw*f)/x;
    zy = (barrelh*f)/y;
    %use function fit from data to estimate distance
    x = STATS(ind).Centroid(2);
    y = STATS(ind).Centroid(1);
    d = round(barrel_dist_func(x_pix), 1);
    text(30,50, [num2str(d), 'm'], 'Color', [0,1,0], 'FontSize', 15);
else
    figure(101);
    subplot(2,1,1);
    imshow(red_mask);
    subplot(2,1,2);
    imshow(test_im_rgb);
    disp('could not find barrel');
end


