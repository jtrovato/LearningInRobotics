function [STATS, RectMetric, AxisMetric ] = rejectCCs(STATS)
%reject connected components that dont look like barrels

bad_inds = [];
for ii=1:length(STATS)
    if STATS(ii).ConvexArea < 3000
        bad_inds = [bad_inds, ii];
    end
    if STATS(ii).Extent < .7
        bad_inds = [bad_inds, ii];
    end
    %if STATS(ii).Orientation < 75 || STATS(ii).Orientation > 115
    %    bad_inds = [bad_inds, ii];
    %end
    ideal_ratio = 40/57;
    axis_ratio = STATS(ii).MinorAxisLength/STATS(ii).MajorAxisLength;
    AxisMetric(ii) = abs(axis_ratio-ideal_ratio)/ideal_ratio;
    RectMetric(ii) = STATS(ii).ConvexArea - prod(STATS(ii).BoundingBox(3:4))/prod(STATS(ii).BoundingBox(3:4));
    if RectMetric < .7
        bad_inds = [bad_inds, ii];
    end
end
STATS(bad_inds) = [];
AxisMetric(bad_inds) = [];
RectMetric(bad_inds) = [];

end

