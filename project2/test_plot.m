figure();
h = [];
for i=1:length(rots)
R = rots(:,:,i);
h = newrotplot(R, h);
title(num2str(vicon_ts(i)-vicon_ts(1)));
end
