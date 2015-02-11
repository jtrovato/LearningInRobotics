figure();
h = [];
for i=1:length(rots)
R = rots(:,:,i);
h = newrotplot(R, h);
i
end
