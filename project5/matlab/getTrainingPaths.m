function paths = getTrainingPaths(N)

paths = cell(N);
for i=1:N
    [cx,cy,~] = improfile;
    paths{i} = [cx, cy];
end

