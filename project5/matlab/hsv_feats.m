roads = Ihsv(:,:,1) <.35 & Ihsv(:,:,1) >.29 & Ihsv(:,:,2) <.11 & Ihsv(:,:,2) >.08 & Ihsv(:,:,3) < 0.45 & Ihsv(:,:,3) > 0.15;
roads = imdilate(roads, ones(9));
roads = imerode(roads, ones(5));
roads = imclose(roads, ones(7));
sidewalks = Ihsv(:,:,1) <.3 & Ihsv(:,:,1) >.14 & Ihsv(:,:,2) <.2 & Ihsv(:,:,2) >.1 & Ihsv(:,:,3) < 0.9 & Ihsv(:,:,3) > 0.6;
trees = Ihsv(:,:,1) <.35 & Ihsv(:,:,1) >.29 & Ihsv(:,:,2) <0.45 & Ihsv(:,:,2) >.25 & Ihsv(:,:,3) < 0.45 & Ihsv(:,:,3) > 0.25;
trees = imclose(trees, ones(15));
