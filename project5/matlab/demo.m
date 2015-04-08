%demo script

%% Load Image
I = imread('../aerial_color.jpg');
I = imresize(I, 0.25);
%% Generate Features (or load it)
verbose = 0;
if verbose
    F = generate_features(I);
else
    load('color_seg_poly.mat');
end
[~,d] = size(F);
%% Generate Training Paths (or load it)
verbose = 0;
if verbose
    imshow(I);
    pause
    train_paths = getTrainingPaths(5);
else
    load('paths.mat')
    num_paths = length(paths);
    for i=1:num_paths
        paths{i} = unique(round(paths{i}/2), 'rows', 'stable');
    end
end

%% Plot training paths
figure();
imshow(I);
hold on
for i=1:num_paths
    path = paths{i};
    plot(path(:,1), path(:,2), 'g', 'MarkerSize', 1);
    
end

%% Generate Cost Map (or load it)
a = ones(d,1)';
C = generateCostMap(a, F);
figure();
imagesc(reshape(C, size(I(:,:,1))));
%% Gradient Descent 
% Cost function = abs(C_desired - C_optimal)
% cost is defined as an exponential so derivative is easy
% parameter being optimized = weights on features
close all;

delta = 1e-3;
lambda = 0;
cost_hist = [];
learning_rate = 1e-4;
diff = 1;
weights = rand(1,d)-0.5;
%weights = (1/d)*ones(1,d);
weights_hist = weights;
[m,n] = size(I(:,:,1));

i=1;
while diff > delta
    tic;
    CMap = double(reshape(generateCostMap(weights, F), size(I(:,:,1))));
    %figure(3);
    %imshow(CMap);
    desired_grad = 0;
    desired_cost = 0;
    optimal_grad = 0;
    optimal_cost = 0;
    for p=1:num_paths
        fprintf('.')
        %the desired path
        path = round(paths{p});
        %generate Dijkstra optimal path using same start and end as desired
        xmax = min(max(path(:,1))+10, n);
        xmin = max(min(path(:,1))-10, 1);
        ymax = min(max(path(:,2))+10, m);
        ymin = max(min(path(:,2))-10, 1);
        cost_to_go = dijkstra_matrix(CMap(ymin:ymax,xmin:xmax),path(end,2)-ymin+1,path(end,1)-xmin+1);
        [yp, xp] = dijkstra_path(cost_to_go, CMap(ymin:ymax,xmin:xmax), path(1,2)-ymin+1, path(1,1)-xmin+1);
        xp_actual = yp+ymin-1;
        yp_actual = xp+xmin-1;
        dijk = [yp_actual, xp_actual];
        
        verbose = 0;
        if i > 0 && mod(p, 10) == 0
            figure(1);
            imshow(CMap(ymin:ymax,xmin:xmax))
            axis equal
            hold on
            plot(path(:,1)-xmin, path(:,2)-ymin, 'g', 'MarkerSize', 1);
            plot(xp, yp, 'r', 'MarkerSize', 1);
            title(['path ', num2str(p)]);
            hold off
            pause(0.025);
        end

        %calculate desired path "gradient" (for each feature)
        for j=1:length(path)
            desired_grad = desired_grad + CMap(path(j,2),path(j,1)).*F(path(j,2)*m +path(j,1),:);
            desired_cost = desired_cost + CMap(path(j,2),path(j,1));
        end
        %calculate optimal path "gradient" 
        for j=1:length(dijk)
            optimal_grad = optimal_grad + CMap(dijk(j,2),dijk(j,1)).*F(dijk(j,2)*m +dijk(j,1), :);
            optimal_cost = optimal_cost + CMap(dijk(j,2),dijk(j,1));
        end
    end
    %calculate cost
    cost = desired_cost - optimal_cost; %scaler
    grad = desired_grad - optimal_grad; % [dx1]
    

    %update the weights
    reg = weights/norm(weights);
    discount = 1;
    weights_new = weights -discount*learning_rate*(grad + lambda*reg) ;
    
    diff = sum(abs(weights_new-weights));
    weights = weights_new;
    weights_hist = [weights_hist; weights];
    cost_hist = [cost_hist; cost];
    
    figure(2);
    subplot(2,1,1);
    for l=1:d
        plot(weights_hist(:,l));
        hold on
    end
    hold off
    subplot(2,1,2);
    plot(cost_hist);
    pause(0.025);
    
    
    fprintf('iteration: %i,   cost: %f \n', i, cost);
    toc
    i = i+1;

end
    

%% Testing

