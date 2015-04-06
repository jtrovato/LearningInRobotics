%demo script
verbose = 0;

%%
I = imread('../aerial_color.jpg');
I = imresize(I, 0.5);
%% Generate Features (or load it)
if verbose
    F = generate_features(I);
else
    load('feats.mat');
end
[~,~,d] = size(F);
%% Generate Training Paths (or load it)
verbose =1;
if verbose
    imshow(I);
    pause
    train_paths = getTrainingPaths(5);
else
    load('paths.mat')
end
num_paths = length(paths);
%% Generate Cost Map (or load it)
a = ones(d,1)';
c = generateCostMap(a, F);

%% Dijkstra Stuff


%% Gradient Descent 
% Cost function = abs(C_desired - C_optimal)
% cost is defined as an exponential so derivative is easy
% parameter being optimized = weights on features
delta = 1e-3;
learning_rate = 1;
diff = 1;
weights = (1/d)*ones(d,1);

i=0;
CMap = reshape(C,size(I));
while diff > delta
    desired_grad = 0;
    optimal_grad = 0;
    for p=1:num_paths
        %the desired path
        path = paths{i};
        %generate Dijkstra optimal path using same start and end as desired
        cost_to_go = dijkstra_matrix(CMap,path(end,1),path(end,2));
        [ip, jp] = dijkstra_path(ctg, costs, path(1,1), path(1,2));

        %calculate desired path "gradient" (for each feature)
        for j=1:length(path)
            desired_grad = desired_grad + C(path(j,1),path(j,2)).*F(path(j,1),path(j,2),:);
        end
        %calculate optimal path "gradient" 
        for j=1:length(ip)
            optimal_grad = optimal_grad + C(ip(j),jp(j)).*F(ip(j),jp(j),:);
        end
    end
    %calculate cost
    cost = desired_grad - optimal_grad;
    fprintf('iteration: %i,   cost: %f', i, cost);

    %update the weights
    weights_new = weights - (1/i)*learning_rate*cost;
    weights_new = weights_new/norm(weights_new);

    diff = sum(abs(weights_new-weights));
    weights = weights_new;
    i = i+1
end
    

%% Testing

