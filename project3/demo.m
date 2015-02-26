% test script
warning('off')
train = 0;
test_dir_name = './test_sample/';

if train == 1
    %train HMMs
    num_models = 6;
    trainHMMs
else
    load('models.mat');
end


figure();
num_models = 6;
dcols = [2 3 4];
test_dir = dir(test_dir_name);
for j=3:length(test_dir)
    
    data = importdata([test_dir_name, test_dir(j).name]);
    data = data(:,dcols);
    %generate observation sequence for different 
    O = knnclassify(data, centroids, [1:K]', 1);

    %maximum likelihood over models
    likelihoods = zeros(num_models, 1);
    for i=1:num_models
        [state_seq, logprob] = MLEstateseq(Pis{i},As{i},Bs{i},O);
        likelihoods(i) = logprob;
        subplot(num_models,1,i);
        plot(state_seq);
    end
    drawnow;
    likelihoods;
    [val, label] = max(likelihoods);
    [confidence, rank] = eval_pred(likelihoods);
    prediction = class_labels{label};
    probability = val;
    fprintf('predicted: %s | confidence: %f | rank: %s, %s, %s, %s, %s, %s \n',...
        prediction, confidence, class_labels{rank(1)}, class_labels{rank(2)},...
        class_labels{rank(3)}, class_labels{rank(4)}, class_labels{rank(5)}, class_labels{rank(6)});
end

