% test script
warning('off')

vs = 1:5;
train_inds = [];
for v=1:10
    fprintf('fold = %i \n', v);
    %train HMMs
    num_models = 6;
    test_inds = round(5*rand());
    train_inds = vs(vs ~= test_inds);
    trainHMMs

    %load test data
    test_files = dir();
    test_classes = [1 2 3 4 5 6];
    %test_inds= [1];
    dcols = [2 3 4];
    test_dir = dir('./train/');
    [test_data, class_labels] = load_test_data(test_dir, test_classes, test_inds, dcols); %should be concatenated

    figure();
    for c=test_classes
        data_cell = test_data(c).trials;
        for trial=test_inds
            data = data_cell{trial};
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
            fprintf('predicted: %s | actual: %s | confidence: %f | rank: %s, %s, %s, %s, %s, %s \n',...
                prediction, class_labels{c}, confidence, class_labels{rank(1)}, class_labels{rank(2)},...
                class_labels{rank(3)}, class_labels{rank(4)}, class_labels{rank(5)}, class_labels{rank(6)});
        end
    end
end