% train the HMMs
close all;

%%  load data (x6)
if ~exist('train_data')
    disp('loading data...');
    train_classes = [1,2,3,4,5,6];
    dcols = [2 3 4];
    train_inds= [1 2 3 4 5];
    train_dir = dir('./train/');
    [train_data, class_labels] = load_data(train_dir, train_classes, train_inds, dcols); %should be concatenated
    
end
%% Kmeans on data to get discrete state
if ~exist('all_data')
    all_data = [];
    for i=1:length(train_classes)
        all_data = [all_data; train_data(i).trials];
    end
%     %standardize the data
%     mean_d = mean(all_data);
%     std_d  = std(all_data);
%     all_data = (all_data - mean_d)/std_d;
%     for i=1:length(train_classes)
%         train_data(i).trials = (train_data(i).trials - mean_d)/std_d;
%     end
    
    K = 45;
    [labels, centroids] = kmeans(all_data, K); %kmeans on the acceleration
end

%% plot kmeans
color_set = varycolor(K);
figure();
hold on;
grid on;
title('KMeans Classification of Dataset');
xlabel('a_x'); ylabel('a_y'); zlabel('a_z');
for i=1:K
    plot3(all_data(labels == i,1), all_data(labels == i,2), all_data(labels == i,3),'.', 'Color', color_set(i,:)');
end
hold off;
%% for all the train
Ns = 8*ones(1,6);
pis = {};
As = {};
Bs = {};
for class=1:length(train_classes)
    disp(['------------------------------', class_labels{class}, '-----------------------------']);
    % discretize raw data into finite output sequence (x6) and label states
    disp('processing raw data');
    O = knnclassify(train_data(class).trials, centroids, [1:K]', 1);
    N = Ns(class);

    %% initialize A, B, Pi
    disp('initializing...');
    A = 0.9*eye(N) + 0.1*circshift(eye(N), -1); %the naive way

    B = 1/N*ones(K,N);

    Pi = zeros(1,N);
    Pi(1) = 1;

    % toy set
    % O = [1,1,2,2,3,3,4,4,1,2,3,4,1]';
    % noise = 0.2*randn(1, 12);
    % %initialize A, B, pi
    % A = [.5 .5 0 0
    %     0 .5 .5 0
    %     0 0 .5 .5
    %     .5 0 0 .5];
    % %A = circshift(eye(4), 1)';
    % B = .25*ones(4);
    % Pi = [1,0,0,0];

    %l = eval_model(Pi, A, B, O);
    %% Baum-Welch to learn model parameters (x6)
    disp('Baum Welch');
    max_iters = 100;
    thres = 1e-4;
    logprob_old = 0;
    for i=1:max_iters
        fprintf('iteration: %i', i);
        [Pi,A,B] = baumwelch(Pi,A,B,O);
        Pi;
        A;
        B;
        [state_seq, logprob] = MLEstateseq(Pi, A, B, O);
        fprintf('    logprob = %f \n', logprob);
        if i > 3 && ((logprob - logprob_old < thres*logprob_old) || logprob < logprob_old)
            fprintf('converged! \n');
            break;
        end
        logprob_old = logprob;
    end
    Pis{class} = Pi;
    As{class} = A;
    Bs{class} = B;
    
end


