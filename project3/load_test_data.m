function [data, class_labels] = load_test_data(train_dir, train_classes, train_inds, dcols)
    %read data from imu
    
    data = [];
    class_labels = {};
    %concatenate data for each gesture
    for i=3:length(train_dir)
        classl = train_dir(i).name;
        if max(i-2 == train_classes)
            class_labels{i-2} = classl;
            data(i-2).('label') = classl;
            class_dir = dir(['./train/', classl]);
            trial_data = {};
            for j=train_inds+2
                d = importdata(['./train/', classl, '/', class_dir(j).name]);
                trial_data{j-2} = d(:, dcols);
            end
            data(i-2).trials = trial_data;
        end
    end
end
