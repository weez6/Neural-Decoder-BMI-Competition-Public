%% testing tSNE classification accuracy

load('tSNE_result.mat');
load('labelsCorrect.mat');

coords = Y;
labels = labelsCorrect;

ixx = randperm(length(Y));
train_set = Y(ixx(1:400),:);
test_set = Y(ixx(401:end),:);
train_labels = labels(ixx(1:400));
test_labels = labels(ixx(401:end));

correct = 0;
incorrect = 0;

for ix = 1:400
    pred_class(ix) = kNN(test_set(ix,:)', train_set', train_labels, 23);
    if pred_class(ix) == test_labels(ix)
        correct = correct + 1;
    else incorrect = incorrect + 1;
    end
end

correct

%% testing PCA/LDA for classification accuracy

template = repmat(1:8, 50, 1);
incorrect_tot = 0;
total = size(trimmed_labels,1) * size(trimmed_labels,2) * size(trimmed_labels,3);
for ix = 1:13
    diff = trimmed_labels(:,:,ix) - template;
    incorrect(ix) = size(find(diff),1);
    incorrect_tot = incorrect(ix) + incorrect_tot;
end
class_accuracy = (total - incorrect_tot)/ total




