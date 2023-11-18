function [x,y,modelParameters]= positionEstimatorTest(past_current_trial, modelParameters, k)
    endTimeTrain = 560;
    bin_size = modelParameters.bin_size;
    endTimeTrial = size(past_current_trial.spikes,2);
    binnedData =  bin_firing(past_current_trial, bin_size,endTimeTrial,1); % preprocessing
    
%        tic     
    % get the relevant parameters from the model
    if endTimeTrial <= endTimeTrain
        
        indx =  (endTimeTrial/bin_size)-(320/bin_size)+1;
        inactNeurons = modelParameters.inactNeurons{1};
        binnedData(inactNeurons,:,:) = [];

        trainLDA = modelParameters.DimReduction(indx).wLDA;
        LDAout = modelParameters.DimReduction(indx).LDAout;
        meanFTrain = modelParameters.DimReduction(indx).meanFiring;
        
%         firingData = mean(binnedData,1)';
%         firingData = reshape(binnedData',[],1);
        firingData = binnedData(:);
       
        testLDA = LDAout'*(firingData -meanFTrain); 
       
       predLabel = kNN(testLDA, trainLDA, k);


%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%SVM%%%%%%%%%%%%%%%%%%%%%
%         tr_lab = repelem([1:8], 50);
%         bayes_model = fitcnb(trainLDA',tr_lab);
%         predLabel = predict(bayes_model, testLDA');
%         
% 
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55

        modelParameters.Label = predLabel;

        xPath = modelParameters.AvPaths.meanX(predLabel,:);
        yPath = modelParameters.AvPaths.meanY(predLabel,:); 

        x =  xPath(endTimeTrial);
        y = yPath(endTimeTrial);

    else 
        % no classification is done, we assume that the direction of the 
        % movement is the same as in previous time indexes

        predLabel = modelParameters.Label;

        xPath = modelParameters.AvPaths.meanX(predLabel,:);
        yPath = modelParameters.AvPaths.meanY(predLabel,:); 

        if endTimeTrial<=size(xPath,2)
            x =  xPath(endTimeTrial);
            y= yPath(endTimeTrial);
        else
            x =  xPath(end);
            y = yPath(end);
        end

    end
%     toc
end     
    
               
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [class] = kNN(test_point, train_data, k)
%alg from: https://www.tutorialspoint.com/machine_learning_with_python/machine_learning_with_python_knn_algorithm_finding_nearest_neighbors.htm
num_points = size(train_data,2)/8;
train_labels = repelem([1:8], num_points);
% Finding euclidean distance between the testing point and each of the
% training data in 4D space

% diff = train_data - test_point;
% distances = diag(diff'*diff);

distances = ones(1, size(train_data,1));
for i = 1:size(train_data, 2)
    distances(i) = sum((test_point - train_data(:,i)) .^ 2);
end  



   % Sort distances and labels by distance
   [~, idx] = sort(distances);
   
   sorted_labels = train_labels(idx);
   
   % Find k nearest neighbors
   
   k_nearest_labels = sorted_labels(1:k);
   class = mode(k_nearest_labels,"all");
   

end


function bin_fir = bin_firing(data, bin_size,endTime,startTime)

% bin_fir=cell(size(data,1),size(data,2));
    dirNum=size(data,2);
    trialNum=size(data,1);
    binNum = ceil((endTime-startTime)/bin_size);
    neuralNum = size(data(1,1).spikes,1);
    bin_sum_arr = zeros(neuralNum,binNum);
%     windowSize = 10; 
%     numer= windowSize:-1:1;
%     b = (1/sum(numer))*numer;
%     a = 1;
    
    bin_fir=zeros(neuralNum,binNum,dirNum*trialNum);
    
    for column = 1:dirNum %going through all 8 movements/degrees
        for row =1:trialNum %going through 100 trials of each movement
            tempData = data(row,column).spikes;
            L=size(tempData,2);
            %binFiring = cell(size(trial,1),ceil)
            
%             for m = 1:neuralNum %looping through each of the 98 rows
                a=1;
                for n =startTime:bin_size:endTime
                    if (n + bin_size) > L % if index + bin_size exceeds data array length
                        bin_array = tempData(1:neuralNum,(n:L)); 
                        bin_sum = sum(bin_array,2); % count the number of spikes
                        bin_sum_arr(1:neuralNum,a)=bin_sum; % add result to the matrix of bin firings
                    else 
                        bin_array = tempData(1:neuralNum,n:(n+bin_size));
                        bin_sum = sum(bin_array,2); % count the number of spikes
                        bin_sum_arr(1:neuralNum,a)=bin_sum; % add result to the matrix of bin firings
                    end
                   
                    
                    a=a+1; % update bin counter (column index for bin_sum_arr)
                end
%             end
%             bin_sum_arr= filter(b,a,bin_sum_arr,[],2);
            bin_fir(:,:,(column-1)*trialNum+row)=bin_sum_arr; % return firings of bins for each trial of each movement
        end
    end

   
end