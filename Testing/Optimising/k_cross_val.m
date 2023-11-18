close all
tic
load monkeydata_training.mat
k=5; % Num for k-fold cross validation
RMSE = zeros(k,1);

% Select training and testing data (you can choose to split your data in a different way if you wish)
trainPerc = 80;
samples=round(100/k); 

for rep = 1:k
    sampling_indx = false(1,100); % I figured logical indexing would come in handy for this
    sampling_indx((rep-1)*samples+1:rep*samples)=true;
    trainingData = trial(~sampling_indx,:);
%     rainingData = trial(ix(1:trainPerc),:);
    testData = trial(sampling_indx,:);
%     testData = trial(ix(trainPerc+1:end),:);
    modelParameters = positionEstimatorTraining(trainingData);
    fprintf('Testing the continuous position estimator...')

    meanSqError = 0;
    n_predictions = 0;  
    
    figure
    hold on
    axis square
    grid
    for tr=1:size(testData,1)
        display(['Decoding block ',num2str(tr),' out of ',num2str(size(testData,1))]);
        pause(0.001)
        for direc=randperm(8) 
            decodedHandPos = [];
    
            times=320:20:size(testData(tr,direc).spikes,2);
            
            for t=times
                past_current_trial.trialId = testData(tr,direc).trialId;
                past_current_trial.spikes = testData(tr,direc).spikes(:,1:t); 
                past_current_trial.decodedHandPos = decodedHandPos;
    
                past_current_trial.startHandPos = testData(tr,direc).handPos(1:2,1); 
                
                if nargout('positionEstimator') == 3
                    [decodedPosX, decodedPosY, newParameters] = positionEstimator(past_current_trial, modelParameters);
                    modelParameters = newParameters;
                elseif nargout('positionEstimator') == 2
                    [decodedPosX, decodedPosY] = positionEstimator(past_current_trial, modelParameters);
                end
                
                decodedPos = [decodedPosX; decodedPosY];
                decodedHandPos = [decodedHandPos decodedPos];
                
                meanSqError = meanSqError + norm(testData(tr,direc).handPos(1:2,t) - decodedPos)^2;
                
            end
            n_predictions = n_predictions+length(times);
            hold on
            plot(decodedHandPos(1,:),decodedHandPos(2,:), 'r');
            plot(testData(tr,direc).handPos(1,times),testData(tr,direc).handPos(2,times),'b')
        end
    end
    
    legend('Decoded Position', 'Actual Position')
    
    
    RMSE(rep) = sqrt(meanSqError/n_predictions)
    close all
    
    toc

end
