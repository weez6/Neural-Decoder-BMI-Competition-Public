function  [modelParameters,firingData] = positionEstimatorTraining(trainingData, pca_dim, lda_dim, bin_size, firing_threshold)

    directionNum = 8;
    pcaNum = pca_dim;
    ldaNum = lda_dim;

    %bin_size = 80;
    startTime = 320; % onset of movement
    endTimeTrain = 560; % size of smallest training set

    trainPerc =  length(trainingData); % how many rows of the initial data we are using
    binnedData =  bin_firing(trainingData, bin_size, endTimeTrain,1); % preprocessing
    
    modelParameters = struct;
    modelParameters.bin_size = bin_size;
    
    neuralNum = size(binnedData,1);
    binIndeces = (startTime:bin_size:endTimeTrain)/bin_size;
    
   
    %% Obtain the average path for each direction
     [xPath,yPath] = getAvPath(trainingData,directionNum,trainPerc);     
     modelParameters.AvPaths.meanX = xPath;
     modelParameters.AvPaths.meanY = yPath;

    %% Detect inactive neural units
    inactNeuron = [];
    for a = 1: neuralNum
        firings = mean(binnedData(a,:,:),"all");
        if firings < firing_threshold
            inactNeuron = [inactNeuron a];
        end
    end
    inactiveNeurons_cell = cell(1,1);
    inactNum = length(inactNeuron);
    inactiveNeurons_cell{1,1} = inactNeuron;
    modelParameters.inactNeurons = inactiveNeurons_cell;

    % supervised labelling for Linear Discrminant Analysis
    dirLabels = repelem(1:8,trainPerc);        
%         dirLabels = repmat(1:8,1,trainPerc);

    binNum = 1;
    for testingIndxs = binIndeces
        % flatten data: rows:(n0 of bins from 1 to current time) columns:trial
        firingData = binnedData(:,1:testingIndxs,:);
        firingData(inactNeuron,:,:) = [];

        firingData = reshape(firingData,[],directionNum*trainPerc);
%         firingData = reshape(firingData,[neuralNum*testingIndxs,directionNum*trainPerc,1]);
%         firingData=squeeze(fft(firingData,1));
    
        %% Do Principal Component Analysis 
%         [princComp, ~, ~] = svd(firingData - mean(firingData, 2));
        [U,S,~] = svd(firingData - mean(firingData, 2), 'econ');
        
        princComp = U*S;
       % princComp = U;
   
        %% Do Linear Discriminant Analysis
        [Sw,Sb] = getLDAscatters(firingData',dirLabels);
    
        %pcaNum = 100; 
      
        % Calculate the Fisher's criterion to use the most Disciriminant
        % Feature method

        fisher=((princComp(:,1:pcaNum)'*Sw*princComp(:,1:pcaNum))^-1 )*(princComp(:,1:pcaNum)'*Sb*princComp(:,1:pcaNum));
        [evectsF, evalsF] = eig(fisher);

        % Obtain the order of eigenvalues from largest to smallest
        [~,sortedInd] = sort(diag(evalsF),'descend');

        % Optimum output
        basis = princComp(:,1:pcaNum)*evectsF(:,sortedInd(1:ldaNum));
        
        w = basis'*(firingData - mean(firingData,2));
%         scatter(w(1,:),w(2,:))

        % Calculating centroids of LDA space for use in nearest centroid
        % classification
        centroids = zeros(ldaNum, 8);
        for ix=1:8
            centroids(:,ix) = mean(w(:,(trainPerc*(ix-1)+1):trainPerc*ix),2);
        end
    
        % store values
        modelParameters.DimReduction(binNum).wLDA = w;
        modelParameters.DimReduction(binNum).PCAdim = pcaNum;
        modelParameters.DimReduction(binNum).LDAdim = ldaNum;
        modelParameters.DimReduction(binNum).LDAout = basis;
        modelParameters.DimReduction(binNum).meanFiring = mean(firingData,2);
        modelParameters.DimReduction(binNum).Centroids = centroids;
    
        binNum = binNum+1;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xPath,yPath] = getAvPath(data,dirNum,trainPerc)

    pathLengths = getPathLengths(data,dirNum,trainPerc);
    maxPath = max(pathLengths);
    % minSize = min(pathLengths);
    
    xPath = zeros(dirNum, maxPath, trainPerc);
    yPath = zeros(dirNum, maxPath, trainPerc);
    pathTraj = zeros(2, maxPath);
    
    % Since trajectories have different length, repeat last element for the
    % shortest trajectories until all trajectories have the same size as
    % maxPath

    for a = 1:dirNum
        for b = 1:trainPerc
            handPos = data((a - 1) * trainPerc + b).handPos(1:2,:);
            handPosEnd = handPos(:, end);
            pathTraj = [handPos repmat(handPosEnd,1, maxPath - pathLengths((a - 1) * trainPerc + b))];
            xPath(a, :, b) = pathTraj(1, :);
            yPath(a, :, b) = pathTraj(2, :);
        end
    end
    xPath=squeeze(mean(xPath,3));
    yPath = squeeze(mean(yPath,3));

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

function [Sw,Sb] = getLDAscatters(X,y)
% X is a data matrix where each row represents an observation and each column represents a feature
% y is a vector of class labels corresponding to each observation in X

    % separate the data into different classes
    classes = unique(y);
    numClasses = length(classes);
%     classData = cell(numClasses, 1);

    % compute the grand mean, within-class and between-class scatter matrix
    GrandMean = mean(X);
    Sw = zeros(size(X, 2));
    Sb = zeros(size(X, 2));
    for a = 1:numClasses
        classData= X(y == classes(a), :);
        classMean = mean(classData);
        classScatter = (classData - classMean)' * (classData - classMean);
        Sw = Sw + classScatter;

        classSize = size(classData, 1);
        betweenScatter = classSize * (classMean - GrandMean)' * (classMean - GrandMean);
        Sb = Sb + betweenScatter;
    end

end


function pathLengths = getPathLengths(data,dirNum,trainPerc)
    pathLengths = zeros(1,trainPerc*dirNum);
    for a=1:dirNum
        for b=1:trainPerc
            pathLengths((a-1)*trainPerc+b) = size(data(b,a).spikes,2);
        end
    end
end


% function [firingData] = flattenData(binnedData)
%     neuralNum = size(binnedData,1);
%     bins = size(binnedData,2);
%     trials = size(binnedData,3);
%     % Flatten data, rows: (neurons x bins) columns: trial
%             for c = 1: trials
%                 firingData(:,c) = reshape(binnedData(:,:,c)',[],1);  
%             end
% end