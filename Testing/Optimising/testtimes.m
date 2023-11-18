%% Optimising parameters
% Can change bin size, pca_dim, lda_dim, or k if using kNN.
% Make sure right classifier is being used in 'PositionEstimator - i.e., 
% kNN or centroids.
% just run section you need


%% Testing firing threshold
seeds = [2013, 2023, 50, 6074, 1998, 2, 12343, 8, 1000000, 1]; %select any random seed list. length of this list = number of trials (splits).
firing_thresholds = [0:0.05:0.5]; %adjust to whatever range you want
k = 35;
pca_dim = 100;
lda_dim = 4;
bin_size = 80;
rmses = ones(size(firing_thresholds,2), size(seeds,2)); %(bins, trials)
times = ones(size(firing_thresholds,2), size(seeds,2)); %(bins, trials)

  for fi = 1:size(firing_thresholds,2)
        for si = 1:size(seeds, 2)
            firing_threshold = firing_thresholds(fi);
            seed = seeds(si);
            disp("firing threshold: " + num2str(firing_threshold) + " trial: " + num2str(si));
            close all
            tic 
            RMSE = boobs(k,seed, pca_dim, lda_dim, bin_size, firing_threshold);
            time_elapsed = toc;
            times(fi,si) = time_elapsed;
            rmses(fi, si) = RMSE;

        end
    end


%% Looking at bin sizes (skip this section if not using)
seeds = [2013, 2023, 50, 6074, 1998, 2, 12343, 8, 1000000, 1]; %select any random seed list. length of this list = number of trials (splits).
bin_sizes = [20, 40, 80]; 
k = 35; %this doesn't do anything if using nearest centroid
pca_dim = 100;
lda_dim = 4;
firing_threshold = 0.1;
rmses = ones(size(bin_sizes,2), size(seeds,2)); %(bins, trials)
times = ones(size(bin_sizes,2), size(seeds,2)); %(bins, trials)

    for bi = 1:size(bin_sizes,2)
        for si = 1:size(seeds, 2)
            bin_size = bin_sizes(bi);
           
            seed = seeds(si);
            disp("bin size: " + num2str(bin_size) + " trial: " + num2str(si));
            close all
            tic 
            RMSE = boobs(k, seed, pca_dim, lda_dim, bin_size, firing_threshold);
            time_elapsed = toc;
            times(bi,si) = time_elapsed;
            rmses(bi, si) = RMSE;

        end
    end

%% Looking at PCA by itself (skip if not using)
seeds = [2013, 2023, 50, 6074, 1998, 2, 12343, 8, 1000000, 1]; %select any random seed list. length of this list = number of trials (splits).
PCA_dims = [20:20:200];
bin_size = 80;
lda_dim = 4;
firing_threshold = 0.1;
k = 35; %doesn't do anything if using centroids
rmses = ones(size(PCA_dims,2), size(seeds,2)); %(pca_dim, trials)
times = ones(size(PCA_dims,2), size(seeds,2)); %(pca_dim, trials)

    for pi = 1:size(PCA_dims,2)
        for si = 1:size(seeds, 2) 
            pca_dim = PCA_dims(pi);
            seed = seeds(si);
            disp("PCA dimension: " + num2str(pca_dim) + " trial: " + num2str(si));
            close all
            tic 
            RMSE = boobs(k, seed, pca_dim, lda_dim, bin_size, firing_threshold);
            time_elapsed = toc;
            times(pi,si) = time_elapsed;
            rmses(pi, si) = RMSE;

        end
    end

%% Looking at LDA by itself (skip if not using)
seeds = [2013, 2023, 50, 6074, 1998, 2, 12343, 8, 1000000, 1]; %select any random seed list. length of this list = number of trials (splits).
LDA_dims = [1:7];
bin_size = 80;
pca_dim = 100;
firing_threshold = 0.1;
k = 35; %doesn't do anything if using centroids
rmses = ones(size(LDA_dims,2), size(seeds,2)); %(pca_dim, trials)
times = ones(size(LDA_dims,2), size(seeds,2)); %(pca_dim, trials)

    for li = 1:size(LDA_dims,2)
        for si = 1:size(seeds, 2) 
            lda_dim = LDA_dims(li);
            seed = seeds(si);
            disp("LDA dimension: " + num2str(lda_dim) + " trial: " + num2str(si));
            close all
            tic 
            RMSE = boobs(k, seed, pca_dim, lda_dim, bin_size, firing_threshold);
            time_elapsed = toc;
            times(li,si) = time_elapsed;
            rmses(li, si) = RMSE;

        end
    end


%% Looking at LDA and PCA together(skip if not using)
seeds = [2013, 2023, 50, 6074, 1998, 2, 12343, 8, 1000000, 1]; %select any random seed list. length of this list = number of trials (splits).
LDA_dims = [1:7];
bin_size = 80;
PCA_dims = [20:20:200];
firing_threshold = 0.1;
k = 35; %doesn't do anything if using centroids
rmses = ones(size(PCA_dims,2), size(LDA_dims,2), size(seeds,2)); %(pca_dim, lda_dim, trials)
times = ones(size(PCA_dims,2), size(LDA_dims,2), size(seeds,2)); %(pca_dim, lda_dim, trials)

    for pi = 1:size(PCA_dims,2)
        for li = 1:size(LDA_dims,2)
            for si = 1:size(seeds, 2) 
                pda_dim = PCA_dims(pi);
                lda_dim = LDA_dims(li);
                seed = seeds(si);
                disp("PCA dimension: " + num2str(pda_dim) + ", LDA dimension: " + num2str(lda_dim) + ", trial: " + num2str(si));
                close all
                tic 
                RMSE = boobs(k, seed, pca_dim, lda_dim, bin_size, firing_threshold);
                time_elapsed = toc;
                times(pi, li,si) = time_elapsed;
                rmses(pi, li, si) = RMSE;
    
            end
        end
    end

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55

%% For plotting the PCA/LDA joint one you'll probably need to do heat map or some 3d thing 

%% Plot errors

figure
bins = [20,40,80];
for ix = [1:10]
subplot(5,2,ix)
%scatter(bins, rmses(:,ix))
plot(bins,rmses(:,ix),'-x')
xlabel('bin size')
ylabel('error')

end
sgtitle('error vs bin size for kNN')

%% Plot times
figure
bins = bin_sizes;
for ix = [1:10]
subplot(5,2,ix)
%scatter(bins, rmses(:,ix))
plot(bins,times(:,ix),'-x')
xlabel('bin size')
ylabel('time')

end
sgtitle('time vs bin size for centroid')

%% Plot averages
figure
bins = [20 40 80];
subplot(1,2,1)
plot(bins, mean(rmses,2));
xlabel('bin size')
ylabel('rmse')
title('average error vs bin size for centroid')
subplot(1,2,2)
plot(bins, mean(times,2))
xlabel('bin size')
ylabel('time')
title('average time vs bin size for centroid')
