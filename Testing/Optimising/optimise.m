%% Looking at LDA and PCA together(skip if not using)
seeds = [2013, 2023, 50, 6074, 1998, 2, 12343, 8, 1000000, 1]; %select any random seed list. length of this list = number of trials (splits).
LDA_dims = [1:7];
bin_size = 80;
PCA_dims = [10:10:300];
k = 35; %doesn't do anything if using centroids
matrix_dimension = length(LDA_dims)*length(PCA_dims);
table = zeros(matrix_dimension,4);
rmses = ones(size(PCA_dims,2), size(LDA_dims,2), size(seeds,2)); %(pca_dim, lda_dim, trials)
times = ones(size(PCA_dims,2), size(LDA_dims,2), size(seeds,2)); %(pca_dim, lda_dim, trials)
row_counter = 1;

    for pix = 1:size(PCA_dims,2)
        for li = 1:size(LDA_dims,2)
            for si = 1:size(seeds, 2) 
                pca_dim = PCA_dims(pix);
                lda_dim = LDA_dims(li);
                seed = seeds(si);
                disp("PCA dimension: " + num2str(pca_dim) + ", LDA dimension: " + num2str(lda_dim) + ", trial: " + num2str(si));
                close all
                tic 
                RMSE = boobs(k, seed, pca_dim, lda_dim, bin_size);
                time_elapsed = toc;
                times(pix, li,si) = time_elapsed;
                rmses(pix, li, si) = RMSE;
                
                table(row_counter, 1) = pca_dim;
                table(row_counter, 2) = lda_dim;
                table(row_counter, 3) = RMSE;
                table(row_counter, 4) = seed;
                row_counter = row_counter + 1;
            end
        end
    end