% 
% 
% ks=[5:5:50];
% seeds = [2013, 2023, 50, 6074, 1998];
% times = ones(10, 5, 5); %(k values, seed values, trials)
% rmses = ones(10, 5, 5);
% 
% for ki=1:10
%     for trial = [1:5]
%         for si = 1:5
% 
%             k = ks(ki);
%             seed = seeds(si);
%             disp("k value: " + num2str(k) + " trial: " + num2str(trial) + " seed: " + num2str(si));
%             close all
%             tic 
%             RMSE = boobs(k,seed);
%             time_elapsed = toc;
%             times(ki, si, trial) = time_elapsed;
%             rmses(ki, si, trial) = RMSE;
% 
%         end
%     end
% end


figure
k = [5:5:50]
for ix = [1:5]
subplot(3,2,ix)
scatter(k, rmses(:,ix,1))
end
