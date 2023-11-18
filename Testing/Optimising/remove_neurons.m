clear
load 'monkeydata_training.mat'

intermediate_matrix = zeros(98, 1);
neuron_matrix = zeros(size(trial, 1), size(trial, 2));
trial_spikes = trial.spikes;

for row = 1:size(trial_spikes, 1)
    for column = 1:size(trial_spikes, 2)
        spikes = trial_spikes(row, column);
        for each_subrow = 1: size(spikes,1)
            sum_elements = sum(spikes(each_subrow, :));
            intermediate_matrix(1, each_subrow) = sum_elements;
        end
        
        end_position = 0;
        for index = 1:length(intermediate_matrix)
            if intermediate_matrix(index, 1) < 5
                spikes_vector(end_position + 1) = intermediate_matrix(index, 1);
                end_position = end_position + 1;
            end
        end
    end
end