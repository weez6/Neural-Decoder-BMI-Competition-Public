function [] = rasterTings(trial, neuronRange)
% Plots raster for given trial, over specified range of the 98 neural units
% e.g. neural units 1 to 10 then neuronRange=[1:10]
% trial is specific trial you want e.g. trial(100)
    
    figure
    hold on
    colors = repmat(['g', 'r', 'b', 'm', 'k'], 1, 20);
    t = [1:size(trial.spikes, 2)];
    for i = neuronRange
        spikes = find(trial.spikes(i,:) ==1);
        for j = 1:length(spikes)
            plot([spikes(j), spikes(j)], [i-0.6, i+0.6], colors(i));
            hold on
        end
    end
    xlabel('Time (ms)')
    ylabel('Neural Unit')
end

