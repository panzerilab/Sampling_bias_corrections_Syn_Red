function [spikes1_with_corrs, spikes2_with_corrs] = corr_neuron_pairs(time,dt,stdJitter,sharedSpikes,spikes1, spikes2)
    %% Jittering shared spikes for the second cell
    jitterTimes = stdJitter*randn(size(sharedSpikes)); % Generate jitter times
    jitteredSpikes2 = zeros(size(spikes2)); % Initialize jittered spikes array for the second cell
    for i = 1:length(sharedSpikes)
        if sharedSpikes(i) > 0
            % Find the appropriate time to add the jittered spike
            jitterTime = jitterTimes(i);
            jitteredIndex = i + round(jitterTime/dt);
            % Ensure the index is within bounds
            if jitteredIndex >= 1 && jitteredIndex <= length(time)
                jitteredSpikes2(jitteredIndex) = jitteredSpikes2(jitteredIndex) + sharedSpikes(i);
            end
        end
    end
    % Adding shared spikes to the first cell directly
    spikes1_with_corrs = spikes1 + sharedSpikes;
    % Adding jittered shared spikes to the second cell
    spikes2_with_corrs = spikes2 + jitteredSpikes2;
    
end