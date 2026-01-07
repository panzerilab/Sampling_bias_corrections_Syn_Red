function [edges_spikes1, edges_spikes2] = simulate_edges(nRbins, n_neurons, n_stimuli, rateStim, sharedRate_Stim, n_trials, iscorr)

%% Generate the simualted un/correlated data of pairs of neurons
if iscorr
    duration = 1; % Duration of spike trains in seconds
    stdJitter = 5e-3; % Standard deviation of jitter in seconds (5 ms)
    dt = 1e-3; % Time step in seconds (1 ms)
    time = 0:dt:duration-dt; % Time vector
    % Initialize spike trains array
    spike_trains = zeros(n_stimuli, n_neurons, n_trials, length(time));
    shared_spike_trains = zeros(n_stimuli, n_neurons, n_trials, length(time));
    meanRates_trials = zeros(n_neurons, n_stimuli, n_trials); % Preallocate for efficiency
    for trial = 1:n_trials
        for idx_stim = 1:n_stimuli
            for idx_cell = 1:n_neurons
                % Generate spikes for each neuron independently across all stimuli
                spike_trains(idx_stim, idx_cell, trial, :) = poissrnd(rateStim{idx_cell}(idx_stim) * dt, [1, length(time)]);
            end
            % Generating shared spikes
            sharedSpikes_Stim = poissrnd(sharedRate_Stim(idx_stim)*dt, [1, length(time)]);

            % Generating the correlated spike trains across neurons
            [shared_spike_trains(idx_stim, 1, trial, :),shared_spike_trains(idx_stim, 2, trial, :)] ...
                = corr_neuron_pairs(time,dt,stdJitter,sharedSpikes_Stim,squeeze(spike_trains(idx_stim, 1, trial, :))', squeeze(spike_trains(idx_stim, 2, trial, :))');
        end
    end
    for trial = 1:n_trials
        for idx_stim = 1:n_stimuli
            for idx_cell = 1:n_neurons
                % Calculate and display mean firing rates for each stimulus
                meanRates_trials(idx_stim,idx_cell,trial) = sum(squeeze(shared_spike_trains(idx_stim, idx_cell, trial, :))) / duration;
            end
        end
    end
    meanRates_spikes_ALL = cell(1, n_neurons);
    for idx_cell = 1:n_neurons
        firing_rate_across_trials_across_stim = [];
        for idx_stim = 1:n_stimuli
            firing_rate_across_trials = squeeze(meanRates_trials(idx_stim,idx_cell,:));
            firing_rate_across_trials_across_stim = [firing_rate_across_trials_across_stim; firing_rate_across_trials];
        end
        meanRates_spikes_ALL{idx_cell} = firing_rate_across_trials_across_stim;
    end
    meanRates_spikes1 = meanRates_spikes_ALL{1};
    meanRates_spikes2 = meanRates_spikes_ALL{2};

    stimuli = [];
    for idx_stim = 1:n_stimuli
        stimuli = [stimuli, idx_stim*ones(1, n_trials)];
    end
else
    duration = 1; % Duration of spike trains in seconds
    dt = 1e-3; % Time step in seconds (1 ms)
    time = 0:dt:duration-dt; % Time vector
    meanRates_trials = zeros(n_stimuli, n_neurons, n_trials); % Preallocate for efficiency
    % Initialize spike trains array
    spike_trains = zeros(n_stimuli, n_neurons, n_trials, length(time));
    for trial = 1:n_trials
        for idx_stim = 1:n_stimuli
            for idx_cell = 1:n_neurons
                % Generate spikes for each neuron independently for each stimulus
                spike_trains(idx_stim, idx_cell, trial, :) = poissrnd(rateStim{idx_cell}(idx_stim) * dt, [1, length(time)]);

                % Calculate and display mean firing rates for each stimulus
                meanRates_trials(idx_stim,idx_cell,trial) = sum(spike_trains(idx_stim, idx_cell, trial, :)) / duration;
            end
        end
    end

    meanRates_spikes_ALL = cell(1, n_neurons);
    for idx_cell = 1:n_neurons
        firing_rate_across_trials_across_stim = [];
        for idx_stim = 1:n_stimuli
            firing_rate_across_trials = squeeze(meanRates_trials(idx_stim,idx_cell,:));
            firing_rate_across_trials_across_stim = [firing_rate_across_trials_across_stim; firing_rate_across_trials];
        end
        meanRates_spikes_ALL{idx_cell} = firing_rate_across_trials_across_stim;
    end

    meanRates_spikes1 = meanRates_spikes_ALL{1};
    meanRates_spikes2 = meanRates_spikes_ALL{2};

    stimuli = [];
    for idx_stim = 1:n_stimuli
        stimuli = [stimuli, idx_stim*ones(1, n_trials)];
    end
end

%% Binning

[~, edges_spikes1] = binr(meanRates_spikes1', nRbins, 'eqpop');
[~, edges_spikes2] = binr(meanRates_spikes2', nRbins, 'eqpop');
edges_spikes1 = edges_spikes1{1};
edges_spikes2 = edges_spikes2{1};
edges_spikes1(1) = 0;
edges_spikes2(1) = 0;
x = nRbins+1;
edges_spikes1(x) = Inf;
edges_spikes2(x) = Inf;

end