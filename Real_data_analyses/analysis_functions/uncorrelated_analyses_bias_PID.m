function [meanRates_spikes1, meanRates_spikes2, meanRates_trials, stimuli, I, ILIN, ISS, ICI, ICD, Ish,...
    Joint_SI_PID, PID_shared, PID_complementary, PID_uniqueX1, PID_uniqueX2,Joint_prob_ind, Joint_minQ] = uncorrelated_analyses_bias_PID(info_amount,nRbins, bias_correction_ALL, n_neurons, n_stimuli, rateStim, n_trials,redundancy_measure, edges_spikes1, edges_spikes2)

%% Generate the simualted uncorrelated data of pairs of neurons
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

    
    %% Binning

    meanRates_spikes1 = (discretize(meanRates_spikes1, edges_spikes1))';
    meanRates_spikes2 = (discretize(meanRates_spikes2, edges_spikes2))';

        
    %% Compute the components of the Information breakdown    
    [I, ILIN, ISS, ICI, ICD, Ish, Joint_SI_PID, PID_shared, PID_complementary, PID_uniqueX1, PID_uniqueX2,Joint_prob_ind, Joint_minQ] = ...
        compute_information_components_bias_PID(info_amount,nRbins, meanRates_spikes1, meanRates_spikes2, stimuli,n_stimuli,bias_correction_ALL,redundancy_measure);
end