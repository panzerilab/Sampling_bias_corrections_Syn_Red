function generate_dataset_neural_activity_and_stimuli_to_analyze(folder, n_datasets,names_datasets, delta)
    %% Explanations from the readme file
    % % % % data related to noise presentatio (different noise in each trial) 
    % % % % there are different noise levels
    % % % % 4 electrodes per experiment - 6 experiments from 1 monkey
    % % % % Audit.Lfp_data{i}{j} is the LFP of electrode i and stimulus type j (each j is a different noise level. J=1 is no noise, and hihgher j index means higher noise level, see Kayser et al Neuron 2009)
    % % % % Audit.Lfp_data{i}{j}(k,t) is the LFP at trial k and time step t b 1 ms time step
    % % % % 
    % % % % Audit.SPk_data{i,u}{j}{k} contains the list of spike times (in ms) at electrode i, unit u, stimulus type j, trial k
     
    %%% Useful for double-checks
    % for idx_noise_level = [2]
    %     alldata = cell(1,6);
    %     for idx_data = 1:6
    %         load([folder,'/',names_datasets(idx_data).name]);
    %         alldata{idx_data} = Audit.Spk_data;
    %     end
    % end
    
    %% Analyses on the datasets

    % # trials for each noise level (4 levels) and dataset (6 datasets)
    nTrials_sessions = {[30,20,28,24,30,30], [30,20,30,24,30,29], [30,20,30,24,30,30], [30,20,30,24,30,30]};
    
    % The levels of noise added to the auditory stimuli 
    noise_levels = 4;
    data_all_electrodes_all_sessions_all_noiselevels = cell(1,noise_levels);
    for idx_noise_level = 1:noise_levels
        idx_noise_level
        data_all_electrodes_all_sessions = cell(1,n_datasets);
        for idx_data = 1:n_datasets
            load([folder,'/',names_datasets(idx_data).name])
            
            % # of electrodes that were used in each session
            n_electrodes = size(Audit.Spk_data, 1);        
            % # of neurons recorded in each electrode
            n_units = size(Audit.Spk_data, 2);
        
            % # of stimuli
            num_time_windows=size(Audit.Lfp_data{1}{1},2);
            num_stimuli=floor(num_time_windows/delta);
        
            data_all_electrodes = cell(1,n_electrodes);
            for idx_electrodes = 1:n_electrodes
                if ~all(cellfun(@isempty, Audit.Spk_data(idx_electrodes,:))) % chack data availability of neuronal units at each electrode
                    data_all_units = cell(1,n_units);
                    for idx_units = 1:n_units
                        if isempty(Audit.Spk_data{idx_electrodes, idx_units}) == 0 % chack data availability for the specified neural unit within the given electrode
                            if length(Audit.Spk_data{idx_electrodes, idx_units}{idx_noise_level}) == nTrials_sessions{idx_noise_level}(idx_data) % select only data with the selected # of trials
                                % Generate the matrix of spiking activity
                                % across time (binIndices) and trials (n_trials) for each neuron
                                % (idx_units)
                                n_trials = size(Audit.Spk_data{idx_electrodes, idx_units}{idx_noise_level},2);
                                spikeTimes = Audit.Spk_data{idx_electrodes, idx_units}{idx_noise_level};
                                spikeMatrix = zeros(n_trials, num_stimuli);            
                                % Iterate over each trial and bin the spike times
                                for trial = 1:n_trials
                                    % Get spike times for the current trial
                                    spikes = spikeTimes{trial};
                                    % Convert spike times into bin indices
                                    binIndices = ceil(spikes / delta);                
                                    % Make sure the bin indices don't exceed the number of bins
                                    binIndices(binIndices > num_stimuli) = [];               
                                    % Increment the corresponding bin for each spike time
                                    for i = 1:length(binIndices)
                                        spikeMatrix(trial, binIndices(i)) = spikeMatrix(trial, binIndices(i)) + 1;
                                    end
                                end
                                data_all_units{idx_units} = spikeMatrix;
                            else
                                disp(['#trials is:',num2str(length(Audit.Spk_data{idx_electrodes, idx_units}{1}))]);
                            end
                        end
                    end
                    data_all_electrodes{idx_electrodes} = data_all_units;
                    clear data_all_units;
                else
                    disp('do not include those cells')
                end
            end
            % Remove some cells with empty fields
            data_cleaned_no_empty_cells = horzcat(data_all_electrodes{:});
            idx_empty_cells = cellfun(@isempty, data_cleaned_no_empty_cells);
            data_cleaned_no_empty_cells(idx_empty_cells) = [];
    
            data_all_units = data_cleaned_no_empty_cells;
            n_units = length(data_all_units);
            nTrials = size(data_all_units{1},1);
            nWindows = size(data_all_units{1},2);
            % Pool the units in a single dataset to compute pairwise PID -
            % The structure of these dataset are n_units x nTrials x nWindows
            data_all_units_pooled = zeros(n_units, nTrials, nWindows);
            for i = 1:n_units
                data_all_units_pooled(i, :, :) = data_all_units{i};
            end
            data_all_electrodes_all_sessions{idx_data} = data_all_units_pooled;
        end
        data_all_electrodes_all_sessions_all_noiselevels{idx_noise_level} = data_all_electrodes_all_sessions;
    end        
    save('A1_dataset_monkey.mat','data_all_electrodes_all_sessions_all_noiselevels','delta');
end