function trial_selection_dataset(data_all_electrodes_all_sessions_all_noiselevels, delta, n_datasets)
%% Full data
    % pool the trials across the noise levels 
    data_pooled_trials_allsessions = cell(1, n_datasets);
    for idx_data = 1:n_datasets
        data_pooled_trials = [];
        for idx_noise_level = 1:4
            current_data = data_all_electrodes_all_sessions_all_noiselevels{idx_noise_level}{idx_data};
            data_pooled_trials = cat(2, data_pooled_trials, current_data);  % If numeric, this concatenates row-wise
        end
        data_pooled_trials_allsessions{idx_data} = data_pooled_trials;
    end
    
    %% Half data
    % pool the trials across the noise levels 
    data_pooled_trials_allsessions_half1 = cell(1, n_datasets);
    data_pooled_trials_allsessions_half2 = cell(1, n_datasets);
    for idx_data = 1:n_datasets
        data_pooled_trials_half1 = [];
        data_pooled_trials_half2 = [];
        for idx_noise_level = 1:4
            current_data = data_all_electrodes_all_sessions_all_noiselevels{idx_noise_level}{idx_data};
            N = size(current_data,2);
            allidx = 1:N;
            halfCount = ceil(N/2);
            idx_half1 = sort(allidx(randperm(N, halfCount)));
            idx_half2 = sort(setdiff(allidx, idx_half1));            
            data_pooled_trials_half1 = cat(2, data_pooled_trials_half1, current_data(:,idx_half1,:));  % If numeric, this concatenates row-wise
            data_pooled_trials_half2 = cat(2, data_pooled_trials_half2, current_data(:,idx_half2,:));  % If numeric, this concatenates row-wise
        end
        data_pooled_trials_allsessions_half1{idx_data} = data_pooled_trials_half1;
        data_pooled_trials_allsessions_half2{idx_data} = data_pooled_trials_half2;
    end  

    save(['selected_trials_windowSize_',num2str(delta),'.mat'],'data_pooled_trials_allsessions','data_pooled_trials_allsessions_half1','data_pooled_trials_allsessions_half2');
end