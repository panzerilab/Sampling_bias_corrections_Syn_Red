function [spike_train] = deconv_DFF(DeltaF_trials_stim,deconv_tyoe)
    nan_trials                              = sum(sum(isnan(DeltaF_trials_stim)),3);
    if sum(sum(sum(isnan(DeltaF_trials_stim(:, nan_trials == 0,:))))) == 0    %%% Sometimes there are NaNs in the data.   
        DeltaF_trials_stim                  = DeltaF_trials_stim(:, nan_trials == 0, :);
        T                                   = size(DeltaF_trials_stim,1);               %%% Timepoints
        trials                              = size(DeltaF_trials_stim,2);               %%% Trials
        N                                   = size(DeltaF_trials_stim,3);               %%% Neurons
        spike_train                         = zeros(T, trials, N);                      % deconvolved neural activity
        
        switch deconv_tyoe
            case 'foopsi_ar1'
                disp('FOOPSI AR1');
                for neurons                         = 1:N
                    parfor tr                       = 1:trials                              
                       [~,spike,~]                  = deconvolveCa(DeltaF_trials_stim(:,tr,neurons),'foopsi','ar1');   
                       spike_train(:,tr,neurons)    = spike;
                    end
                end
                
                
            case 'foopsi_ar2'
                disp('FOOPSI AR2');
                for neurons                         = 1:N
                    parfor tr                       = 1:trials                              
                       [~,spike,~]                  = deconvolveCa(DeltaF_trials_stim(:,tr,neurons),'foopsi','ar2');   
                       spike_train(:,tr,neurons)    = spike;
                    end
                end
      
            
        end
        
    else
        delete(['Data/data_for_information_analysis/', char(processed_original_data(ind).name)]);
    end 
                
end
            