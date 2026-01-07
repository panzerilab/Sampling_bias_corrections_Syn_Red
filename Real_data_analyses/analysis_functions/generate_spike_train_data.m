function generate_spike_train_data
    %% DFF deconvolution (only brightcells; remove NaN trials)
    processed_original_data                     = dir(['Data/data_for_information_analysis/*.mat']);  %%% Consider all experiments and mice
    type_data_name                              = {'corrected','smoothed','FIR','sgolay'};
    DFF_alltrials                               = struct;
    spike_train_alltrials                       = struct;
    deconv_type_name                            = {'foopsi_ar1','foopsi_ar2'};
    for deconv_type                             = 1:length(deconv_type_name)
        deconv_type
        for ind                                 = 1:length(processed_original_data)  
            ind
            for type_data                       = 1:length(type_data_name)
                load(['Data/data_for_information_analysis/', char(processed_original_data(ind).name)]);       %%% Load file with Fluorescence responses, stimuli and choices
                DFF_data_Beh                                                = {DFF_corrected_Beh, DFF_smoothed_Beh, DFF_FIR_Beh, DFF_sgolay_Beh};
                nan_trials                                                                 = sum(sum(isnan(DFF_data_Beh{type_data}(:,:,brightcells == 1))),3);
                if sum(sum(sum(isnan(DFF_data_Beh{type_data}(:, nan_trials == 0, brightcells == 1))))) == 0  
                    
                    DFF_alltrials(ind,type_data).DFF_alltrials                             = DFF_data_Beh{type_data}(:, nan_trials == 0, brightcells == 1);
                    T                                                                      = size(DFF_alltrials(ind,type_data).DFF_alltrials,1);               %%% Timepoints
                    N                                                                      = size(DFF_alltrials(ind,type_data).DFF_alltrials,3);               %%% Neurons
                    DFF_alltrials(ind,type_data).name_exp_mice                             = name_exp_mice;
                    DFF_alltrials(ind,type_data).date_exp_mice                             = date_exp_mice;
                    DFF_alltrials(ind,type_data).N                                         = N;
                    DFF_alltrials(ind,type_data).type_data_name                            = type_data_name(type_data);

                    spike_train_alltrials(ind,type_data,deconv_type).spike_train_alltrials = deconv_DFF(DFF_alltrials(ind,type_data).DFF_alltrials, deconv_type_name{deconv_type});
                    spike_train_alltrials(ind,type_data,deconv_type).name_exp_mice         = name_exp_mice;
                    spike_train_alltrials(ind,type_data,deconv_type).date_exp_mice         = date_exp_mice;
                    spike_train_alltrials(ind,type_data,deconv_type).N                     = N;
                    spike_train_alltrials(ind,type_data,deconv_type).deconv_type           = deconv_type_name{deconv_type};
                end
                
            end
        end
    end
    mkdir(['Results/']);
    save(['Results/current_analysis/DFF_spktrain_alltrials.mat'],'DFF_alltrials','spike_train_alltrials','-v7.3'); 
end