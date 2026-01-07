function [info_results_pairwise_PEAK, bin_size, selected_trials] = pairwise_neuron_PEAK_information_correct_vs_incorrect(spike_train_alltrials, deconv_type, type_data, type_trials)
    %% This function works only after performing the analysis to extract the neurons_class_GC
    load([pwd,'/DATA/metadata.mat']);
    load([pwd,'/DATA/information_metadata.mat']);
    bin_size                                                             = information_metadata(1, type_data, deconv_type).bin_size;
    type_sig_cells                                                    = 'Sig_cells';    
    load(['DATA/Neurons_class_GC_sliding_window_',char(type_sig_cells),'_',num2str(bin_size),'tp.mat'],'Neurons_class','II_CLASS_PEAK_value_mice'); 
    load([pwd,'/DATA/significant_neurons_sliding_window_',num2str(bin_size),'tp.mat']);
    info_results_pairwise_PEAK                               = struct;
    selected_trials                                                    = cell(1, 34);
    for ind                                                                = 1:34
        if length(significant_neurons{ind})               > 19
            ind
            bin_size                                                     = information_metadata(ind, type_data, deconv_type).bin_size;
            R_bins_pairwise                                         = information_metadata(ind, type_data, deconv_type).R_bins_pairwise;
            t_end                                                         = information_metadata(ind, type_data, deconv_type).t_end;                               
            stimulus_vector                                         = metadata(ind,type_data).stimulus_vector;
            choice_vector                                            = metadata(ind,type_data).choice_vector;
            c_vec                                                         = metadata(ind,type_data).c_vec;
            spike_train                                                 = spike_train_alltrials(ind,type_data,deconv_type).spike_train_alltrials(:,c_vec'~=0,:);
            n_bins_trials                                              = t_end - bin_size;
            S_R_info_windows_PEAK                           = cell(1,n_bins_trials-1);
            C_R_info_windows_PEAK                           = cell(1,n_bins_trials-1);        
            S_R_info_btsp_windows_PEAK                  = cell(1,n_bins_trials-1);
            C_R_info_btsp_windows_PEAK                  = cell(1,n_bins_trials-1);
            % Added this three lines of code to consider only the incorrect trials
            switch type_trials
                case 'Incorrect'
                    incorrect1                                                     = stimulus_vector == 1 & choice_vector == 2;
                    incorrect2                                                     = stimulus_vector == 2 & choice_vector == 1;
                    selected_trials{ind}                                              =  incorrect1 | incorrect2;
                case 'Correct'
                    incorrect1                                                     = stimulus_vector == 1 & choice_vector == 2;
                    incorrect2                                                     = stimulus_vector == 2 & choice_vector == 1;
                    n_sample_trials                                            = sum(incorrect1 | incorrect2);      
                    target1                                                        = stimulus_vector == 1 & choice_vector == 1;
                    target2                                                            = stimulus_vector == 2 & choice_vector == 2;
                    correct_trials                                                   = target1 | target2;
                    index_correct_trials                                         =  find(correct_trials == 1);
                    index_randsample                                       = randsample(length(index_correct_trials), n_sample_trials);
                    correct_trials_sampled                                =  zeros(length(stimulus_vector),1);
                    correct_trials_sampled(index_correct_trials(index_randsample)) = 1;    
                    selected_trials{ind}                                             = logical(correct_trials_sampled);
            end
            N_trials                                                     = size(spike_train,2);
            N_neurons                                                = length(significant_neurons{ind}); %% Significant neurons
            S                                                               = stimulus_vector';
            S                                                               = S(selected_trials{ind});
            C                                                               = choice_vector';
            C                                                               = C(selected_trials{ind});
            index_pairwise                                          = nchoosek(1:N_neurons,2);
            N                                                               = length(index_pairwise);
            R_pairwise                                                 = zeros(length(index_pairwise), N_trials);
            for i                                                           = 1:length(index_pairwise)
                N1                                                         = significant_neurons{ind}(index_pairwise(i,1));
                N2                                                         = significant_neurons{ind}(index_pairwise(i,2));
                R1                                                         = mean(spike_train(II_CLASS_PEAK_value_mice{ind}(N1): II_CLASS_PEAK_value_mice{ind}(N1)+bin_size, :, N1),1);
                R2                                                          = mean(spike_train(II_CLASS_PEAK_value_mice{ind}(N2): II_CLASS_PEAK_value_mice{ind}(N2)+bin_size, :, N2),1);
                R1(R1>0)                                               = 1;
                R2(R2>0)                                               = 1;
                R1                                                          = string(R1(:)');
                R2                                                          = string(R2(:)');
                R1(strcmp(R1(1,:),'1'))                           = 'a';
                R2(strcmp(R2(1,:),'1'))                           = 'b';               
                R_pairwise(i, strcmp(R1,'0') & strcmp(R2,'0'))                  = 0;
                R_pairwise(i, strcmp(R1,'a') & strcmp(R2,'0'))                  = 1;
                R_pairwise(i, strcmp(R1,'0') & strcmp(R2,'b'))                  = 2;
                R_pairwise(i, strcmp(R1,'a') & strcmp(R2,'b'))                  = 3;
            end
            R_pairwise                                                  = R_pairwise(:, selected_trials{ind});
            R1                                                              = R1(:, selected_trials{ind});
            R2                                                              = R2(:, selected_trials{ind});
            [S_R_info_windows_PEAK, C_R_info_windows_PEAK, S_R_info_btsp_windows_PEAK, C_R_info_btsp_windows_PEAK]       = SI_CI_code(S,R_pairwise,C,N);
            info_results_pairwise_PEAK(ind, type_data, deconv_type).S_R_info_windows                            = S_R_info_windows_PEAK;
            info_results_pairwise_PEAK(ind, type_data, deconv_type).C_R_info_windows                            = C_R_info_windows_PEAK;
            info_results_pairwise_PEAK(ind, type_data, deconv_type).S_R_info_btsp_windows                   = S_R_info_btsp_windows_PEAK;
            info_results_pairwise_PEAK(ind, type_data, deconv_type).C_R_info_btsp_windows                   = C_R_info_btsp_windows_PEAK;
        end
    end
end