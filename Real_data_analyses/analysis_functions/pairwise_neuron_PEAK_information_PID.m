function [info_results_pairwise_PEAK_PID, bin_size] = pairwise_neuron_PEAK_information_PID(Rbins, info_type,type_trials,spike_train_alltrials, deconv_type, type_data, bias_correction_ALL, redundancy_measure)
%% This function works only after performing the analysis to extract the neurons_class_GC
    significant_neurons = cell(1,34);
    load([pwd,'/DATA/metadata.mat']);
    load([pwd,'/DATA/information_metadata.mat']);
    bin_size = information_metadata(1, type_data, deconv_type).bin_size;
    type_sig_cells = 'Sig_cells';    
    load([pwd,'/DATA/Neurons_class_GC_sliding_window_',char(type_sig_cells),'_',num2str(bin_size),'tp.mat'],'II_CLASS_PEAK_value_mice'); 
    load([pwd,'/DATA/significant_neurons_sliding_window_',num2str(bin_size),'tp.mat']);
    info_results_pairwise_PEAK_PID  = struct;
    for ind  = 1:34
        if length(significant_neurons{ind}) > 19
            stimulus_vector  = metadata(ind,type_data).stimulus_vector;
            choice_vector  = metadata(ind,type_data).choice_vector;            
            bin_size  = information_metadata(ind, type_data, deconv_type).bin_size;
            switch info_type 
                case 'stimulus_information'
                    input_info_type  = stimulus_vector';
                case 'choice_information'
                    input_info_type  = choice_vector';
                otherwise
                    disp('ERROR');
            end

            % Added this three lines of code to consider only the incorrect trials
            switch type_trials
                case 'ALL'
                    selected_trials{ind} = linspace(1,length(input_info_type),length(input_info_type));
                case 'Incorrect'
                    incorrect1                                                     = stimulus_vector == 1 & choice_vector == 2;
                    incorrect2                                                     = stimulus_vector == 2 & choice_vector == 1;
                    selected_trials{ind}                                              =  incorrect1 | incorrect2;
                    input_info_type = input_info_type(selected_trials{ind});
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
                    input_info_type = input_info_type(selected_trials{ind});
            end
            c_vec            = metadata(ind,type_data).c_vec;
            spike_train      = spike_train_alltrials(ind,type_data,deconv_type).spike_train_alltrials(:,c_vec'~=0,:);
            N_neurons  = length(significant_neurons{ind}); %% Significant neurons
            index_pairwise  = nchoosek(1:N_neurons,2);

            n_trials_sessions = length(selected_trials{ind});
            opts_PID = struct;
            opts_PID.n_binsX1  = Rbins; 
            opts_PID.n_binsX2  = Rbins; 
            opts_PID.n_binsY = 2;
            opts_PID.bin_methodX1 = 'none';
            opts_PID.bin_methodX2 = 'none';
            opts_PID.bin_methodY = 'none';
            opts_PID.old_output = 1;
            if strcmp(redundancy_measure,'I_min') == 1
                opts_PID.function = @pidimin;
            end    
            if strcmp(redundancy_measure,'I_MMI') == 1
                opts_PID.function = @pidimmi;
            end        
            opts_PID.redundancy_measure = redundancy_measure;    
            PID_ALL = cell(1,length(index_pairwise));
            PID_shared = cell(1,length(index_pairwise));
            PID_uniqueX1 = cell(1,length(index_pairwise));
            PID_uniqueX2 = cell(1,length(index_pairwise));
            PID_complementary = cell(1,length(index_pairwise));
            switch bias_correction_ALL
                case 'naive'
                    opts_PID.bias  = 'naive'; 
                    tic
                    switch Rbins
                        case 3
                            parfor i           = 1:length(index_pairwise)
                                N1          = significant_neurons{ind}(index_pairwise(i,1));
                                N2          = significant_neurons{ind}(index_pairwise(i,2));
                                R1          = mean(spike_train(II_CLASS_PEAK_value_mice{ind}(N1): II_CLASS_PEAK_value_mice{ind}(N1)+bin_size, :, N1),1);
                                R2          = mean(spike_train(II_CLASS_PEAK_value_mice{ind}(N2): II_CLASS_PEAK_value_mice{ind}(N2)+bin_size, :, N2),1);
                                R1(R1 > 1) = 2;
                                R1(R1 > 0 & R1 <= 1) = 1;
                                R2(R2 > 1) = 2;
                                R2(R2 > 0 & R2 <= 1) = 1;                        
                                R1          = R1(:, selected_trials{ind});
                                R2          = R2(:, selected_trials{ind});
                                [PID_allpairs_S_R1R2] = PID(input_info_type, R1, R2, opts_PID);
                                PID_shared{i} = PID_allpairs_S_R1R2.biased.shared.value;
                                PID_uniqueX1{i} = PID_allpairs_S_R1R2.biased.uniqueX1.value;
                                PID_uniqueX2{i} = PID_allpairs_S_R1R2.biased.uniqueX2.value;
                                PID_complementary{i} = PID_allpairs_S_R1R2.biased.complementary.value;
                            end
                        case 2
                            parfor i           = 1:length(index_pairwise)
                                N1          = significant_neurons{ind}(index_pairwise(i,1));
                                N2          = significant_neurons{ind}(index_pairwise(i,2));
                                R1          = mean(spike_train(II_CLASS_PEAK_value_mice{ind}(N1): II_CLASS_PEAK_value_mice{ind}(N1)+bin_size, :, N1),1);
                                R2          = mean(spike_train(II_CLASS_PEAK_value_mice{ind}(N2): II_CLASS_PEAK_value_mice{ind}(N2)+bin_size, :, N2),1);
                                R1(R1>0)    = 1;
                                R2(R2>0)    = 1;                       
                                R1          = R1(:, selected_trials{ind});
                                R2          = R2(:, selected_trials{ind});
                                [PID_allpairs_S_R1R2] = PID(input_info_type, R1, R2, opts_PID);
                                PID_shared{i} = PID_allpairs_S_R1R2.biased.shared.value;
                                PID_uniqueX1{i} = PID_allpairs_S_R1R2.biased.uniqueX1.value;
                                PID_uniqueX2{i} = PID_allpairs_S_R1R2.biased.uniqueX2.value;
                                PID_complementary{i} = PID_allpairs_S_R1R2.biased.complementary.value;
                            end
                    end
                    toc
                case 'qe'
                    opts_PID.bias  = 'qe'; 
                    opts_PID.xtrp = 1;
                    tic
                    switch Rbins
                        case 3
                            parfor i           = 1:length(index_pairwise)
                                try
                                    N1          = significant_neurons{ind}(index_pairwise(i,1));
                                    N2          = significant_neurons{ind}(index_pairwise(i,2));
                                    R1          = mean(spike_train(II_CLASS_PEAK_value_mice{ind}(N1): II_CLASS_PEAK_value_mice{ind}(N1)+bin_size, :, N1),1);
                                    R2          = mean(spike_train(II_CLASS_PEAK_value_mice{ind}(N2): II_CLASS_PEAK_value_mice{ind}(N2)+bin_size, :, N2),1);
                                    R1(R1 > 1) = 2;
                                    R1(R1 > 0 & R1 <= 1) = 1;
                                    R2(R2 > 1) = 2;
                                    R2(R2 > 0 & R2 <= 1) = 1;  
                                    R1          = R1(:, selected_trials{ind});
                                    R2          = R2(:, selected_trials{ind});
                                    [PID_allpairs_S_R1R2] = PID(input_info_type, R1, R2, opts_PID);
                                    PID_shared{i} = PID_allpairs_S_R1R2.unbiased.shared.value;
                                    PID_uniqueX1{i} = PID_allpairs_S_R1R2.unbiased.uniqueX1.value;
                                    PID_uniqueX2{i} = PID_allpairs_S_R1R2.unbiased.uniqueX2.value;
                                    PID_complementary{i} = PID_allpairs_S_R1R2.unbiased.complementary.value;       
                                catch
                                    continue;
                                end
                            end                    
                        case 2
                            parfor i           = 1:length(index_pairwise)
                                try
                                    N1          = significant_neurons{ind}(index_pairwise(i,1));
                                    N2          = significant_neurons{ind}(index_pairwise(i,2));
                                    R1          = mean(spike_train(II_CLASS_PEAK_value_mice{ind}(N1): II_CLASS_PEAK_value_mice{ind}(N1)+bin_size, :, N1),1);
                                    R2          = mean(spike_train(II_CLASS_PEAK_value_mice{ind}(N2): II_CLASS_PEAK_value_mice{ind}(N2)+bin_size, :, N2),1);
                                    R1(R1>0)    = 1;
                                    R2(R2>0)    = 1; 
                                    R1          = R1(:, selected_trials{ind});
                                    R2          = R2(:, selected_trials{ind});
                                    [PID_allpairs_S_R1R2] = PID(input_info_type, R1, R2, opts_PID);
                                    PID_shared{i} = PID_allpairs_S_R1R2.unbiased.shared.value;
                                    PID_uniqueX1{i} = PID_allpairs_S_R1R2.unbiased.uniqueX1.value;
                                    PID_uniqueX2{i} = PID_allpairs_S_R1R2.unbiased.uniqueX2.value;
                                    PID_complementary{i} = PID_allpairs_S_R1R2.unbiased.complementary.value;       
                                catch
                                    continue;
                                end
                            end  
                    end
                    toc
                    
                case 'btsp_correction'
                    opts_PID.bias  = 'naive'; 
                    tic
                    parfor i           = 1:length(index_pairwise)
                        N1          = significant_neurons{ind}(index_pairwise(i,1));
                        N2          = significant_neurons{ind}(index_pairwise(i,2));
                        R1          = mean(spike_train(II_CLASS_PEAK_value_mice{ind}(N1): II_CLASS_PEAK_value_mice{ind}(N1)+bin_size, :, N1),1);
                        R2          = mean(spike_train(II_CLASS_PEAK_value_mice{ind}(N2): II_CLASS_PEAK_value_mice{ind}(N2)+bin_size, :, N2),1);

                        R1(R1 > 1) = 2;
                        R1(R1 > 0 & R1 <= 1) = 1;
                        R2(R2 > 1) = 2;
                        R2(R2 > 0 & R2 <= 1) = 1;                         
                        R1          = R1(:, selected_trials{ind});
                        R2          = R2(:, selected_trials{ind});          

                        [PID_allpairs_S_R1R2] = PID(input_info_type, R1, R2, opts_PID);
                        PID_shared_naive = PID_allpairs_S_R1R2.biased.shared.value;
                        PID_uniqueX1_naive = PID_allpairs_S_R1R2.biased.uniqueX1.value;
                        PID_uniqueX2_naive = PID_allpairs_S_R1R2.biased.uniqueX2.value;
                        PID_complementary_naive = PID_allpairs_S_R1R2.biased.complementary.value;
                        
                        input_info_type_rnd = input_info_type(randperm(length(input_info_type)));
                        [PID_allpairs_S_R1R2] = PID(input_info_type_rnd, R1, R2, opts_PID);
                        PID_shared_btsp = PID_allpairs_S_R1R2.biased.shared.value;
                        PID_uniqueX1_btsp = PID_allpairs_S_R1R2.biased.uniqueX1.value;
                        PID_uniqueX2_btsp = PID_allpairs_S_R1R2.biased.uniqueX2.value;
                        PID_complementary_btsp = PID_allpairs_S_R1R2.biased.complementary.value;     

                        PID_shared{i} = PID_shared_naive - PID_shared_btsp;
                        PID_uniqueX1{i} = PID_uniqueX1_naive - PID_uniqueX1_btsp;
                        PID_uniqueX2{i} = PID_uniqueX2_naive - PID_uniqueX2_btsp;
                        PID_complementary{i} = PID_complementary_naive - PID_complementary_btsp;                       
                    end         
                    toc
                    
                case 'bootstrapped_qe'
                    opts_PID.bias  = 'qe'; 
                    opts_PID.xtrp = 1;
                    tic
                    parfor i        = 1:length(index_pairwise)
                        N1          = significant_neurons{ind}(index_pairwise(i,1));
                        N2          = significant_neurons{ind}(index_pairwise(i,2));
                        R1          = mean(spike_train(II_CLASS_PEAK_value_mice{ind}(N1): II_CLASS_PEAK_value_mice{ind}(N1)+bin_size, :, N1),1);
                        R2          = mean(spike_train(II_CLASS_PEAK_value_mice{ind}(N2): II_CLASS_PEAK_value_mice{ind}(N2)+bin_size, :, N2),1);
                        R1(R1 > 1) = 2;
                        R1(R1 > 0 & R1 <= 1) = 1;
                        R2(R2 > 1) = 2;
                        R2(R2 > 0 & R2 <= 1) = 1;                         

                        R1          = R1(:, selected_trials{ind});
                        R2          = R2(:, selected_trials{ind});
                        
                        [PID_allpairs_S_R1R2] = PID(input_info_type, R1, R2, opts_PID);
                        PID_shared_raw = PID_allpairs_S_R1R2.unbiased.shared.value;
                        PID_uniqueX1_raw = PID_allpairs_S_R1R2.unbiased.uniqueX1.value;
                        PID_uniqueX2_raw = PID_allpairs_S_R1R2.unbiased.uniqueX2.value;
                        PID_complementary_raw = PID_allpairs_S_R1R2.unbiased.complementary.value; 
                        
                        input_info_type_rnd = input_info_type(randperm(length(input_info_type)));
                        [PID_allpairs_S_R1R2] = PID(input_info_type_rnd, R1, R2, opts_PID);
                        PID_shared_btsp = PID_allpairs_S_R1R2.unbiased.shared.value;
                        PID_uniqueX1_btsp = PID_allpairs_S_R1R2.unbiased.uniqueX1.value;
                        PID_uniqueX2_btsp = PID_allpairs_S_R1R2.unbiased.uniqueX2.value;
                        PID_complementary_btsp = PID_allpairs_S_R1R2.unbiased.complementary.value; 
                        
                        PID_shared{i} = PID_shared_raw - PID_shared_btsp;
                        PID_uniqueX1{i} = PID_uniqueX1_raw - PID_uniqueX1_btsp;
                        PID_uniqueX2{i} = PID_uniqueX2_raw - PID_uniqueX2_btsp;
                        PID_complementary{i} = PID_complementary_raw - PID_complementary_btsp;  
                                         
                    end                
                    toc
                                  
            end
                    
            info_results_pairwise_PEAK_PID(ind, type_data, deconv_type).n_trials_sessions = n_trials_sessions;
            info_results_pairwise_PEAK_PID(ind, type_data, deconv_type).PID_shared = PID_shared;
            info_results_pairwise_PEAK_PID(ind, type_data, deconv_type).PID_uniqueX1 = PID_uniqueX1;
            info_results_pairwise_PEAK_PID(ind, type_data, deconv_type).PID_uniqueX2 = PID_uniqueX2;
            info_results_pairwise_PEAK_PID(ind, type_data, deconv_type).PID_complementary = PID_complementary;
        end
    end
end