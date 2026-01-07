function bias_correction_on_data_hipp_et_al_2022(neuro_struct, bias_correction_ALL, redundancy_measure)
%% Bias correction 
    pid_opts.bin_methodY = 'eqpop';
    pid_opts.n_binsY = 12;
    pid_opts.bin_methodX1 = 'none';
    pid_opts.bin_methodX2 = 'none';
    pid_opts.old_output = 1;
    if strcmp(redundancy_measure,'I_min') == 1
        pid_opts.function = @pidimin;
    end    
    if strcmp(redundancy_measure,'I_MMI') == 1
        pid_opts.function = @pidimmi;
    end        
    pid_opts.redundancy_measure = redundancy_measure;    
    
%     pid_opts.btsp = 1;
%     pid_opts.btsp_variables = {'Y'};
%     pid_opts.btsp_type = {'all'};    
    subject_names = fieldnames(neuro_struct);
    info_results_pairwise_PEAK_PID  = struct;
    for ind = 1:length(subject_names)
        ind
        data = neuro_struct.(subject_names{ind});
        n_neurons = length(data.ROI);
        n_trials_sessions = length(data.Positions);
        [allpairs] = nchoosek(1:n_neurons,2);
        PID_shared = cell(1,length(allpairs));
        PID_uniqueX1 = cell(1,length(allpairs));
        PID_uniqueX2 = cell(1,length(allpairs));
        PID_complementary = cell(1,length(allpairs));   
        switch bias_correction_ALL
            case 'naive'
                pid_opts.bias  = 'naive'; 
                tic
                parfor pair = 1 : length(allpairs)
                    N1 = allpairs(pair,1);
                    N2 = allpairs(pair,2);
                    PID_allpairs_S_R1R2 = PID(data.Positions,data.Binary(N1,:),data.Binary(N2,:),pid_opts);
                    PID_shared{pair} = PID_allpairs_S_R1R2.biased.shared.value;
                    PID_uniqueX1{pair} = PID_allpairs_S_R1R2.biased.uniqueX1.value;
                    PID_uniqueX2{pair} = PID_allpairs_S_R1R2.biased.uniqueX2.value;
                    PID_complementary{pair} = PID_allpairs_S_R1R2.biased.complementary.value;              
                end                
                toc
            case 'qe'
                pid_opts.bias  = 'qe'; 
                pid_opts.xtrp = 1;
%                 pid_opts.max_draws_per_split_number = 5;
                tic
                parfor pair = 1 : length(allpairs)
                    N1 = allpairs(pair,1);
                    N2 = allpairs(pair,2);
                    PID_allpairs_S_R1R2 = PID(data.Positions,data.Binary(N1,:),data.Binary(N2,:),pid_opts);
                    PID_shared{pair} = PID_allpairs_S_R1R2.unbiased.shared.value;
                    PID_uniqueX1{pair} = PID_allpairs_S_R1R2.unbiased.uniqueX1.value;
                    PID_uniqueX2{pair} = PID_allpairs_S_R1R2.unbiased.uniqueX2.value;
                    PID_complementary{pair} = PID_allpairs_S_R1R2.unbiased.complementary.value;              
                end                    
                toc
            case 'btsp_correction'
                pid_opts.bias  = 'naive'; 
                input_info_type = data.Positions;
                tic
                parfor pair = 1:length(allpairs)
                    N1 = allpairs(pair,1);
                    N2 = allpairs(pair,2);
                    PID_allpairs_S_R1R2 = PID(input_info_type,data.Binary(N1,:),data.Binary(N2,:),pid_opts);
                    PID_shared_naive = PID_allpairs_S_R1R2.biased.shared.value;
                    PID_uniqueX1_naive = PID_allpairs_S_R1R2.biased.uniqueX1.value;
                    PID_uniqueX2_naive = PID_allpairs_S_R1R2.biased.uniqueX2.value;
                    PID_complementary_naive = PID_allpairs_S_R1R2.biased.complementary.value;   

                    input_info_type_rnd = input_info_type(randperm(length(input_info_type)));
                    [PID_allpairs_S_R1R2] = PID(input_info_type_rnd,data.Binary(N1,:),data.Binary(N2,:),pid_opts);
                    PID_shared_btsp = PID_allpairs_S_R1R2.biased.shared.value;
                    PID_uniqueX1_btsp = PID_allpairs_S_R1R2.biased.uniqueX1.value;
                    PID_uniqueX2_btsp = PID_allpairs_S_R1R2.biased.uniqueX2.value;
                    PID_complementary_btsp = PID_allpairs_S_R1R2.biased.complementary.value;     

                    PID_shared{pair} = PID_shared_naive - PID_shared_btsp;
                    PID_uniqueX1{pair} = PID_uniqueX1_naive - PID_uniqueX1_btsp;
                    PID_uniqueX2{pair} = PID_uniqueX2_naive - PID_uniqueX2_btsp;
                    PID_complementary{pair} = PID_complementary_naive - PID_complementary_btsp;                       
                end         
                toc
            case 'bootstrapped_qe'
                pid_opts.bias  = 'qe'; 
                input_info_type = data.Positions;
                pid_opts.xtrp = 1;
                tic
                parfor pair = 1:length(allpairs)
                    N1 = allpairs(pair,1);
                    N2 = allpairs(pair,2);                
   
                    PID_allpairs_S_R1R2 = PID(input_info_type,data.Binary(N1,:),data.Binary(N2,:),pid_opts);
                    PID_shared_qe = PID_allpairs_S_R1R2.unbiased.shared.value;
                    PID_uniqueX1_qe = PID_allpairs_S_R1R2.unbiased.uniqueX1.value;
                    PID_uniqueX2_qe = PID_allpairs_S_R1R2.unbiased.uniqueX2.value;
                    PID_complementary_qe = PID_allpairs_S_R1R2.unbiased.complementary.value;  

                    input_info_type_rnd = input_info_type(randperm(length(input_info_type)));
                    PID_allpairs_S_R1R2 = PID(input_info_type_rnd,data.Binary(N1,:),data.Binary(N2,:),pid_opts);
                    PID_shared_btsp = PID_allpairs_S_R1R2.unbiased.shared.value;
                    PID_uniqueX1_btsp = PID_allpairs_S_R1R2.unbiased.uniqueX1.value;
                    PID_uniqueX2_btsp = PID_allpairs_S_R1R2.unbiased.uniqueX2.value;
                    PID_complementary_btsp = PID_allpairs_S_R1R2.unbiased.complementary.value;  

                    PID_shared{pair} = PID_shared_qe - PID_shared_btsp;
                    PID_uniqueX1{pair} = PID_uniqueX1_qe - PID_uniqueX1_btsp;
                    PID_uniqueX2{pair} = PID_uniqueX2_qe - PID_uniqueX2_btsp;
                    PID_complementary{pair} = PID_complementary_qe - PID_complementary_btsp;                
                end
                toc 
                
                tic
                case 'venkatesh'
                    pid_opts.bias  = 'naive'; 
                    input_info_type = data.Positions;
                    parfor pair = 1 : length(allpairs)
                        N1 = allpairs(pair,1);
                        N2 = allpairs(pair,2);      

                        opts = []; opts.bias = 'naive';opts.method = 'dr'; opts.bin_methodX = 'none'; opts.bin_methodY = 'eqpop'; opts.verbose = false ; opts.n_binsX = 2; opts.n_binsY = 12;
                        outputs = information([data.Binary(N1,:); data.Binary(N2,:)], input_info_type, opts, {'I'});
                        Ijoint = outputs{1}; 

                        opts = []; opts.bias = 'qe';opts.method = 'dr'; opts.bin_methodX = 'none'; opts.bin_methodY = 'eqpop'; opts.verbose = false ; opts.n_binsX = 2; opts.n_binsY = 12;
                        outputs = information([data.Binary(N1,:); data.Binary(N2,:)], input_info_type, opts, {'I'});
                        Ijoint_qe = outputs{1};        
                        Ijoint_qe = max([0, Ijoint_qe]);            

                        debias_factor = Ijoint_qe/Ijoint;

                        opts = []; opts.bias = 'qe';opts.method = 'dr'; opts.bin_methodX = 'none'; opts.bin_methodY = 'eqpop'; opts.verbose = false ; opts.n_binsX = 2; opts.n_binsY = 12;
                        outputs = information([data.Binary(N1,:)], input_info_type, opts, {'I'});
                        imx_qe = outputs{1};  
                        imx_qe = max([0, imx_qe]);            

                        opts = []; opts.bias = 'qe';opts.method = 'dr'; opts.bin_methodX = 'none'; opts.bin_methodY = 'eqpop'; opts.verbose = false ; opts.n_binsX = 2; opts.n_binsY = 12;
                        outputs = information([data.Binary(N2,:)], input_info_type, opts, {'I'});
                        imy_qe = outputs{1};     
                        imy_qe = max(0, imy_qe);                        

                        Ijoint_qe = max([Ijoint_qe, imx_qe, imy_qe]);

                        PID_allpairs_S_R1R2 = PID(input_info_type, data.Binary(N1,:), data.Binary(N2,:), pid_opts);
                        PID_shared_raw = PID_allpairs_S_R1R2.biased.shared.value;    
                        PID_uniqueX1_raw = PID_allpairs_S_R1R2.biased.uniqueX1.value;
                        PID_uniqueX2_raw = PID_allpairs_S_R1R2.biased.uniqueX2.value;

                        union_info = ( PID_shared_raw + PID_uniqueX1_raw + PID_uniqueX2_raw ) * debias_factor;
                        union_info = max([union_info, imx_qe, imy_qe]);
                        union_info = min([union_info, imx_qe + imy_qe, Ijoint_qe]);

                        PID_uniqueX1{pair} = union_info - imy_qe;
                        PID_uniqueX2{pair} = union_info - imx_qe;
                        PID_shared{pair} = imx_qe + imy_qe - union_info;
                        PID_complementary{pair} = Ijoint_qe - union_info;
                    end
                    toc                              
        end  
        info_results_pairwise_PEAK_PID(ind).n_trials_sessions = n_trials_sessions;
        info_results_pairwise_PEAK_PID(ind).PID_shared = PID_shared;
        info_results_pairwise_PEAK_PID(ind).PID_uniqueX1 = PID_uniqueX1;
        info_results_pairwise_PEAK_PID(ind).PID_uniqueX2 = PID_uniqueX2;
        info_results_pairwise_PEAK_PID(ind).PID_complementary = PID_complementary;        
    end
    save(['Results/info_results_pairwise_PEAK_PID_hipp',char(bias_correction_ALL),'.mat'],'info_results_pairwise_PEAK_PID');
end