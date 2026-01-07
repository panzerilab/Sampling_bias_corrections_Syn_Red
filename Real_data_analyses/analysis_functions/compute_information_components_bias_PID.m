function [I, ILIN, ISS, ICI, ICD, Ish, Joint_SI_PID, PID_shared, PID_complementary, PID_uniqueX1, PID_uniqueX2,Joint_prob_ind, Joint_minQ] ...
    = compute_information_components_bias_PID(info_amount,nRbins, meanRates_spikes1, meanRates_spikes2, stimuli,n_stimuli,bias_correction_ALL,redundancy_measure)
    %% Compute the components of the Information breakdown
    opts = [];
    switch bias_correction_ALL
        case {'naive','venkatesh'}    
            opts.bias = 'naive';
            opts.method = 'dr';  opts.bin_methodX = 'none'; opts.bin_methodY = 'none';  opts.verbose = false;
            opts.n_binsX = nRbins; opts.n_binsY = n_stimuli;
            outputs = information([meanRates_spikes1; meanRates_spikes2], stimuli, opts, {'I', 'ILIN', 'ISS', 'ICI', 'ICD','Ish'});
            I = outputs{1}; ILIN = outputs{2}; ISS = outputs{3}; ICI = outputs{4}; ICD = outputs{5}; Ish = outputs{6};  
        case {'qe'}
            opts.bias = bias_correction_ALL;
            opts.xtrp = 100;
            opts.method = 'dr';  opts.bin_methodX = 'none'; opts.bin_methodY = 'none';  opts.verbose = false;
            opts.n_binsX = nRbins; opts.n_binsY = n_stimuli;
            outputs = information([meanRates_spikes1; meanRates_spikes2], stimuli, opts, {'I', 'ILIN', 'ISS', 'ICI', 'ICD','Ish'});
            I = outputs{1}; ILIN = outputs{2}; ISS = outputs{3}; ICI = outputs{4}; ICD = outputs{5}; Ish = outputs{6};  
        case {'bootstrapped','btsp_correction'}
            opts.bias = 'naive';
            opts.method = 'dr';  opts.bin_methodX = 'none'; opts.bin_methodY = 'none';  opts.verbose = false;
            opts.n_binsX = nRbins; opts.n_binsY = n_stimuli;
            stimuli_rnd = stimuli(randperm(length(stimuli)));
            outputs = information([meanRates_spikes1; meanRates_spikes2], stimuli_rnd, opts, {'I', 'ILIN', 'ISS', 'ICI', 'ICD','Ish'});
            I = outputs{1}; ILIN = outputs{2}; ISS = outputs{3}; ICI = outputs{4}; ICD = outputs{5}; Ish = outputs{6};              
        case {'bootstrapped_qe'}
            opts.bias = 'qe';
            opts.method = 'dr';  opts.bin_methodX = 'none'; opts.bin_methodY = 'none';  opts.verbose = false;
            opts.n_binsX = nRbins; opts.n_binsY = n_stimuli;
            stimuli_rnd = stimuli(randperm(length(stimuli)));
            outputs = information([meanRates_spikes1; meanRates_spikes2], stimuli_rnd, opts, {'I', 'ILIN', 'ISS', 'ICI', 'ICD','Ish'});
            I = outputs{1}; ILIN = outputs{2}; ISS = outputs{3}; ICI = outputs{4}; ICD = outputs{5}; Ish = outputs{6};              
    end


    %% Compute the components of the PID
    opts_PID = [];
    opts_PID.n_binsX1  = nRbins; 
    opts_PID.n_binsX2  = nRbins; 
    opts_PID.n_binsY = n_stimuli;
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
    switch bias_correction_ALL
        case 'naive'
            opts_PID.bias  = bias_correction_ALL; 
            [PID_allpairs_S_R1R2, ~, ~,q_distr, pind_distr,~] =  PID(stimuli, meanRates_spikes1, meanRates_spikes2, opts_PID);

            PID_shared = PID_allpairs_S_R1R2.biased.shared.value;    
            PID_uniqueX1 = PID_allpairs_S_R1R2.biased.uniqueX1.value;
            PID_uniqueX2 = PID_allpairs_S_R1R2.biased.uniqueX2.value;
            PID_complementary = PID_allpairs_S_R1R2.biased.complementary.value; 
            Joint_prob_ind = mutualInformationXYZ(pind_distr);
            Joint_minQ = mutualInformationXYZ(q_distr);
        case 'qe'
            opts_PID.bias  = bias_correction_ALL; 
            opts_PID.xtrp = 1;
            switch info_amount
                case {'low','high'}                
                    n_extrp = 20;
                case 'parameter_space'
                    n_extrp = 20;
                case 'parameter_space_asymptotic'
                    n_extrp = 1;
            end
            PID_shared_extrp = zeros(1,n_extrp);
            PID_uniqueX1_extrp = zeros(1,n_extrp);
            PID_uniqueX2_extrp = zeros(1,n_extrp);
            PID_complementary_extrp = zeros(1,n_extrp);
            for extrp = 1:n_extrp
                PID_allpairs_S_R1R2 = PID(stimuli, meanRates_spikes1, meanRates_spikes2, opts_PID);
                PID_shared_extrp(extrp) = PID_allpairs_S_R1R2.unbiased.shared.value;
                PID_uniqueX1_extrp(extrp) = PID_allpairs_S_R1R2.unbiased.uniqueX1.value;
                PID_uniqueX2_extrp(extrp) = PID_allpairs_S_R1R2.unbiased.uniqueX2.value;
                PID_complementary_extrp(extrp) = PID_allpairs_S_R1R2.unbiased.complementary.value;  
            end
            PID_shared = mean(PID_shared_extrp);
            PID_uniqueX1 = mean(PID_uniqueX1_extrp);
            PID_uniqueX2 = mean(PID_uniqueX2_extrp);
            PID_complementary = mean(PID_complementary_extrp);
            
            Joint_prob_ind = zeros(1,1);
            Joint_minQ = zeros(1,1);
         
        case 'bootstrapped'
            opts_PID.bias  = 'naive'; 
            stimuli_rnd = stimuli(randperm(length(stimuli)));
            [PID_allpairs_S_R1R2, ~, ~,q_distr, pind_distr,~] = PID(stimuli_rnd, meanRates_spikes1, meanRates_spikes2, opts_PID);
            PID_shared = PID_allpairs_S_R1R2.biased.shared.value;    
            PID_uniqueX1 = PID_allpairs_S_R1R2.biased.uniqueX1.value;
            PID_uniqueX2 = PID_allpairs_S_R1R2.biased.uniqueX2.value;
            PID_complementary = PID_allpairs_S_R1R2.biased.complementary.value;       
            Joint_prob_ind = mutualInformationXYZ(pind_distr);
            Joint_minQ = mutualInformationXYZ(q_distr);     
            
        case 'btsp_correction'
            opts_PID.bias  = 'naive'; 
            PID_allpairs_S_R1R2 = PID(stimuli, meanRates_spikes1, meanRates_spikes2, opts_PID);
            PID_shared = PID_allpairs_S_R1R2.biased.shared.value;    
            PID_uniqueX1 = PID_allpairs_S_R1R2.biased.uniqueX1.value;
            PID_uniqueX2 = PID_allpairs_S_R1R2.biased.uniqueX2.value;
            PID_complementary = PID_allpairs_S_R1R2.biased.complementary.value;                         
            switch info_amount
                case {'low','high'}                
                    n_btsp = 20;
                case 'parameter_space'
                    n_btsp = 20;
                case 'parameter_space_asymptotic'
                    n_btsp = 1;
            end
            PID_shared_btsp = zeros(1,n_btsp);
            PID_uniqueX1_btsp = zeros(1,n_btsp);
            PID_uniqueX2_btsp = zeros(1,n_btsp);
            PID_complementary_btsp = zeros(1,n_btsp);
            for btsp = 1:n_btsp
                stimuli_rnd = stimuli(randperm(length(stimuli)));
                PID_allpairs_S_R1R2 = PID(stimuli_rnd, meanRates_spikes1, meanRates_spikes2, opts_PID);
                PID_shared_btsp(btsp) = PID_allpairs_S_R1R2.biased.shared.value;    
                PID_uniqueX1_btsp(btsp) = PID_allpairs_S_R1R2.biased.uniqueX1.value;
                PID_uniqueX2_btsp(btsp) = PID_allpairs_S_R1R2.biased.uniqueX2.value;
                PID_complementary_btsp(btsp) = PID_allpairs_S_R1R2.biased.complementary.value;       
            end
            PID_shared = PID_shared - mean(PID_shared_btsp);
            PID_uniqueX1 = PID_uniqueX1 - mean(PID_uniqueX1_btsp);
            PID_uniqueX2 = PID_uniqueX2 - mean(PID_uniqueX2_btsp);
            PID_complementary = PID_complementary - mean(PID_complementary_btsp);
            
            Joint_prob_ind = zeros(1,1);
            Joint_minQ = zeros(1,1);
            
        case 'bootstrapped_qe'
            opts_PID.bias  = 'qe'; 
            opts_PID.xtrp = 1;
            switch info_amount
                case {'low','high'}                
                    n_extrp = 20;
                case 'parameter_space'
                    n_extrp = 20;
                case 'parameter_space_asymptotic'
                    n_extrp = 1;                    
            end                        
            PID_shared_extrp = zeros(1,n_extrp);
            PID_uniqueX1_extrp = zeros(1,n_extrp);
            PID_uniqueX2_extrp = zeros(1,n_extrp);
            PID_complementary_extrp = zeros(1,n_extrp);
            for extrp = 1:n_extrp
                PID_allpairs_S_R1R2 = PID(stimuli, meanRates_spikes1, meanRates_spikes2, opts_PID);
                PID_shared_extrp(extrp) = PID_allpairs_S_R1R2.unbiased.shared.value;
                PID_uniqueX1_extrp(extrp) = PID_allpairs_S_R1R2.unbiased.uniqueX1.value;
                PID_uniqueX2_extrp(extrp) = PID_allpairs_S_R1R2.unbiased.uniqueX2.value;
                PID_complementary_extrp(extrp) = PID_allpairs_S_R1R2.unbiased.complementary.value;  
            end
            PID_shared = mean(PID_shared_extrp);
            PID_uniqueX1 = mean(PID_uniqueX1_extrp);
            PID_uniqueX2 = mean(PID_uniqueX2_extrp);
            PID_complementary = mean(PID_complementary_extrp);         

            opts_PID.xtrp = 1;            
            PID_shared_extrp = zeros(1,n_extrp);
            PID_uniqueX1_extrp = zeros(1,n_extrp);
            PID_uniqueX2_extrp = zeros(1,n_extrp);
            PID_complementary_extrp = zeros(1,n_extrp);
            for extrp = 1:n_extrp
                stimuli_rnd = stimuli(randperm(length(stimuli)));
                PID_allpairs_S_R1R2 = PID(stimuli_rnd, meanRates_spikes1, meanRates_spikes2, opts_PID);
                PID_shared_extrp(extrp) = PID_allpairs_S_R1R2.unbiased.shared.value;
                PID_uniqueX1_extrp(extrp) = PID_allpairs_S_R1R2.unbiased.uniqueX1.value;
                PID_uniqueX2_extrp(extrp) = PID_allpairs_S_R1R2.unbiased.uniqueX2.value;
                PID_complementary_extrp(extrp) = PID_allpairs_S_R1R2.unbiased.complementary.value;  
            end
            PID_shared_btsp = mean(PID_shared_extrp);
            PID_uniqueX1_btsp = mean(PID_uniqueX1_extrp);
            PID_uniqueX2_btsp = mean(PID_uniqueX2_extrp);
            PID_complementary_btsp = mean(PID_complementary_extrp);
            
            PID_shared = PID_shared - PID_shared_btsp;
            PID_uniqueX1 = PID_uniqueX1 - PID_uniqueX1_btsp;
            PID_uniqueX2 = PID_uniqueX2 - PID_uniqueX2_btsp;
            PID_complementary = PID_complementary - PID_complementary_btsp;
            
            Joint_prob_ind = zeros(1,1);
            Joint_minQ = zeros(1,1);
                        
        case 'venkatesh'
            opts = []; opts.bias = 'naive';opts.method = 'dr'; opts.bin_methodX = 'none'; opts.bin_methodY = 'none'; 
            opts.verbose = false ; opts.n_binsX = nRbins; opts.n_binsY = n_stimuli;
            outputs = information([meanRates_spikes1; meanRates_spikes2], stimuli, opts, {'I'});
            Ijoint = outputs{1}; 
            
            opts = []; opts.bias = 'qe';opts.method = 'dr'; opts.bin_methodX = 'none'; opts.bin_methodY = 'none'; 
            opts.verbose = false ; opts.n_binsX = nRbins; opts.n_binsY = n_stimuli;
            outputs = information([meanRates_spikes1; meanRates_spikes2], stimuli, opts, {'I'});
            Ijoint_qe = outputs{1};        
            Ijoint_qe = max([0, Ijoint_qe]);            
            
            debias_factor = Ijoint_qe/Ijoint;
                        
            opts = []; opts.bias = 'qe';opts.method = 'dr'; opts.bin_methodX = 'none'; opts.bin_methodY = 'none'; 
            opts.verbose = false ; opts.n_binsX = nRbins; opts.n_binsY = n_stimuli;
            outputs = information([meanRates_spikes1], stimuli, opts, {'I'});
            imx_qe = outputs{1};  
            imx_qe = max([0, imx_qe]);            
            
            opts = []; opts.bias = 'qe';opts.method = 'dr'; opts.bin_methodX = 'none'; opts.bin_methodY = 'none'; 
            opts.verbose = false ; opts.n_binsX = nRbins; opts.n_binsY = n_stimuli;
            outputs = information([meanRates_spikes2], stimuli, opts, {'I'});
            imy_qe = outputs{1};     
            imy_qe = max(0, imy_qe);                        
            
            Ijoint_qe = max([Ijoint_qe, imx_qe, imy_qe]);
            
            opts_PID.bias  = 'naive'; 
            PID_allpairs_S_R1R2 = PID(stimuli, meanRates_spikes1, meanRates_spikes2, opts_PID);
            PID_shared = PID_allpairs_S_R1R2.biased.shared.value;    
            PID_uniqueX1 = PID_allpairs_S_R1R2.biased.uniqueX1.value;
            PID_uniqueX2 = PID_allpairs_S_R1R2.biased.uniqueX2.value;
            
            union_info = ( PID_shared + PID_uniqueX1 + PID_uniqueX2 ) * debias_factor;
            union_info = max([union_info, imx_qe, imy_qe]);
            union_info = min([union_info, imx_qe + imy_qe, Ijoint_qe]);
            
            PID_uniqueX1 = union_info - imy_qe;
            PID_uniqueX2 = union_info - imx_qe;
            PID_shared = imx_qe + imy_qe - union_info;
            PID_complementary = Ijoint_qe - union_info;
            
            Joint_prob_ind = zeros(1,1);
            Joint_minQ = zeros(1,1);
    end
    Joint_SI_PID = PID_shared + PID_uniqueX1 + PID_uniqueX2 + PID_complementary;
end