function plot_ALL_correct_incorrect_auditory_cortex_data(info_results_pairwise_PEAK_PID, type_trials, bias_correction_ALL, bin_size, type_data, deconv_type)
%%
    subject_names = 1:34;
    load(['DATA/significant_neurons_sliding_window_',num2str(bin_size),'tp.mat']);
    pooled_n_trials_sessions = cell(1,length(subject_names));
    mean_PID_shared_pairs = zeros(1,length(subject_names));
    mean_PID_uniqueX1_pairs = zeros(1,length(subject_names));
    mean_PID_uniqueX2_pairs = zeros(1,length(subject_names));
    mean_PID_complementary_pairs = zeros(1,length(subject_names));
    mean_joint_info_pairs = zeros(1,length(subject_names));
    pooled_PID_shared_pairs = cell(1,length(subject_names));
    pooled_PID_uniqueX1_pairs = cell(1,length(subject_names));
    pooled_PID_uniqueX2_pairs = cell(1,length(subject_names));
    pooled_PID_uniqueX1_and_uniqueX2_pairs = cell(1,length(subject_names));
    pooled_PID_complementary_pairs = cell(1,length(subject_names));
    considered_sessions = zeros(1,length(subject_names));
    for ind = 1:length(subject_names)
        if length(significant_neurons{ind}) > 19
% % % %             pooled_n_trials_sessions{ind} = info_results_pairwise_PEAK_PID(ind, type_data, deconv_type).n_trials_sessions;
            PID_shared_pairs  = info_results_pairwise_PEAK_PID(ind,type_data, deconv_type).PID_shared;  
            PID_shared_pairs = vertcat(PID_shared_pairs{:});
            pooled_PID_shared_pairs{ind} = PID_shared_pairs;
            mean_PID_shared_pairs(ind) = mean(PID_shared_pairs);

            PID_uniqueX1_pairs  = info_results_pairwise_PEAK_PID(ind,type_data, deconv_type).PID_uniqueX1;
            PID_uniqueX1_pairs = vertcat(PID_uniqueX1_pairs{:});
            pooled_PID_uniqueX1_pairs{ind} = PID_uniqueX1_pairs;
            mean_PID_uniqueX1_pairs(ind) = mean(PID_uniqueX1_pairs);

            PID_uniqueX2_pairs  = info_results_pairwise_PEAK_PID(ind,type_data, deconv_type).PID_uniqueX2;
            PID_uniqueX2_pairs = vertcat(PID_uniqueX2_pairs{:});
            pooled_PID_uniqueX2_pairs{ind} = PID_uniqueX2_pairs;
            mean_PID_uniqueX2_pairs(ind) = mean(PID_uniqueX2_pairs);

            pooled_PID_uniqueX1_and_uniqueX2_pairs{ind} = PID_uniqueX1_pairs + PID_uniqueX2_pairs;
            
            PID_complementary_pairs  = info_results_pairwise_PEAK_PID(ind,type_data, deconv_type).PID_complementary;
            PID_complementary_pairs = vertcat(PID_complementary_pairs{:});
            pooled_PID_complementary_pairs{ind} = PID_complementary_pairs;
            mean_PID_complementary_pairs(ind) = mean(PID_complementary_pairs);

            joint_info_pairs = PID_shared_pairs+PID_uniqueX1_pairs+PID_uniqueX2_pairs+PID_complementary_pairs;
            mean_joint_info_pairs(ind) = mean(joint_info_pairs);
            
            considered_sessions(ind) = 1;
        end
    end

%     ALL_pooled_PID_shared_pairs = mean_PID_shared_pairs(considered_sessions==1);
%     ALL_pooled_PID_uniqueX1_pairs = mean_PID_uniqueX1_pairs(considered_sessions==1);
%     ALL_pooled_PID_uniqueX2_pairs = mean_PID_uniqueX2_pairs(considered_sessions==1);
%     ALL_pooled_PID_complementary_pairs = mean_PID_complementary_pairs(considered_sessions==1);
%     ALL_pooled_joint_info_pairs = ALL_pooled_PID_shared_pairs+ALL_pooled_PID_uniqueX1_pairs+ALL_pooled_PID_uniqueX2_pairs+ALL_pooled_PID_complementary_pairs;
    ALL_pooled_n_trials_sessions = vertcat(pooled_n_trials_sessions{:});
    ALL_pooled_PID_shared_pairs = vertcat(pooled_PID_shared_pairs{:});
    ALL_pooled_PID_uniqueX1_pairs = vertcat(pooled_PID_uniqueX1_pairs{:});
    ALL_pooled_PID_uniqueX2_pairs = vertcat(pooled_PID_uniqueX2_pairs{:});
    ALL_pooled_PID_uniqueX1_and_uniqueX2_pairs = vertcat(pooled_PID_uniqueX1_and_uniqueX2_pairs{:});
    ALL_pooled_PID_complementary_pairs = vertcat(pooled_PID_complementary_pairs{:});
    ALL_pooled_joint_info_pairs = ALL_pooled_PID_shared_pairs+ALL_pooled_PID_uniqueX1_and_uniqueX2_pairs+ALL_pooled_PID_complementary_pairs;
    
    %%    
    figure('Visible','on'); set(gcf,'Position',[00, 00, 300, 200]); hold on;
    customColors_info_PID = ...
        [0.9290, 0.6940, 0.1250;  % Color for 'Joint'
        0, 0.4470, 0.7410;       % Color for 'Red'
        0.4660, 0.6740, 0.1880; % Color for 'Syn'
        0.5, 0.5, 0.5; % Color for 'U1'
        0.4, 0.26, 0.13]; % Color for 'U2'
    
    ALL_DATA = {ALL_pooled_joint_info_pairs, ALL_pooled_PID_shared_pairs, ALL_pooled_PID_complementary_pairs, ALL_pooled_PID_uniqueX1_and_uniqueX2_pairs}; 
    mean_values = [mean(ALL_pooled_joint_info_pairs) mean(ALL_pooled_PID_shared_pairs) mean(ALL_pooled_PID_complementary_pairs) mean(ALL_pooled_PID_uniqueX1_and_uniqueX2_pairs)]
    std_error_values = [2*std(ALL_pooled_joint_info_pairs,'omitnan')/sqrt(length(ALL_pooled_joint_info_pairs)) ...
        2*std(ALL_pooled_PID_shared_pairs,'omitnan')/sqrt(length(ALL_pooled_PID_shared_pairs)) ...
        2*std(ALL_pooled_PID_complementary_pairs,'omitnan')/sqrt(length(ALL_pooled_PID_complementary_pairs)) ...
        2*std(ALL_pooled_PID_uniqueX1_and_uniqueX2_pairs,'omitnan')/sqrt(length(ALL_pooled_PID_uniqueX1_and_uniqueX2_pairs))];
    
    bar_handles = bar(mean_values, 'FaceColor', 'flat'); 
    for i = 1:length(mean_values)
        bar_handles.CData(i,:) = customColors_info_PID(i,:); 
        eb = errorbar(i, mean_values(i), std_error_values(i), std_error_values(i));
        eb.Color = 'k'; 
        eb.LineStyle = 'none'; 
    end
    categories = categorical({'Joint', 'RED', 'SYN', 'Unq'});
    set(gca, 'xtick', 1:length(categories), 'xticklabel', categories); 
    
    p_ttest_onesample = zeros(1, length(mean_values));
    for i = 1:length(mean_values)
        [~, p_ttest_onesample(i)] = ttest(ALL_DATA{i}, 0, 'Tail','right');        
        if p_ttest_onesample(i) <= 0.05
            sigstar({[i, i]},0.05);
        end                     
    end
    
    n_PID_terms = length(ALL_DATA);    
    allpairs = nchoosek(2:n_PID_terms,2);
    p_values = zeros(1, length(allpairs));
    for pair = 1:length(allpairs)
        [~, p_values(pair), ~, ~] = ttest(ALL_DATA{allpairs(pair,1)}, ALL_DATA{allpairs(pair,2)});
    end
    for pair = 1:length(allpairs)
        if p_values(pair) <= 1E-3
            sigstar({[allpairs(pair,1) allpairs(pair,2)]},1E-3);
        elseif p_values(pair) <= 1E-2
            sigstar({[allpairs(pair,1),allpairs(pair,2)]},1E-2);
        elseif p_values(pair) <= 0.05
            sigstar({[allpairs(pair,1),allpairs(pair,2)]},0.05);
        end                
    end
    ylim([-0.01, 0.25]);
    mkdir(['Figures/Figure_4/']); filename = ['Figures/Figure_4/auditory_cortex_data_',char(bias_correction_ALL),'_',char(type_trials),'.pdf']; exportgraphics(gcf,filename,'ContentType','vector');       
end