function plot_bias_correction_on_data_hipp_et_al_2022(neuro_struct, bias_correction_ALL)
%%
    subject_names = fieldnames(neuro_struct);
    load(['Results/info_results_pairwise_PEAK_PID_hipp',char(bias_correction_ALL),'.mat'],'info_results_pairwise_PEAK_PID');
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
    for ind = 1:length(subject_names)
        
%         pooled_n_trials_sessions{ind} = info_results_pairwise_PEAK_PID(ind).n_trials_sessions;
        
        PID_shared_pairs  = info_results_pairwise_PEAK_PID(ind).PID_shared;  
        PID_shared_pairs = vertcat(PID_shared_pairs{:});
        pooled_PID_shared_pairs{ind} = PID_shared_pairs;
        mean_PID_shared_pairs(ind) = mean(PID_shared_pairs);

        PID_uniqueX1_pairs  = info_results_pairwise_PEAK_PID(ind).PID_uniqueX1;
        PID_uniqueX1_pairs = vertcat(PID_uniqueX1_pairs{:});
        pooled_PID_uniqueX1_pairs{ind} = PID_uniqueX1_pairs;
        mean_PID_uniqueX1_pairs(ind) = mean(PID_uniqueX1_pairs);

        PID_uniqueX2_pairs  = info_results_pairwise_PEAK_PID(ind).PID_uniqueX2;
        PID_uniqueX2_pairs = vertcat(PID_uniqueX2_pairs{:});
        pooled_PID_uniqueX2_pairs{ind} = PID_uniqueX2_pairs;
        mean_PID_uniqueX2_pairs(ind) = mean(PID_uniqueX2_pairs);

        pooled_PID_uniqueX1_and_uniqueX2_pairs{ind} = PID_uniqueX1_pairs + PID_uniqueX2_pairs;
        
        PID_complementary_pairs  = info_results_pairwise_PEAK_PID(ind).PID_complementary;
        PID_complementary_pairs = vertcat(PID_complementary_pairs{:});
        pooled_PID_complementary_pairs{ind} = PID_complementary_pairs;
        mean_PID_complementary_pairs(ind) = mean(PID_complementary_pairs);

        joint_info_pairs = PID_shared_pairs+PID_uniqueX1_pairs+PID_uniqueX2_pairs+PID_complementary_pairs;
        mean_joint_info_pairs(ind) = mean(joint_info_pairs);
    end

%     ALL_pooled_PID_shared_pairs = mean_PID_shared_pairs;
%     ALL_pooled_PID_uniqueX1_pairs = mean_PID_uniqueX1_pairs;
%     ALL_pooled_PID_uniqueX2_pairs = mean_PID_uniqueX2_pairs;
%     ALL_pooled_PID_complementary_pairs = mean_PID_complementary_pairs;
%     
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
    
    ALL_DATA = {ALL_pooled_joint_info_pairs, ALL_pooled_PID_shared_pairs, ALL_pooled_PID_complementary_pairs, ALL_pooled_PID_uniqueX1_and_uniqueX2_pairs}; % Example ALL_DATA cell array
    mean_values = [mean(ALL_pooled_joint_info_pairs) mean(ALL_pooled_PID_shared_pairs) mean(ALL_pooled_PID_complementary_pairs) mean(ALL_pooled_PID_uniqueX1_and_uniqueX2_pairs)]
    std_error_values = [2*std(ALL_pooled_joint_info_pairs,'omitnan')/sqrt(length(ALL_pooled_joint_info_pairs)) ...
        2*std(ALL_pooled_PID_shared_pairs,'omitnan')/sqrt(length(ALL_pooled_PID_shared_pairs)) ...
        2*std(ALL_pooled_PID_complementary_pairs,'omitnan')/sqrt(length(ALL_pooled_PID_complementary_pairs)) ...
        2*std(ALL_pooled_PID_uniqueX1_and_uniqueX2_pairs,'omitnan')/sqrt(length(ALL_pooled_PID_uniqueX1_and_uniqueX2_pairs))]
    
%     std_error_values = [std(ALL_pooled_joint_info_pairs,'omitnan') ...
%         std(ALL_pooled_PID_shared_pairs,'omitnan')...
%         std(ALL_pooled_PID_complementary_pairs,'omitnan')...
%         std(ALL_pooled_PID_uniqueX1_pairs,'omitnan')...
%         std(ALL_pooled_PID_uniqueX2_pairs,'omitnan')];    

    bar_handles = bar(mean_values, 'FaceColor', 'flat'); % Create a bar object to hold the handles
%     errorbar(1:5, mean_values, std_error_values);

    for i = 1:length(mean_values)
        bar_handles.CData(i,:) = customColors_info_PID(i,:); % Apply colors
        % Add error bars - one for each bar, with matching color
        eb = errorbar(i, mean_values(i), std_error_values(i), std_error_values(i));
        eb.Color = 'k'; % Set the error bar color
        eb.LineStyle = 'none'; % Remove the line connecting error bars 
    end
    categories = categorical({'Joint', 'RED', 'SYN', 'Unq'});
    set(gca, 'xtick', 1:length(categories), 'xticklabel', categories); % Set the categorical names manually
    
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

    % Number of tests
    m = length(p_values);
    
    % Sort p-values and keep track of the original indices
    [sorted_p_values, sort_idx] = sort(p_values);
    
    % Initialize an array to store the adjusted p-values
    adjusted_p_values = zeros(1, m);
    
    % Apply the Holm-Bonferroni adjustment
    for i = 1:m
        adjusted_p_values(i) = sorted_p_values(i) * (m - i + 1);
    end
    
    % Ensure that the adjusted p-values are not greater than 1
    adjusted_p_values = min(adjusted_p_values, 1);
    
    % Reverse the sorting to get the adjusted p-values in the original order
    final_adjusted_p_values = zeros(1, m);
    final_adjusted_p_values(sort_idx) = adjusted_p_values;

    for pair = 1:length(allpairs)
        if final_adjusted_p_values(pair) <= 1E-3    
            sigstar({[allpairs(pair,1) allpairs(pair,2)]},1E-3);
        elseif final_adjusted_p_values(pair) <= 1E-2
            sigstar({[allpairs(pair,1),allpairs(pair,2)]},1E-2);
        elseif final_adjusted_p_values(pair) <= 0.05
            sigstar({[allpairs(pair,1),allpairs(pair,2)]},0.05);
        end                
    end
    ylim([-0.01, 0.12]);
    
    mkdir(['Figures/Figure_4/']); filename = ['Figures/Figure_4/Curreli_etal_',char(bias_correction_ALL),'.pdf']; exportgraphics(gcf,filename,'ContentType','vector');       
    
%         mean_values = [mean(ALL_pooled_joint_info_pairs) mean(ALL_pooled_PID_shared_pairs) mean(ALL_pooled_PID_complementary_pairs) mean(ALL_pooled_PID_uniqueX1_pairs) mean(ALL_pooled_PID_uniqueX2_pairs)];
%     bar_handles = bar(mean_values, 'FaceColor', 'flat'); % Create a bar object to hold the handles
%     for i = 1:length(mean_values)
%         bar_handles.CData(i,:) = customColors_info_PID(i,:); % Apply colors
%     end
%     categories = categorical({'Joint', 'RED', 'SYN', 'U1','U2'});
%     set(gca, 'xtick', 1:length(categories), 'xticklabel', categories); % Set the categorical names manually
% 
%     ALL_DATA = {ALL_pooled_joint_info_pairs, ALL_pooled_PID_shared_pairs, ALL_pooled_PID_complementary_pairs, ALL_pooled_PID_uniqueX1_pairs, ALL_pooled_PID_uniqueX2_pairs}; % Example ALL_DATA cell array
%     n_PID_terms = length(ALL_DATA);    
%     allpairs = nchoosek(2:n_PID_terms,2);
%     p_allpairs = zeros(size(allpairs, 1), 1);
%     max_y = max(mean_values) + 0.; % Calculate max for plotting significance lines
%     y_offset = 0.05; % Offset for drawing the significance lines
% 
%     for pair_idx = 1:size(allpairs, 1)
%         pair = allpairs(pair_idx, :);
%         [~, p_allpairs(pair_idx), ~, ~] = ttest2(ALL_DATA{pair(1)}, ALL_DATA{pair(2)});
%         x1 = pair(1);
%         x2 = pair(2);
%         y = max_y + (y_offset * pair_idx); % Incremental y level for each test
% 
%         % Annotate significance
%         if p_allpairs(pair_idx) <= 0.05
%             sig = '*';
%             if p_allpairs(pair_idx) <= 0.01
%                 sig = '**';
%                 if p_allpairs(pair_idx) <= 0.001
%                     sig = '***';
%                 end
%             end
%             line([x1, x2], [y, y], 'Color', 'k'); % Draw line
%             text((x1+x2)/2, y, sig, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'Color', 'k');
%         end
%     end
% 
%     ylim([-0.01, max_y + y_offset * size(allpairs, 1)]); % Adjust ylim to accommodate annotations
%     hold off;
%     ylabel('bits');

%     title(['Curreli etal (',char(bias_correction_ALL),')']);

end