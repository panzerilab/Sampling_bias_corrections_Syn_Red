clc, clear all
% Load the stimuli and choices across trials
load('Data/Data_A1/metadata.mat');
deconv_type = 1; type_data = 2; bin_size = 10;
load(['Data/Data_A1/significant_neurons_sliding_window_',num2str(bin_size),'tp.mat'],'significant_neurons');
load(['Data/Data_A1/DFF_spktrain_alltrials.mat'],'spike_train_alltrials');
load(['Data/Data_A1/preprocessed_spk_input_data.mat'],'spike_train_sig_cells_peak_info_ind');
load(['Data/Data_A1/preprocessed_spk_input_data.mat'],'spike_train_sig_cells_prestim_ind');
results_folder = fullfile('../Results/Real_data_analysis');
pop_size_list = 1:19;
max_combinations = 100;
nSessions = 34;
nan_trials_in_spike_trains = cell(1,nSessions);
for ind = 1:nSessions
    if  length(significant_neurons{ind}) > 20
        c_vec            = metadata(ind,type_data).c_vec;
        spike_train = spike_train_alltrials(ind,type_data,deconv_type).spike_train_alltrials(:,c_vec' ~= 0,:);
        nan_trials_in_spike_trains{ind} = find((sum(sum(isnan(spike_train),1),3)));
    end
end


% -------------------------------------------------------------------------
% Define optional PID
% ------------------------------------------------------------------------
% information opts
opts.bias = 'plugin';
opts.bin_method = {'none'};
opts.suppressWarnings = true;
opts.computeNulldist = false;

opts_qe.bias = 'qe';
opts_qe.xtrp = 20;
opts_qe.bin_method = {'none'};
opts_qe.suppressWarnings = true;
opts_qe.computeNulldist = false;

opts_shuffSub.bias = 'shuffSub';
opts_shuffSub.shuff = 20;
opts_shuffSub.bin_method = {'none'};
opts_shuffSub.suppressWarnings = true;
opts_shuffSub.computeNulldist = false;

Results_plugin  = struct();
Results_shuffSub   = struct();
Results_qe     = struct();
Results_qeshuffSub = struct();

field_names = {};
for ind = 1:nSessions
    if  length(significant_neurons{ind}) > 20
        field_name = sprintf('session_%d', ind);
        field_names = [field_names, {field_name}];
        disp(field_name);
        S = metadata(ind,type_data).stimulus_vector;
        S(nan_trials_in_spike_trains{ind}) = [];
        R_all = spike_train_sig_cells_prestim_ind{ind};
        R_all(nan_trials_in_spike_trains{ind},:) = [];
        R = ((R_all > 0) + (R_all > 1))';

        [nNeurons, nTrials] = size(R);
        nPairs = nNeurons*(nNeurons-1)/2;
        % -------------------------------------------------------------------------
        % Step 2: Compute PID 
        % -------------------------------------------------------------------------
        pairlist = nchoosek(1:nNeurons,2);
        numPairs = length(pairlist);

        % Full Data
        PID_v_plugin_all     = zeros(numPairs, 4);
        PID_v_qe_all         = zeros(numPairs, 4);
        PID_v_qeShuffSub_all = zeros(numPairs, 4);
        PID_v_shuffSub_all   = zeros(numPairs, 4);
        
        MI_v_plugin_all     = zeros(numPairs, 4);
        MI_v_qe_all         = zeros(numPairs, 4);
        MI_v_shuffSub_all   = zeros(numPairs, 4);
        MI_v_qeShuffSub_all = zeros(numPairs, 4);
        
        for pairi = 1:length(pairlist)
            cell1 = pairlist(pairi,1);
            cell2 = pairlist(pairi,2);
            resp1 = R(cell1,:);
            resp2 = R(cell2,:);
            % PID_v_plugin_tmp   = PID({resp1,resp2,S}, {'PID_atoms'}, opts);
            % PID_v_qe_tmp       = PID({resp1,resp2,S}, {'PID_atoms'}, opts_qe);
            % PID_v_shuffSub_tmp = PID({resp1,resp2,S}, {'PID_atoms'}, opts_shuffSub);

            % MI_v_plugin_tmp   = MI({[resp1;resp2],S}, {'I(A;B)x'}, opts);
            MI_v_qe_tmp       = MI({[resp1;resp2],S}, {'I(A;B)'}, opts_qe);
            MI_v_shuffSub_tmp = MI({[resp1;resp2],S}, {'I(A;B)'}, opts_shuffSub);
            
            plugin     = [  PID_v_plugin_tmp{1}, PID_v_plugin_tmp{2}, PID_v_plugin_tmp{3}, PID_v_plugin_tmp{4}];
            qe         = [      PID_v_qe_tmp{1}, PID_v_qe_tmp{2}, PID_v_qe_tmp{3}, PID_v_qe_tmp{4}];
            shuffSub   = [PID_v_shuffSub_tmp{1}, PID_v_shuffSub_tmp{2}, PID_v_shuffSub_tmp{3}, PID_v_shuffSub_tmp{4}];
            qeShuffSub = mean([ PID_v_qe_tmp{1}, PID_v_qe_tmp{2}, PID_v_qe_tmp{3}, PID_v_qe_tmp{4}; ...
                 PID_v_shuffSub_tmp{1}, PID_v_shuffSub_tmp{2}, PID_v_shuffSub_tmp{3}, PID_v_shuffSub_tmp{4}],1);
            PID_v_plugin_all(pairi, :)      = [                     MI_v_plugin_tmp{1},   plugin(1),   plugin(2),   sum(plugin(3:4))];
            PID_v_qe_all(pairi,:)           = [                         MI_v_qe_tmp{1},       qe(1),       qe(2),       sum(qe(3:4))];
            PID_v_shuffSub_all(pairi, :)    = [                   MI_v_shuffSub_tmp{1}, shuffSub(1), shuffSub(2), sum(shuffSub(3:4))];
            PID_v_qeShuffSub_all(pairi, :)  = [(MI_v_qe_tmp{1}+MI_v_shuffSub_tmp{1})/2, qeShuffSub(1), qeShuffSub(2), sum(qeShuffSub(3:4))];
        end
        Results_plugin(ind).FullData.Result = PID_v_plugin_all;
        Results_qe(ind).FullData.Result = PID_v_qe_all;
        Results_shuffSub(ind).FullData.Result = PID_v_shuffSub_all;
        Results_qeshuffSub(ind).FullData.Result = PID_v_qeShuffSub_all;
    end
end
save(fullfile(results_folder,'A1_PID_prestim.mat'), "Results_plugin", "Results_qe", "Results_shuffSub", "Results_qeshuffSub");



%% Helper

function [S, R_all] = split_trials(S, R_all)
nTrials = length(S);
nHalf1 = floor(nTrials/2);
nHalf2 = nTrials - nHalf1;

classes = unique(S);
idx1 = [];
idx2 = [];

for i = 1:length(classes)
    current_class = classes(i);
    class_idx = find(S == current_class);
    class_idx = class_idx(randperm(length(class_idx)));
    nClass = length(class_idx);
    nClass1 = round(nClass/2);
    idx1 = [idx1, class_idx(1:nClass1)];
    idx2 = [idx2, class_idx(nClass1+1:end)];
end
idx1 = idx1(randperm(length(idx1)));
idx2 = idx2(randperm(length(idx2)));

S{1} = S(idx1);
S{2}= S(idx2);
R_all{1} = R_all(:, idx1);
R_all{2} = R_all(:, idx2);
end

function plotJointInfo(meanJoint, semJoint, subplotIdx)
    barLabels = {'Plugin', 'qe', 'shuffSub', 'qe-shuffSub'};
    colorIJoint = [0.9290, 0.6940, 0.1250];
    nexttile(subplotIdx);
    b = bar(meanJoint, 'FaceColor', colorIJoint, 'EdgeColor', 'k');
    hold on;
    errorbar(1:length(meanJoint), meanJoint, semJoint, 'k', 'LineStyle', 'none', 'LineWidth', 1.5);
    xticks(1:length(barLabels));
    xticklabels(barLabels);
    ylabel('Information [bits]');
    grid on;
    box on;
    hold off;
end

function plotPIDterms(meanJoint, semJoint, subplotIdx, plotTitle, minV,maxLim, isPercent)
    barLabels = {'Red', 'Syn'};
    colors = [0, 0.4470, 0.7410;  
              0.4660, 0.6740, 0.1880]; 
    nexttile(subplotIdx);
    b = bar(meanJoint, 'FaceColor', 'flat', 'EdgeColor', 'k');
    for i = 1:length(meanJoint)
        b.CData(i, :) = colors(i, :);
    end
    hold on;
    errorbar(1:length(meanJoint), meanJoint, semJoint, 'k', 'LineStyle', 'none', 'LineWidth', 1.5);
    xticks(1:length(barLabels));
    xticklabels(barLabels);
    if isPercent
        ylabel('% of Joint Info');
    else
        ylabel('Information [bits]');
    end
    title(plotTitle, 'Interpreter', 'none');
    ylim([minV, maxLim])
    grid on;
    box on;
    hold off;
end


function add_pvalue(p_values, x1, x2, text_heights, line_heights)
for i = 1:numel(p_values)
    stars = get_stars(p_values(i));
    line([x1(i), x2(i)], [line_heights(i),line_heights(i)], 'Color', 'black', 'LineWidth', 1);
    text((x1(i) + x2(i)) / 2, line_heights(i) + text_heights, stars, 'HorizontalAlignment', 'center', 'Color', 'black', 'FontSize', 12);
end
end

function stars = get_stars(p_value)
if p_value < 0.001
    stars = '***';
elseif p_value < 0.01
    stars = '**';
elseif p_value < 0.05
    stars = '*';
else
    stars = 'n.s';
end
end