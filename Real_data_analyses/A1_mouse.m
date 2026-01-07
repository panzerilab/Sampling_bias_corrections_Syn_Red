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
        R_all = spike_train_sig_cells_peak_info_ind{ind};
        R_all(nan_trials_in_spike_trains{ind},:) = [];
        R = ((R_all > 0) + (R_all > 1))';

        R_prestim = spike_train_sig_cells_peak_info_ind{ind};
        R_prestim(nan_trials_in_spike_trains{ind},:) = [];
        Rp = ((R_prestim > 0) + (R_prestim > 1))';
        
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
        parfor pairi = 1:length(pairlist)
            cell1 = pairlist(pairi,1);
            cell2 = pairlist(pairi,2);
            resp1 = R(cell1,:);
            resp2 = R(cell2,:);
            PID_v_plugin_tmp   = PID({resp1,resp2,S}, {'PID_atoms'}, opts);
            PID_v_qe_tmp       = PID({resp1,resp2,S}, {'PID_atoms'}, opts_qe);
            PID_v_shuffSub_tmp = PID({resp1,resp2,S}, {'PID_atoms'}, opts_shuffSub);
            
            plugin     = [PID_v_plugin_tmp{1}, PID_v_plugin_tmp{2}, PID_v_plugin_tmp{3}, PID_v_plugin_tmp{4}];
            qe         = [PID_v_qe_tmp{1}, PID_v_qe_tmp{2}, PID_v_qe_tmp{3}, PID_v_qe_tmp{4}];
            shuffSub   = [PID_v_shuffSub_tmp{1}, PID_v_shuffSub_tmp{2}, PID_v_shuffSub_tmp{3}, PID_v_shuffSub_tmp{4}];
            qeShuffSub = mean([PID_v_qe_tmp{1}, PID_v_qe_tmp{2}, PID_v_qe_tmp{3}, PID_v_qe_tmp{4}; ...
                 PID_v_shuffSub_tmp{1}, PID_v_shuffSub_tmp{2}, PID_v_shuffSub_tmp{3}, PID_v_shuffSub_tmp{4}],2);
            PID_v_plugin_all(pairi, :)      = [sum(plugin), plugin(1), plugin(2), sum(plugin(3:4))];
            PID_v_qe_all(pairi,:)           = [sum(qe), qe(1), qe(2), sum(qe(3:4))];
            PID_v_shuffSub_all(pairi, :)    = [sum(shuffSub), shuffSub(1), shuffSub(2), sum(shuffSub(3:4))];
            PID_v_qeShuffSub_all(pairi, :)  = [sum(qeShuffSub), qeShuffSub(1), qeShuffSub(2), sum(qeShuffSub(3:4))];
        end
        Results_plugin(ind).FullData.Result = PID_v_plugin_all;
        Results_qe(ind).FullData.Result = PID_v_qe_all;
        Results_shuffSub(ind).FullData.Result = PID_v_shuffSub_all;
        Results_qeshuffSub(ind).FullData.Result = PID_v_qeShuffSub_all;
    
        % Half Data
        [S_part, R_part] = split_trials(S, R);
        numIterations = length(S_part);
        PID_v_plugin_all     = zeros(numPairs, numIterations, 4);
        PID_v_qe_all         = zeros(numPairs, numIterations, 4);
        PID_v_qeShuffSub_all = zeros(numPairs, numIterations, 4);
        PID_v_shuffSub_all   = zeros(numPairs, numIterations, 4);

        for iter = 1:length(S_part)
            S_iter = S_part{iter};
            R_iter = R_part{iter};
            parfor pairi = 1:length(pairlist)
                cell1 = pairlist(pairi,1);
                cell2 = pairlist(pairi,2);

                resp1 = R_iter(cell1,:);
                resp2 = R_iter(cell2,:);
                PID_v_plugin_tmp   = PID({resp1,resp2,S_iter}, {'PID_atoms'}, opts);
                PID_v_qe_tmp       = PID({resp1,resp2,S_iter}, {'PID_atoms'}, opts_qe);
                PID_v_shuffSub_tmp = PID({resp1,resp2,S_iter}, {'PID_atoms'}, opts_shuffSub);

                plugin     = [PID_v_plugin_tmp{1}, PID_v_plugin_tmp{2}, PID_v_plugin_tmp{3}, PID_v_plugin_tmp{4}];
                qe         = [PID_v_qe_tmp{1}, PID_v_qe_tmp{2}, PID_v_qe_tmp{3}, PID_v_qe_tmp{4}];
                shuffSub   = [PID_v_shuffSub_tmp{1}, PID_v_shuffSub_tmp{2}, PID_v_shuffSub_tmp{3}, PID_v_shuffSub_tmp{4}];
                qeShuffSub = mean([PID_v_qe_tmp{1}, PID_v_qe_tmp{2}, PID_v_qe_tmp{3}, PID_v_qe_tmp{4}; ...
                    PID_v_shuffSub_tmp{1}, PID_v_shuffSub_tmp{2}, PID_v_shuffSub_tmp{3}, PID_v_shuffSub_tmp{4}],2);

                PID_v_plugin_all(pairi, iter, :)     = [sum(plugin), plugin(1), plugin(2), sum(plugin(3:4))];
                PID_v_qe_all(pairi, iter, :)         = [sum(qe), qe(1), qe(2), sum(qe(3:4))];
                PID_v_shuffSub_all(pairi, iter, :)   = [sum(shuffSub), shuffSub(1), shuffSub(2), sum(shuffSub(3:4))];
                PID_v_qeShuffSub_all(pairi, iter, :) = [sum(qeShuffSub), qeShuffSub(1), qeShuffSub(2), sum(qeShuffSub(3:4))];
            end
        end
        Results_plugin(ind).HalfData.Result = PID_v_plugin_all;
        Results_qe(ind).HalfData.Result = PID_v_qe_all;      
        Results_shuffSub(ind).HalfData.Result = PID_v_shuffSub_all;      
        Results_qeshuffSub(ind).HalfData.Result = PID_v_qeShuffSub_all;
    end
end
save(fullfile(results_folder,'A1_PID.mat'), "Results_plugin", "Results_qe", "Results_shuffSub", "Results_qeshuffSub");

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

