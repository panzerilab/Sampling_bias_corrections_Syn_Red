clc, clear all
data_folder = fullfile('Data');
results_folder = fullfile('../Results/Real_data_analysis');
if ~isfolder(results_folder)
    mkdir(results_folder)
end
figures_folder = fullfile('../Figures/Figure_4');
if ~isfolder(figures_folder)
    mkdir(figures_folder)
end

load(fullfile(data_folder,'CA1_data.mat'));

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

widxow_frame_number = 5;

Results_plugin  = struct();
Results_shuffSub   = struct();
Results_qe     = struct();
Results_qeshuffSub = struct();

subs = fieldnames(neuro_struct);
numel(subs)
for idx = 1:numel(subs)
    idx
    data    = neuro_struct.(subs{idx});
    S       = data.Positions;
    opts_bin.bin_method = {'eqpop'};
    opts_bin.n_bins = {12};
    S = cell2mat(binning({S}, opts_bin));

    R_all   = data.Binary;
    [nNeurons,nTrials] = size(R_all);
    % -------------------------------------------------------------------------
    % Step 2: Compute PID
    % -------------------------------------------------------------------------
    pairlist = nchoosek(1:nNeurons,2);
    numPairs = length(pairlist);

    % Full Data
    disp('Full Data')
    PID_v_plugin_all     = zeros(numPairs, 4);
    PID_v_qe_all         = zeros(numPairs, 4);
    PID_v_qeShuffSub_all = zeros(numPairs, 4);
    PID_v_shuffSub_all   = zeros(numPairs, 4);
    parfor pairi = 1:length(pairlist)
        cell1 = pairlist(pairi,1);
        cell2 = pairlist(pairi,2);

        resp1 = R_all(cell1,:);
        resp2 = R_all(cell2,:);
      
        PID_v_plugin_tmp   = PID({resp1,resp2,S}, {'PID_atoms'}, opts);
        PID_v_qe_tmp       = PID({resp1,resp2,S}, {'PID_atoms'}, opts_qe);
        PID_v_shuffSub_tmp = PID({resp1,resp2,S}, {'PID_atoms'}, opts_shuffSub);

        plugin     = [PID_v_plugin_tmp{1}, PID_v_plugin_tmp{2}, PID_v_plugin_tmp{3}, PID_v_plugin_tmp{4}];
        qe         = [PID_v_qe_tmp{1}, PID_v_qe_tmp{2}, PID_v_qe_tmp{3}, PID_v_qe_tmp{4}];
        shuffSub   = [PID_v_shuffSub_tmp{1}, PID_v_shuffSub_tmp{2}, PID_v_shuffSub_tmp{3}, PID_v_shuffSub_tmp{4}];
        qeShuffSub = mean([PID_v_qe_tmp{1}, PID_v_qe_tmp{2}, PID_v_qe_tmp{3}, PID_v_qe_tmp{4}; ...
            PID_v_shuffSub_tmp{1}, PID_v_shuffSub_tmp{2}, PID_v_shuffSub_tmp{3}, PID_v_shuffSub_tmp{4}],1);

        PID_v_plugin_all(pairi, :)      = [sum(plugin), plugin(1), plugin(2), sum(plugin(3:4))];
        PID_v_qe_all(pairi,:)           = [sum(qe), qe(1), qe(2), sum(qe(3:4))];
        PID_v_shuffSub_all(pairi, :)    = [sum(shuffSub), shuffSub(1), shuffSub(2), sum(shuffSub(3:4))];
        PID_v_qeShuffSub_all(pairi, :)  = [sum(qeShuffSub), qeShuffSub(1), qeShuffSub(2), sum(qeShuffSub(3:4))];
    end
    Results_plugin(idx).FullData.Result = PID_v_plugin_all;
    Results_qe(idx).FullData.Result = PID_v_qe_all;
    Results_shuffSub(idx).FullData.Result = PID_v_shuffSub_all;
    Results_qeshuffSub(idx).FullData.Result = PID_v_qeShuffSub_all;

    % Half Data
    disp('Half Data')
    [S_part, R_part] = split_trials(S, R_all);
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
                PID_v_shuffSub_tmp{1}, PID_v_shuffSub_tmp{2}, PID_v_shuffSub_tmp{3}, PID_v_shuffSub_tmp{4}],1);

            PID_v_plugin_all(pairi, iter, :)     = [sum(plugin), plugin(1), plugin(2), sum(plugin(3:4))];
            PID_v_qe_all(pairi, iter, :)         = [sum(qe), qe(1), qe(2), sum(qe(3:4))];
            PID_v_shuffSub_all(pairi, iter, :)   = [sum(shuffSub), shuffSub(1), shuffSub(2), sum(shuffSub(3:4))];
            PID_v_qeShuffSub_all(pairi, iter, :) = [sum(qeShuffSub), qeShuffSub(1), qeShuffSub(2), sum(qeShuffSub(3:4))];
        end
    end
    Results_plugin(idx).HalfData.Result = PID_v_plugin_all;
    Results_qe(idx).HalfData.Result = PID_v_qe_all;
    Results_shuffSub(idx).HalfData.Result = PID_v_shuffSub_all;
    Results_qeshuffSub(idx).HalfData.Result = PID_v_qeShuffSub_all;
end
save(fullfile(results_folder,'CA1_PID.mat'),  "Results_plugin", "Results_qe", "Results_shuffSub", "Results_qeshuffSub");

%% Helper

function [S_part, R_part] = split_trials(S, R_all)
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

S_part{1} = S(idx1);
S_part{2}= S(idx2);
R_part{1} = R_all(:, idx1);
R_part{2} = R_all(:, idx2);
end


