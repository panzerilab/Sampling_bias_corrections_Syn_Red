%% PID lattice + three corrections: demo
% Make sure the following are on your path:
%   pid_lattice.m
%   qele.m     % uses pid.is_gaussian internally
%   shuffsub.m     % uses pid.is_gaussian internally
%   qe.m           % uses pid.is_gaussian internally

rng(0);
clc;
clear;
% addpath(genpath('...')); Add Path to data
load('subjects.mat', 'subj');
if ~iscell(subj)
    subj = cellstr(subj);
end
load('behav_data_idx.mat', 'behav_data_idx');
conditions = {'attn'};
SAVE_DIR = pwd;


for c_idx = 1:length(conditions)
    condition = conditions{c_idx};

    fprintf('Processing condition: %s\n', upper(condition));
    condition_save_dir = fullfile(SAVE_DIR, condition);
    if ~exist(condition_save_dir, 'dir')
        mkdir(condition_save_dir);
    end

    for s_idx = 1:length(subj)
        s = subj{s_idx};
        behav_data = behav_data_idx.(condition).(s);
        epochs = behav_data_idx.(condition).epochs;
        process_subject_matched(s, behav_data, epochs, condition, condition_save_dir);
    end
end
function process_subject_matched(s, behav, epochs, condition, condition_save_dir)
sf = 200;
etimep = -0.5 : 1/sf : 4.5;

epochs_to_process = [1];
d0s = [0.05, 0.20, 0.35, 0.50, 0.65, 0.80, 2.70, 2.85, 3.00, 3.15, 3.30, 3.45];
d0 = arrayfun(@(x) find(abs(etimep - x) == min(abs(etimep - x)), 1), d0s);

% Visual cortex (indices 1–5)
areas1 = {'V1'};
areas2 = {'V2'};
areas3 = {'V3'};


PID_plugin = cell(1,length(epochs_to_process));
PID_res = cell(1,length(epochs_to_process));
PID_shuff = cell(1,length(epochs_to_process));

trials_idx = (behav.upTrials | behav.downTrials); %& behav.genMean4Filt;
did = find(trials_idx);

data_area1 = zeros(length(trials_idx), length(areas1), length(etimep));
data_area2 = zeros(length(trials_idx), length(areas2), length(etimep));
data_area3 = zeros(length(trials_idx), length(areas3), length(etimep));

save_path = fullfile(condition_save_dir, sprintf('PID3Source_%s_%s.mat', s, condition));

for p = 1:length(areas1)
    area_name = ['HCPMMP1_' areas1{p}];
    %pca_filename = fullfile('',s, sprintf('PCA_%s_%s_%s.mat', s, condition, area_name)); ADD PATH TO DATA
    if ~exist(pca_filename, 'file')
        fprintf('Missing file: %s\n', pca_filename);
        continue;
    end
    loaded_data = load(pca_filename, 'score');
    score = loaded_data.score;
    data = cellfun(@(x) permute(x, [3 1 2]), score, 'UniformOutput', false);
    data = cell2mat(data);
    data_tmp = permute(data, [2 3 1]);
    assert(size(data_tmp, 3)==length(etimep))
    data_tmp = zscore(data_tmp, 0, 3); 
    data_area1(:,p,:) =  data_tmp(:,1,:);
end
for p = 1:length(areas2)
    area_name = ['HCPMMP1_' areas2{p}];
    %pca_filename = fullfile('',s, sprintf('PCA_%s_%s_%s.mat', s, condition, area_name)); ADD PATH TO DATA
    if ~exist(pca_filename, 'file')
        fprintf('Missing file: %s\n', pca_filename);
        continue;
    end
    loaded_data = load(pca_filename, 'score');
    score = loaded_data.score;
    data = cellfun(@(x) permute(x, [3 1 2]), score, 'UniformOutput', false);
    data = cell2mat(data);
    data_tmp = permute(data, [2 3 1]);
    assert(size(data_tmp, 3)==length(etimep))
    data_tmp = zscore(data_tmp, 0, 3); 
    data_area2(:,p,:) =  data_tmp(:,1,:);
end

for p = 1:length(areas3)
    area_name = ['HCPMMP1_' areas3{p}];
    %pca_filename = fullfile('',s, sprintf('PCA_%s_%s_%s.mat', s, condition, area_name)); % ADD PATH TO DATA
    if ~exist(pca_filename, 'file')
        fprintf('Missing file: %s\n', pca_filename);
        continue;
    end
    loaded_data = load(pca_filename, 'score');
    score = loaded_data.score;
    data = cellfun(@(x) permute(x, [3 1 2]), score, 'UniformOutput', false);
    data = cell2mat(data);
    data_tmp = permute(data, [2 3 1]);
    assert(size(data_tmp, 3)==length(etimep))
    data_tmp = zscore(data_tmp, 0, 3); 
    data_area3(:,p,:) =  data_tmp(:,1,:);
end

for d_idx = 1:length(epochs_to_process)
    d = epochs_to_process(d_idx);

    R1 = data_area1(did, :, epochs(d, 1):epochs(d, 2));
    R2 = data_area2(did, :, epochs(d, 1):epochs(d, 2));
    R3 = data_area3(did, :, epochs(d, 1):epochs(d, 2));
    S = behav.sample(:, d);S = S(did);

    [~, N_reg1, Nt] = size(R1);
    [~, N_reg2, ~] = size(R2);
    [~, N_reg3, ~] = size(R3);

    n_blocks = floor(Nt / 4);

    PIDttmp_plugin = nan(18,n_blocks);
    PIDttmp_shuff  = nan(18,n_blocks);
    PIDttmp_res = nan(18,n_blocks);

    for b = 1:n_blocks 
        t_start = (b-1) * 4 + 1;
        t_end = t_start + 3;
        R1_avg = mean(R1(:, :, t_start:t_end), 3);
        R2_avg = mean(R2(:, :, t_start:t_end), 3);
        R3_avg = mean(R3(:, :, t_start:t_end), 3);

        opts_bin.bin_method = {'eqpop'};
        opts_bin.n_bins = {3};

        binnedData = binning({R1_avg', R2_avg', R3_avg', S'}, opts_bin);
        data_PID = cat(1, binnedData{1}, binnedData{2}, binnedData{3}, binnedData{4});
        pid_lat = pid_lattice(3, 'IMMI', false);

        idx1 = 1:N_reg1;
        idx2 = N_reg1 + (1:N_reg2);
        idx3 = N_reg1 + N_reg2 + (1:N_reg3);
        idy  = N_reg1 + N_reg2 + N_reg3 + (1:1);
        pid_lat.source_dims = {idx1, idx2, idx3};
        pid_lat.target_dims = idy;

        obsD   = data_PID.';                         % N x V
        nb_vec = max(obsD, [], 1);
        pD     = accumarray(obsD, 1, nb_vec);
        pD     = pD / sum(pD(:));
        lat_emp_discrete = pid_lat.calculate_latvals(pD);

        lat_res_discrete    = qe(data_PID, 10, pid_lat);      % vector
        lat_sh_discrete     = shuffsub(data_PID, 10, pid_lat);      % atoms x nshuff

        PIDttmp_plugin(:,b) = lat_emp_discrete;
        PIDttmp_res(:,b) = lat_res_discrete;
        PIDttmp_shuff(:,b) = lat_sh_discrete;
    end
    PID_plugin{d_idx} = PIDttmp_plugin;
    PID_res{d_idx} = PIDttmp_res;
    PID_shuff{d_idx} = PIDttmp_shuff;
end
% Save results
save(save_path, 'PID_plugin', 'PID_res', 'PID_shuff', 's', 'condition');
fprintf('Saved: %s\n', save_path);
end