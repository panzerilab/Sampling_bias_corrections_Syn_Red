%% PID lattice + three corrections: demo
% Make sure the following are on your path:
%   pid_lattice.m
%   qele.m     % uses pid.is_gaussian internally
%   shuffsub.m     % uses pid.is_gaussian internally
%   qe.m           % uses pid.is_gaussian internally

rng(0);
clc;
clear;
addpath(genpath('/Volumes/Nicola/Data_HamePark'));
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
    pca_filename = fullfile('/Volumes/Nicola/Data_HamePark',s, sprintf('PCA_%s_%s_%s.mat', s, condition, area_name));
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
    pca_filename = fullfile('/Volumes/Nicola/Data_HamePark',s, sprintf('PCA_%s_%s_%s.mat', s, condition, area_name));
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
    pca_filename = fullfile('/Volumes/Nicola/Data_HamePark',s, sprintf('PCA_%s_%s_%s.mat', s, condition, area_name));
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

%%
clc; close all; clear;
load('subjects.mat', 'subj');
plugin_all = nan(length(subj), 18, 77);
qe_all = nan(length(subj), 18, 77);
shuff_all = nan(length(subj), 18, 77);
for subjIdx = 1:length(subj)
    path = fullfile(pwd, 'attn', sprintf('PID3Source_%s_attn.mat', subj{subjIdx}));
    res = load(path);
    plugin_all(subjIdx,:,:) = res.PID_plugin{1};
    qe_all(subjIdx,:,:) = res.PID_res{1};
    shuff_all(subjIdx,:,:) = res.PID_shuff{1};
end 

% Calculate mean and SEM, ignoring any NaN values
mean_plugin = squeeze(mean(plugin_all, 1, 'omitnan'));
sem_plugin  = squeeze(std(plugin_all, [], 1, 'omitnan') ./ sqrt(sum(~isnan(plugin_all),1)));

mean_qe = squeeze(mean(qe_all, 1, 'omitnan'));
sem_qe  = squeeze(std(qe_all, [], 1, 'omitnan') ./ sqrt(sum(~isnan(qe_all),1)));

mean_shuff  = squeeze(mean(shuff_all, 1, 'omitnan'));
sem_shuff   = squeeze(std(shuff_all, [], 1, 'omitnan') ./ sqrt(sum(~isnan(shuff_all),1)));

% Calculate joint measures, ignoring any NaN values that might have resulted
joint_mean_plugin = sum(mean_plugin(1:18,:), 1, 'omitnan');
joint_sem_plugin  = sqrt(sum(sem_plugin(1:18,:).^2, 1, 'omitnan'));  % SEM of sum

joint_mean_qe = sum(mean_qe(1:18,:), 1, 'omitnan');
joint_sem_qe  = sqrt(sum(sem_qe(1:18,:).^2, 1, 'omitnan'));

joint_mean_shuff  = sum(mean_shuff(1:18,:), 1, 'omitnan');
joint_sem_shuff   = sqrt(sum(sem_shuff(1:18,:).^2, 1, 'omitnan'));

% Extract synergy (syn) and redundancy (red) components
syn_mean_plugin = mean_plugin(18,:);
syn_sem_plugin  = sem_plugin(18,:);
syn_mean_qe = mean_qe(18,:);
syn_sem_qe  = sem_qe(18,:);
syn_mean_shuff  = mean_shuff(18,:);
syn_sem_shuff   = sem_shuff(18,:);

red_mean_plugin = mean_plugin(8,:);
red_sem_plugin  = sem_plugin(8,:);
red_mean_qe = mean_qe(8,:);
red_sem_qe  = sem_qe(8,:);
red_mean_shuff  = mean_shuff(8,:);
red_sem_shuff   = sem_shuff(8,:);

%%
fs_original = 200;  % original sampling rate (Hz)
frames = 1:311;
time_original = frames / fs_original;  % 1/200 s steps

% Downsample to 1/50 s
ds_factor = 4; % 200/50 = 4
time  = time_original(1:ds_factor:end);
time = time(1:77);


% Colors
c_plugin = [0 0.4470 0.7410];  % blue
c_qe = [0.8500 0.3250 0.0980]; % orange
c_shuff  = [0.4660 0.6740 0.1880]; % green

timesSEM =3;
alphaVal = 0.4;
figure;

% Example: alpha value for shading
alphaVal = 0.3;  % adjust for visibility

% Subplot 1: Joint
subplot(3,1,1); hold on;

% Plugin
Xfill = [time, fliplr(time)];
Yfill = [joint_mean_plugin + timesSEM*joint_sem_plugin, fliplr(joint_mean_plugin - timesSEM*joint_sem_plugin)];
fill(Xfill, Yfill, c_plugin, 'FaceAlpha', alphaVal, 'EdgeColor','none', 'HandleVisibility','off');
plot(time, joint_mean_plugin, 'Color', c_plugin, 'LineWidth', 1.5);

% Resamp
Yfill = [joint_mean_qe + timesSEM*joint_sem_qe, fliplr(joint_mean_qe - timesSEM*joint_sem_qe)];
fill(Xfill, Yfill, c_qe, 'FaceAlpha', alphaVal, 'EdgeColor','none', 'HandleVisibility','off');
plot(time, joint_mean_qe, 'Color', c_qe, 'LineWidth', 1.5);

% Shuff
Yfill = [joint_mean_shuff + timesSEM*joint_sem_shuff, fliplr(joint_mean_shuff - timesSEM*joint_sem_shuff)];
fill(Xfill, Yfill, c_shuff, 'FaceAlpha', alphaVal, 'EdgeColor','none', 'HandleVisibility','off');
plot(time, joint_mean_shuff, 'Color', c_shuff, 'LineWidth', 1.5);

title('Joint'); xlabel('Time (s)'); ylabel('Signal');
legend({'Plugin','Resamp','Shuff'});

% Subplot 2: Syn
subplot(3,1,2); hold on;

% Plugin
Yfill = [syn_mean_plugin + timesSEM*syn_sem_plugin, fliplr(syn_mean_plugin - timesSEM*syn_sem_plugin)];
fill(Xfill, Yfill, c_plugin, 'FaceAlpha', alphaVal, 'EdgeColor','none', 'HandleVisibility','off');
plot(time, syn_mean_plugin, 'Color', c_plugin, 'LineWidth', 1.5);

% Resamp
Yfill = [syn_mean_qe + timesSEM*syn_sem_qe, fliplr(syn_mean_qe - timesSEM*syn_sem_qe)];
fill(Xfill, Yfill, c_qe, 'FaceAlpha', alphaVal, 'EdgeColor','none', 'HandleVisibility','off');
plot(time, syn_mean_qe, 'Color', c_qe, 'LineWidth', 1.5);

% Shuff
Yfill = [syn_mean_shuff + timesSEM*syn_sem_shuff, fliplr(syn_mean_shuff - timesSEM*syn_sem_shuff)];
fill(Xfill, Yfill, c_shuff, 'FaceAlpha', alphaVal, 'EdgeColor','none', 'HandleVisibility','off');
plot(time, syn_mean_shuff, 'Color', c_shuff, 'LineWidth', 1.5);

title('Syn'); xlabel('Time (s)'); ylabel('Signal');

% Subplot 3: Red
subplot(3,1,3); hold on;

% Plugin
Yfill = [red_mean_plugin + timesSEM*red_sem_plugin, fliplr(red_mean_plugin - timesSEM*red_sem_plugin)];
fill(Xfill, Yfill, c_plugin, 'FaceAlpha', alphaVal, 'EdgeColor','none', 'HandleVisibility','off');
plot(time, red_mean_plugin, 'Color', c_plugin, 'LineWidth', 1.5);

% Resamp
Yfill = [red_mean_qe + timesSEM*red_sem_qe, fliplr(red_mean_qe - timesSEM*red_sem_qe)];
fill(Xfill, Yfill, c_qe, 'FaceAlpha', alphaVal, 'EdgeColor','none', 'HandleVisibility','off');
plot(time, red_mean_qe, 'Color', c_qe, 'LineWidth', 1.5);

% Shuff
Yfill = [red_mean_shuff + timesSEM*red_sem_shuff, fliplr(red_mean_shuff - timesSEM*red_sem_shuff)];
fill(Xfill, Yfill, c_shuff, 'FaceAlpha', alphaVal, 'EdgeColor','none', 'HandleVisibility','off');
plot(time, red_mean_shuff, 'Color', c_shuff, 'LineWidth', 1.5);

title('Red'); xlabel('Time (s)'); ylabel('Signal');

%legend({'Plugin','Resamp','Shuff'});
save_path = 'results_3SourcePID_discrete.png';

exportgraphics(gcf, save_path, 'Resolution', 600);  % 300 dpi



