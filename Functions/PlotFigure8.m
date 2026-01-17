clc, clear, close all;
load("subjects.mat");
PID_plugin_all = zeros(length(subj), 18, 77);
PID_shuff_all =  zeros(length(subj), 18, 77);
PID_res_all =  zeros(length(subj), 18, 77);
for subjIdx = 1:length(subj)
    subjId = subj{subjIdx};
    load(sprintf("attn_discrete/PID3Source_%s_attn.mat", subjId));
    PID_res_all(subjIdx,:,:)  = PID_res{1, 1};
    PID_plugin_all(subjIdx,:,:)  = PID_plugin{1, 1};
    PID_shuff_all(subjIdx,:,:)  = PID_shuff{1, 1};
end
idx_pre = 21;
idx_post = 36;
joint_qe_pre = sum(PID_res_all(:, 1:18, idx_pre), 2);
joint_plugin_pre = sum(PID_plugin_all(:, 1:18, idx_pre), 2);
joint_shuff_pre = sum(PID_shuff_all(:, 1:18, idx_pre), 2);
sf_down = 50;
nt = length(joint_qe_pre);
tvec = -0.55 : 1/sf_down : -0.55 + (nt - 1)/sf_down;
joint_qe_post = sum(PID_res_all(:, 1:18, idx_post), 2);
joint_plugin_post = sum(PID_plugin_all(:, 1:18, idx_post), 2);
joint_shuff_post = sum(PID_shuff_all(:, 1:18, idx_post), 2);
% Syn Pre/Post
syn_qe_pre = PID_res_all(:, 18, idx_pre);
syn_plugin_pre = PID_plugin_all(:, 18, idx_pre);
syn_shuff_pre = PID_shuff_all(:, 18, idx_pre);
syn_qe_post = PID_res_all(:, 18, idx_post);
syn_plugin_post = PID_plugin_all(:, 18, idx_post);
syn_shuff_post = PID_shuff_all(:, 18, idx_post);
% Red Pre/Post
red_qe_pre = PID_res_all(:, 8, idx_pre);
red_plugin_pre = PID_plugin_all(:, 8, idx_pre);
red_shuff_pre = PID_shuff_all(:, 8, idx_pre);
red_qe_post = PID_res_all(:, 8, idx_post);
red_plugin_post = PID_plugin_all(:, 8, idx_post);
red_shuff_post = PID_shuff_all(:, 8, idx_post);
% --- Statistical calculations (Mean and SEM) ---
% Note: MATLAB's standard deviation is divided by N-1, not N.
% To calculate the SEM, we use std(data) / sqrt(length(data)).
function [means, sems] = calculate_stats(data_plugin, data_shuff, data_qe)
    calc_sem = @(data) std(data) / sqrt(length(data));
    means = [mean(data_plugin), mean(data_shuff), mean(data_qe)];
    sems = [calc_sem(data_plugin), calc_sem(data_shuff), calc_sem(data_qe)];
end
% Pre-stimulus data
[pre_means_joint, pre_sems_joint] = calculate_stats(joint_plugin_pre, joint_shuff_pre, joint_qe_pre);
[pre_means_syn, pre_sems_syn] = calculate_stats(syn_plugin_pre, syn_shuff_pre, syn_qe_pre);
[pre_means_red, pre_sems_red] = calculate_stats(red_plugin_pre, red_shuff_pre, red_qe_pre);
% Post-stimulus data
[post_means_joint, post_sems_joint] = calculate_stats(joint_plugin_post, joint_shuff_post, joint_qe_post);
[post_means_syn, post_sems_syn] = calculate_stats(syn_plugin_post, syn_shuff_post, syn_qe_post);
[post_means_red, post_sems_red] = calculate_stats(red_plugin_post, red_shuff_post, red_qe_post);
% Second data set from attn_gauss
PID_plugin_all_gauss = zeros(length(subj), 18, 77);
PID_shuff_all_gauss =  zeros(length(subj), 18, 77);
PID_res_all_gauss =  zeros(length(subj), 18, 77);
for subjIdx = 1:length(subj)
    subjId = subj{subjIdx};
    load(sprintf("attn_gauss/PID3Source_%s_attn.mat", subjId));
    PID_res_all_gauss(subjIdx,:,:)  = PID_res{1, 1};
    PID_plugin_all_gauss(subjIdx,:,:)  = PID_plugin{1, 1};
    PID_shuff_all_gauss(subjIdx,:,:)  = PID_shuff{1, 1};
end
joint_qe_pre_gauss = sum(PID_res_all_gauss(:, 1:18, idx_pre), 2);
joint_plugin_pre_gauss = sum(PID_plugin_all_gauss(:, 1:18, idx_pre), 2);
joint_shuff_pre_gauss = sum(PID_shuff_all_gauss(:, 1:18, idx_pre), 2);
joint_qe_post_gauss = sum(PID_res_all_gauss(:, 1:18, idx_post), 2);
joint_plugin_post_gauss = sum(PID_plugin_all_gauss(:, 1:18, idx_post), 2);
joint_shuff_post_gauss = sum(PID_shuff_all_gauss(:, 1:18, idx_post), 2);
% Syn Pre/Post
syn_qe_pre_gauss = PID_res_all_gauss(:, 18, idx_pre);
syn_plugin_pre_gauss = PID_plugin_all_gauss(:, 18, idx_pre);
syn_shuff_pre_gauss = PID_shuff_all_gauss(:, 18, idx_pre);
syn_qe_post_gauss = PID_res_all_gauss(:, 18, idx_post);
syn_plugin_post_gauss = PID_plugin_all_gauss(:, 18, idx_post);
syn_shuff_post_gauss = PID_shuff_all_gauss(:, 18, idx_post);
% Red Pre/Post
red_qe_pre_gauss = PID_res_all_gauss(:, 8, idx_pre);
red_plugin_pre_gauss = PID_plugin_all_gauss(:, 8, idx_pre);
red_shuff_pre_gauss = PID_shuff_all_gauss(:, 8, idx_pre);
red_qe_post_gauss = PID_res_all_gauss(:, 8, idx_post);
red_plugin_post_gauss = PID_plugin_all_gauss(:, 8, idx_post);
red_shuff_post_gauss = PID_shuff_all_gauss(:, 8, idx_post);
% Statistical calculations (Mean and SEM) for Gaussian data
[pre_means_joint_gauss, pre_sems_joint_gauss] = calculate_stats(joint_plugin_pre_gauss, joint_shuff_pre_gauss, joint_qe_pre_gauss);
[pre_means_syn_gauss, pre_sems_syn_gauss] = calculate_stats(syn_plugin_pre_gauss, syn_shuff_pre_gauss, syn_qe_pre_gauss);
[pre_means_red_gauss, pre_sems_red_gauss] = calculate_stats(red_plugin_pre_gauss, red_shuff_pre_gauss, red_qe_pre_gauss);
[post_means_joint_gauss, post_sems_joint_gauss] = calculate_stats(joint_plugin_post_gauss, joint_shuff_post_gauss, joint_qe_post_gauss);
[post_means_syn_gauss, post_sems_syn_gauss] = calculate_stats(syn_plugin_post_gauss, syn_shuff_post_gauss, syn_qe_post_gauss);
[post_means_red_gauss, post_sems_red_gauss] = calculate_stats(red_plugin_post_gauss, red_shuff_post_gauss, red_qe_post_gauss);
% --- Plotting ---
figure('Position', [100, 100, 1200, 600]); % Taller figure for two rows
set(0, 'DefaultAxesFontSize', 8);
% Use new, visually distinct colors
colors = [
    0.5843, 0.7333, 0.8431; % Light Blue
    0.8980, 0.6902, 0.5059; % Light Orange
    0.5804, 0.8118, 0.6941  % Light Green
];
tolNeg = 0.005;
tolPos = 0.05;
method_names = {'Plugin', 'Shuff', 'QE'};
% Determine consistent y-axis limits for each condition
discrete_pre_min = min([pre_means_joint, pre_means_syn, pre_means_red]) - tolNeg;
discrete_pre_max = max([pre_means_joint, pre_means_syn, pre_means_red]) + tolPos;
discrete_post_min = min([post_means_joint, post_means_syn, post_means_red]) - tolNeg;
discrete_post_max = max([post_means_joint, post_means_syn, post_means_red]) + tolPos;
gauss_pre_min = min([pre_means_joint_gauss, pre_means_syn_gauss, pre_means_red_gauss]) - tolNeg;
gauss_pre_max = max([pre_means_joint_gauss, pre_means_syn_gauss, pre_means_red_gauss]) + tolPos;
gauss_post_min = min([post_means_joint_gauss, post_means_syn_gauss, post_means_red_gauss]) - tolNeg;
gauss_post_max = max([post_means_joint_gauss, post_means_syn_gauss, post_means_red_gauss]) + tolPos;
sgtitle('Discrete Data', 'FontSize', 12, 'FontWeight', 'bold');
% First Row (Discrete Data)
subplot(2, 6, 1);
bar_pre_joint = bar(pre_means_joint, 'FaceColor', 'flat');
bar_pre_joint.CData = colors;
hold on;
er = errorbar(1:3, pre_means_joint, 2 * pre_sems_joint);
er.Color = [0 0 0];
er.LineStyle = 'none';
er.CapSize = 5;
hold off;
title('Joint');
set(gca, 'Box', 'off');
set(gca, 'XColor', 'k', 'YColor', 'k');
set(gca, 'XTick', 1:3, 'XTickLabel', method_names);
ylabel('Info [bits]');
ylim([discrete_pre_min, discrete_pre_max]);
subplot(2, 6, 2);
bar_pre_syn = bar(pre_means_syn, 'FaceColor', 'flat');
bar_pre_syn.CData = colors;
hold on;
er = errorbar(1:3, pre_means_syn, 2 * pre_sems_syn);
er.Color = [0 0 0];
er.LineStyle = 'none';
er.CapSize = 5;
hold off;
title('Syn');
set(gca, 'Box', 'off');
set(gca, 'XColor', 'k', 'YColor', 'k');
set(gca, 'XTick', 1:3, 'XTickLabel', method_names);
set(gca, 'YTickLabel', []);
ylim([discrete_pre_min, discrete_pre_max]);
subplot(2, 6, 3);
bar_pre_red = bar(pre_means_red, 'FaceColor', 'flat');
bar_pre_red.CData = colors;
hold on;
er = errorbar(1:3, pre_means_red, 2 * pre_sems_red);
er.Color = [0 0 0];
er.LineStyle = 'none';
er.CapSize = 5;
hold off;
title('Red');
set(gca, 'Box', 'off');
set(gca, 'XColor', 'k', 'YColor', 'k');
set(gca, 'XTick', 1:3, 'XTickLabel', method_names);
set(gca, 'YTickLabel', []);
ylim([discrete_pre_min, discrete_pre_max]);
subplot(2, 6, 4);
bar_post_joint = bar(post_means_joint, 'FaceColor', 'flat');
bar_post_joint.CData = colors;
hold on;
er = errorbar(1:3, post_means_joint, 2 * post_sems_joint);
er.Color = [0 0 0];
er.LineStyle = 'none';
er.CapSize = 5;
hold off;
title('Joint');
set(gca, 'Box', 'off');
set(gca, 'XColor', 'k', 'YColor', 'k');
set(gca, 'XTick', 1:3, 'XTickLabel', method_names);
ylabel('Info [bits]');
ylim([discrete_post_min, discrete_post_max]);
subplot(2, 6, 5);
bar_post_syn = bar(post_means_syn, 'FaceColor', 'flat');
bar_post_syn.CData = colors;
hold on;
er = errorbar(1:3, post_means_syn, 2 * post_sems_syn);
er.Color = [0 0 0];
er.LineStyle = 'none';
er.CapSize = 5;
hold off;
title('Syn');
set(gca, 'Box', 'off');
set(gca, 'XColor', 'k', 'YColor', 'k');
set(gca, 'XTick', 1:3, 'XTickLabel', method_names);
set(gca, 'YTickLabel', []);
ylim([discrete_post_min, discrete_post_max]);
subplot(2, 6, 6);
bar_post_red = bar(post_means_red, 'FaceColor', 'flat');
bar_post_red.CData = colors;
hold on;
er = errorbar(1:3, post_means_red, 2 * post_sems_red);
er.Color = [0 0 0];
er.LineStyle = 'none';
er.CapSize = 5;
hold off;
title('Red');
set(gca, 'Box', 'off');
set(gca, 'XColor', 'k', 'YColor', 'k');
set(gca, 'XTick', 1:3, 'XTickLabel', method_names);
set(gca, 'YTickLabel', []);
ylim([discrete_post_min, discrete_post_max]);
sgtitle('Gaussian Data', 'FontSize', 12, 'FontWeight', 'bold');
% Second Row (Gaussian Data)
subplot(2, 6, 7);
bar_pre_joint_gauss = bar(pre_means_joint_gauss, 'FaceColor', 'flat');
bar_pre_joint_gauss.CData = colors;
hold on;
er = errorbar(1:3, pre_means_joint_gauss, 2 * pre_sems_joint_gauss);
er.Color = [0 0 0];
er.LineStyle = 'none';
er.CapSize = 5;
hold off;
title('Joint');
set(gca, 'Box', 'off');
set(gca, 'XColor', 'k', 'YColor', 'k');
set(gca, 'XTick', 1:3, 'XTickLabel', method_names);
ylabel('Info [bits]');
ylim([gauss_pre_min, gauss_pre_max]);
subplot(2, 6, 8);
bar_pre_syn_gauss = bar(pre_means_syn_gauss, 'FaceColor', 'flat');
bar_pre_syn_gauss.CData = colors;
hold on;
er = errorbar(1:3, pre_means_syn_gauss, 2 * pre_sems_syn_gauss);
er.Color = [0 0 0];
er.LineStyle = 'none';
er.CapSize = 5;
hold off;
title('Syn');
set(gca, 'Box', 'off');
set(gca, 'XColor', 'k', 'YColor', 'k');
set(gca, 'XTick', 1:3, 'XTickLabel', method_names);
set(gca, 'YTickLabel', []);
ylim([gauss_pre_min, gauss_pre_max]);
subplot(2, 6, 9);
bar_pre_red_gauss = bar(pre_means_red_gauss, 'FaceColor', 'flat');
bar_pre_red_gauss.CData = colors;
hold on;
er = errorbar(1:3, pre_means_red_gauss, 2 * pre_sems_red_gauss);
er.Color = [0 0 0];
er.LineStyle = 'none';
er.CapSize = 5;
hold off;
title('Red');
set(gca, 'Box', 'off');
set(gca, 'XColor', 'k', 'YColor', 'k');
set(gca, 'XTick', 1:3, 'XTickLabel', method_names);
set(gca, 'YTickLabel', []);
ylim([gauss_pre_min, gauss_pre_max]);
subplot(2, 6, 10);
bar_post_joint_gauss = bar(post_means_joint_gauss, 'FaceColor', 'flat');
bar_post_joint_gauss.CData = colors;
hold on;
er = errorbar(1:3, post_means_joint_gauss, 2 * post_sems_joint_gauss);
er.Color = [0 0 0];
er.LineStyle = 'none';
er.CapSize = 5;
hold off;
title('Joint');
set(gca, 'Box', 'off');
set(gca, 'XColor', 'k', 'YColor', 'k');
set(gca, 'XTick', 1:3, 'XTickLabel', method_names);
ylabel('Info [bits]');
ylim([gauss_post_min, gauss_post_max]);
subplot(2, 6, 11);
bar_post_syn_gauss = bar(post_means_syn_gauss, 'FaceColor', 'flat');
bar_post_syn_gauss.CData = colors;
hold on;
er = errorbar(1:3, post_means_syn_gauss, 2 * post_sems_syn_gauss);
er.Color = [0 0 0];
er.LineStyle = 'none';
er.CapSize = 5;
hold off;
title('Syn');
set(gca, 'Box', 'off');
set(gca, 'XColor', 'k', 'YColor', 'k');
set(gca, 'XTick', 1:3, 'XTickLabel', method_names);
set(gca, 'YTickLabel', []);
ylim([gauss_post_min, gauss_post_max]);
subplot(2, 6, 12);
bar_post_red_gauss = bar(post_means_red_gauss, 'FaceColor', 'flat');
bar_post_red_gauss.CData = colors;
hold on;
er = errorbar(1:3, post_means_red_gauss, 2 * post_sems_red_gauss);
er.Color = [0 0 0];
er.LineStyle = 'none';
er.CapSize = 5;
hold off;
title('Red');
set(gca, 'Box', 'off');
set(gca, 'XColor', 'k', 'YColor', 'k');
set(gca, 'XTick', 1:3, 'XTickLabel', method_names);
set(gca, 'YTickLabel', []);
ylim([gauss_post_min, gauss_post_max]);


outDir = fullfile(fileparts(pwd), 'Figures_mat'); 
if ~exist(outDir, 'dir')
    mkdir(outDir);
end
exportgraphics(gcf, fullfile(outDir, 'Figure_8.svg'), 'ContentType','vector');