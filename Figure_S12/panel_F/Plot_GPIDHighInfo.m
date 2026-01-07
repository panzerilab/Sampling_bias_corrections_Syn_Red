clc, clear, close all;

IGt = 1.628989209509945;
mean_Joint = [];
mean_RED =[];
mean_UNQ1 = [];
mean_UNQ2 =[];
mean_SYN =[];

load('Sim_results/data_trials_16.mat')
mean_Joint(end+1) = mean(Ijoint_vals);
mean_RED(end+1)= mean(red_venk_vals);
mean_UNQ1(end+1) = mean(unique_info1_vals);
mean_UNQ2(end+1) = mean(unique_info2_vals);
mean_SYN(end+1) = mean(syn_venk_vals);

load('Sim_results/data_trials_32.mat')
mean_Joint(end+1) = mean(Ijoint_vals);
mean_RED(end+1)= mean(red_venk_vals);
mean_UNQ1(end+1) = mean(unique_info1_vals);
mean_UNQ2(end+1) = mean(unique_info2_vals);
mean_SYN(end+1) = mean(syn_venk_vals);
load('Sim_results/data_trials_64.mat')
mean_Joint(end+1) = mean(Ijoint_vals);
mean_RED(end+1)= mean(red_venk_vals);
mean_UNQ1(end+1) = mean(unique_info1_vals);
mean_UNQ2(end+1) = mean(unique_info2_vals);
mean_SYN(end+1) = mean(syn_venk_vals);
load('Sim_results/data_trials_128.mat')
mean_Joint(end+1) = mean(Ijoint_vals);
mean_RED(end+1)= mean(red_venk_vals);
mean_UNQ1(end+1) = mean(unique_info1_vals);
mean_UNQ2(end+1) = mean(unique_info2_vals);
mean_SYN(end+1) = mean(syn_venk_vals);
load('Sim_results/data_trials_256.mat')
mean_Joint(end+1) = mean(Ijoint_vals);
mean_RED(end+1)= mean(red_venk_vals);
mean_UNQ1(end+1) = mean(unique_info1_vals);
mean_UNQ2(end+1) = mean(unique_info2_vals);
mean_SYN(end+1) = mean(syn_venk_vals);
load('Sim_results/data_trials_512.mat')
mean_Joint(end+1) = mean(Ijoint_vals);
mean_RED(end+1)= mean(red_venk_vals);
mean_UNQ1(end+1) = mean(unique_info1_vals);
mean_UNQ2(end+1) = mean(unique_info2_vals);
mean_SYN(end+1) = mean(syn_venk_vals);
load('Sim_results/data_trials_1024.mat')
mean_Joint(end+1) = mean(Ijoint_vals);
mean_RED(end+1)= mean(red_venk_vals);
mean_UNQ1(end+1) = mean(unique_info1_vals);
mean_UNQ2(end+1) = mean(unique_info2_vals);
mean_SYN(end+1) = mean(syn_venk_vals);
load('Sim_results/data_trials_2048.mat')
mean_Joint(end+1) = mean(Ijoint_vals);
mean_RED(end+1)= mean(red_venk_vals);
mean_UNQ1(end+1) = mean(unique_info1_vals);
mean_UNQ2(end+1) = mean(unique_info2_vals);
mean_SYN(end+1) = mean(syn_venk_vals);

%%
set(0, 'DefaultTextFontName', 'Arial');
set(0, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultUicontrolFontName', 'Arial');
set(0, 'DefaultUitableFontName', 'Arial');
set(0, 'DefaultUipanelFontName', 'Arial');
set(0, 'DefaultLegendFontName', 'Arial');

customColors_info_PID = [
    0.9290, 0.6940, 0.1250;  % Color for 'Joint'
    0, 0.4470, 0.7410;       % Color for 'Red'
    0.8500, 0.3250, 0.0980;  % Color for 'Union'
    0.4940, 0.1840, 0.5560;  % Color for 'SR'
    0.4660, 0.6740, 0.1880;  % Color for 'Syn'
    0.5, 0.5, 0.5;           % Color for 'U1'
    0.4, 0.26, 0.13;         % Color for 'U2'
    0.4940, 0.1840, 0.5560];  % Color for 'PInd'


trial_categories =  [16,32,64,128,256,512, 1024, 2048];
axisLineWidth = 2; 
linewidth = 2;

figure_handle = figure('Units', 'centimeters', 'Position', [1, 1, 10, 10]);
hold on;
plot(trial_categories, mean_Joint, 'Color', customColors_info_PID(1,:), 'LineWidth', linewidth);
plot(trial_categories, mean_RED, 'Color', customColors_info_PID(2,:), 'LineWidth',linewidth);
plot(trial_categories, mean_UNQ1, 'Color', customColors_info_PID(6,:), 'LineWidth', linewidth);
plot(trial_categories, mean_UNQ2, 'Color', customColors_info_PID(3,:), 'LineWidth', linewidth);
plot(trial_categories, mean_SYN, 'Color', customColors_info_PID(5,:), 'LineWidth', linewidth);
plot(trial_categories, IGt * ones(size(trial_categories)), 'k--','Color', customColors_info_PID(1,:), 'LineWidth', linewidth, 'DisplayName', 'JointInf GroundTruth');

yline(0, 'LineWidth',1.2);


tickFontSize = 18;
FontSize = tickFontSize;
legendFontSize = tickFontSize;
titleFontWeight = 'bold';
set(gca, 'XTickLabel', get(gca, 'XTickLabel'), 'FontSize', tickFontSize);
set(gca, 'XScale', 'log', 'LineWidth',1.2);
set(gca, 'FontSize', tickFontSize, 'LineWidth', axisLineWidth);
xticks([100, 1000]);
xticklabels({'10^2', '10^3'});
xlim([16 2048]);
ylim([0 1.7]);
lgd = legend({'Joint', 'Red', 'U1', 'U2', 'Syn'}, 'Box', 'off', 'FontSize', legendFontSize,'Location', 'northeast');
 xlabel('Trials per stimulus', 'FontSize', FontSize);
  ylabel('Information [bits]', 'FontSize', FontSize);

  saveas(figure_handle, 'GPID_highInf.svg');
  