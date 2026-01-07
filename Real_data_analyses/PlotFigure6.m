clc, clear, close all;

set(0, 'DefaultTextFontName', 'Arial');
set(0, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultUicontrolFontName', 'Arial');
set(0, 'DefaultUitableFontName', 'Arial');
set(0, 'DefaultUipanelFontName', 'Arial');
set(0, 'DefaultLegendFontName', 'Arial');

customColors_info_PID = [
    0.9290, 0.6940, 0.1250;  % Color for 'Joint'
    0, 0.4470, 0.7410;       % Color for 'Red'
    0.4660, 0.6740, 0.1880;  % Color for 'Syn'
    0.5, 0.5, 0.5;           % Color for 'Unq'
];

% Joint, Syn, Red, Unq
chosenAtom = 'Red';
%% ---------- Adaptive weighting helpers (QE ⟷ shuff-sub) ----------
% Load once (file used in your first function)
load('Functions/adaptive_weight_matrix_allcases.mat', 'weight_matrix', 'info_levels_fine');

% Build a 1D "weight curve" vs. info fraction by averaging across the
% trial-count axis of your 3D weight tensor. The first index (=3) matches
% how you selected weights for Red in your first function.
weight_curve_vs_info = squeeze(mean(weight_matrix(3,:,:), 2));  % [nInfoLevels x 1]

% Simple accessor: map info fraction in [0,1] → weight w in [0,1]
get_weight = @(info_frac) interp1(info_levels_fine(:), weight_curve_vs_info(:), ...
                                  max(0,min(1,info_frac)), 'linear', 'extrap');
%% Prepare A1 data
filename = ['Results/Real_data_analysis/A1_PID_' chosenAtom '.mat'];
load(filename);

plugin_all      = [];
qe_all          = [];
qeshuffSub_all  = [];
shuffSub_all    = [];
for dataSet = 1:length(Results_plugin)
    if ~isempty(Results_plugin(dataSet).FullData)
        plugin_all     = [plugin_all    ; Results_plugin(dataSet).FullData.Result];
        qe_all         = [qe_all        ; Results_qe(dataSet).FullData.Result];
        qeshuffSub_all = [qeshuffSub_all; Results_qeshuffSub(dataSet).FullData.Result];
        shuffSub_all   = [shuffSub_all  ; Results_shuffSub(dataSet).FullData.Result];
    end
end 

[h, p_A1_plugin, ci, stats] = ttest(plugin_all(:,3), plugin_all(:,2));
[h, p_A1_qe, ci, stats] = ttest(qe_all(:,3), qe_all(:,2));
[h, p_A1_shuffSub, ci, stats] = ttest(shuffSub_all(:,3), shuffSub_all(:,2));
[h, p_A1_qeshuffSub, ci, stats] = ttest(qeshuffSub_all(:,3), qeshuffSub_all(:,2));


% Calculate joint means and SEM (column 5)
meanJointA1 = [mean(plugin_all(:,1)),  mean(qe_all(:,1)), ...
               mean(shuffSub_all(:,1)), mean(qeshuffSub_all(:,1))];
semJointA1  = [std(plugin_all(:,1))/sqrt(size(plugin_all,1)), ...
               std(qe_all(:,1))/sqrt(size(qe_all,1)), ...
               std(shuffSub_all(:,1))/sqrt(size(shuffSub_all,1)), ...
               std(qeshuffSub_all(:,1))/sqrt(size(qeshuffSub_all,1))];
% ow 
meanpluginA1 = [mean(plugin_all(:,3)), ...
                mean(plugin_all(:,2))];
sempluginA1  =  [std(plugin_all(:,3))/sqrt(size(plugin_all,1)), ...
                 std(plugin_all(:,2))/sqrt(size(plugin_all,1))];

meanqeA1 = [mean(qe_all(:,3)), ...
            mean(qe_all(:,2))];
semqeA1  =  [std(qe_all(:,3))/sqrt(size(qe_all,1)), ...
             std(qe_all(:,2))/sqrt(size(qe_all,1))];

meanshuffsubA1 = [mean(shuffSub_all(:,3)), ...
                  mean(shuffSub_all(:,2))];
semshuffsubA1  =  [std(shuffSub_all(:,3))/sqrt(size(shuffSub_all,1)), ...
                   std(shuffSub_all(:,2))/sqrt(size(shuffSub_all,1))];


meanqeShuffsubA1 = [mean(qeshuffSub_all(:,3)), ...
                    mean(qeshuffSub_all(:,2))];
semqeShuffsubA1  =  [std(qeshuffSub_all(:,3))/sqrt(size(qeshuffSub_all,1)), ...
                     std(qeshuffSub_all(:,2))/sqrt(size(qeshuffSub_all,1))];

% Adaptive weighted blend (replaces qeshuffSub)
max_info = log2(2);  % same as in your first function
[weighted_all, wA1] = apply_adaptive_weighting(qe_all, shuffSub_all, max_info, get_weight);

% Stats using the same columns you t-test (Red vs Syn are cols 3 vs 2 in your code)
[~, p_A1_weighted] = ttest(weighted_all(:,3), weighted_all(:,2));

% Joint means/SEMs (col 1) for bars
meanJointA1_weighted = mean(weighted_all(:,1), 'omitnan');
semJointA1_weighted  = std(weighted_all(:,1),  'omitnan')/sqrt(size(weighted_all,1));

% Component bars (Red, Syn) for the two-bar subplots
meanWeightedA1 = [mean(weighted_all(:,3),'omitnan'), mean(weighted_all(:,2),'omitnan')];
semWeightedA1  = [std(weighted_all(:,3),'omitnan')/sqrt(size(weighted_all,1)), ...
                  std(weighted_all(:,2),'omitnan')/sqrt(size(weighted_all,1))];



%% Prepare CA1 data
filename = ['Results/Real_data_analysis/CA1_PID_' chosenAtom '.mat'];
load(filename);
plugin_all      = [];
qe_all          = [];
qeshuffSub_all  = [];
shuffSub_all    = [];
for dataSet = 1:length(Results_plugin)
    if ~isempty(Results_plugin(dataSet).FullData)
        plugin_all     = [plugin_all    ; Results_plugin(dataSet).FullData.Result];
        qe_all         = [qe_all        ; Results_qe(dataSet).FullData.Result];
        qeshuffSub_all = [qeshuffSub_all; Results_qeshuffSub(dataSet).FullData.Result];
        shuffSub_all   = [shuffSub_all  ; Results_shuffSub(dataSet).FullData.Result];
    end
end 
[h, p_CA1_plugin, ci, stats] = ttest(plugin_all(:,3), plugin_all(:,2));
[h, p_CA1_qe, ci, stats] = ttest(qe_all(:,3), qe_all(:,2));
[h, p_CA1_shuffSub, ci, stats] = ttest(shuffSub_all(:,3), shuffSub_all(:,2));
[h, p_CA1_qeshuffSub, ci, stats] = ttest(qeshuffSub_all(:,3), qeshuffSub_all(:,2));

% Calculate joint means and SEM (column 5)
meanJointCA1 = [mean(plugin_all(:,1)),  mean(qe_all(:,1)), ...
               mean(shuffSub_all(:,1)), mean(qeshuffSub_all(:,1))];
semJointCA1  = [std(plugin_all(:,1))/sqrt(size(plugin_all,1)), ...
               std(qe_all(:,1))/sqrt(size(qe_all,1)), ...
               std(shuffSub_all(:,1))/sqrt(size(shuffSub_all,1)), ...
               std(qeshuffSub_all(:,1))/sqrt(size(qeshuffSub_all,1))];

meanpluginCA1 = [mean(plugin_all(:,3)), ...
                mean(plugin_all(:,2))];
sempluginCA1  =  [std(plugin_all(:,3))/sqrt(size(plugin_all,1)), ...
                 std(plugin_all(:,2))/sqrt(size(plugin_all,1))];

meanqeCA1 = [mean(qe_all(:,3)), ...
            mean(qe_all(:,2))];
semqeCA1  =  [std(qe_all(:,3))/sqrt(size(qe_all,1)), ...
             std(qe_all(:,2))/sqrt(size(qe_all,1))];

meanshuffsubCA1 = [mean(shuffSub_all(:,3)), ...
                  mean(shuffSub_all(:,2))];
semshuffsubCA1  =  [std(shuffSub_all(:,3))/sqrt(size(shuffSub_all,1)), ...
                   std(shuffSub_all(:,2))/sqrt(size(shuffSub_all,1))];


meanqeShuffsubCA1 = [mean(qeshuffSub_all(:,3)), ...
                    mean(qeshuffSub_all(:,2))];
semqeShuffsubCA1  =  [std(qeshuffSub_all(:,3))/sqrt(size(qeshuffSub_all,1)), ...
                     std(qeshuffSub_all(:,2))/sqrt(size(qeshuffSub_all,1))];

max_info = log2(4);
[weighted_all, wCA1] = apply_adaptive_weighting(qe_all, shuffSub_all, max_info, get_weight);
[~, p_CA1_weighted] = ttest(weighted_all(:,3), weighted_all(:,2));
meanJointCA1_weighted = mean(weighted_all(:,1),'omitnan');
semJointCA1_weighted  = std(weighted_all(:,1),'omitnan')/sqrt(size(weighted_all,1));
meanWeightedCA1 = [mean(weighted_all(:,3),'omitnan'), mean(weighted_all(:,2),'omitnan')];
semWeightedCA1  = [std(weighted_all(:,3),'omitnan')/sqrt(size(weighted_all,1)), ...
                   std(weighted_all(:,2),'omitnan')/sqrt(size(weighted_all,1))];


%% Prepare A1_monkey data
filename = ['Results/Real_data_analysis/KayserA1_PID_' chosenAtom '.mat'];
load(filename);
plugin_all      = [];
qe_all          = [];
qeshuffSub_all  = [];
shuffSub_all    = [];
for dataSet = 1:length(Results_plugin)
    if ~isempty(Results_plugin(dataSet).FullData)
        plugin_all     = [plugin_all    ; Results_plugin(dataSet).FullData.Result];
        qe_all         = [qe_all        ; Results_qe(dataSet).FullData.Result];
        qeshuffSub_all = [qeshuffSub_all; Results_qeshuffSub(dataSet).FullData.Result];
        shuffSub_all   = [shuffSub_all  ; Results_shuffSub(dataSet).FullData.Result];
    end
end 
[h, p_A1_monkey_plugin, ci, stats] = ttest(plugin_all(:,3), plugin_all(:,2));
[h, p_A1_monkey_qe, ci, stats] = ttest(qe_all(:,3), qe_all(:,2));
[h, p_A1_monkey_shuffSub, ci, stats] = ttest(shuffSub_all(:,3), shuffSub_all(:,2));
[h, p_A1_monkey_qeshuffSub, ci, stats] = ttest(qeshuffSub_all(:,3), qeshuffSub_all(:,2));

[h, p_A1_qe1, ci, stats] = ttest(qe_all(:,2), 0);
[h, p_A1_qe2, ci, stats] = ttest(qe_all(:,3), 0);

[h, p_A1_1, ci, stats] = ttest(shuffSub_all(:,2), 0);
[h, p_A1_2, ci, stats] = ttest(shuffSub_all(:,3), 0);

% Calculate joint means and SEM (column 5)
meanJointA1_monkey = [mean(plugin_all(:,1)),  mean(qe_all(:,1)), ...
               mean(shuffSub_all(:,1)), mean(qeshuffSub_all(:,1))];
semJointA1_monkey  = [std(plugin_all(:,1))/sqrt(size(plugin_all,1)), ...
               std(qe_all(:,1))/sqrt(size(qe_all,1)), ...
               std(shuffSub_all(:,1))/sqrt(size(shuffSub_all,1)), ...
               std(qeshuffSub_all(:,1))/sqrt(size(qeshuffSub_all,1))];

meanpluginA1_monkey = [mean(plugin_all(:,3)), ...
                mean(plugin_all(:,2))];
sempluginA1_monkey  =  [std(plugin_all(:,3))/sqrt(size(plugin_all,1)), ...
                 std(plugin_all(:,2))/sqrt(size(plugin_all,1))];

meanqeA1_monkey = [mean(qe_all(:,3)), ...
            mean(qe_all(:,2))];
semqeA1_monkey  =  [std(qe_all(:,3))/sqrt(size(qe_all,1)), ...
             std(qe_all(:,2))/sqrt(size(qe_all,1))];

meanshuffsubA1_monkey = [mean(shuffSub_all(:,3)), ...
                  mean(shuffSub_all(:,2))];
semshuffsubA1_monkey  =  [std(shuffSub_all(:,3))/sqrt(size(shuffSub_all,1)), ...
                   std(shuffSub_all(:,2))/sqrt(size(shuffSub_all,1))];


meanqeShuffsubA1_monkey = [mean(qeshuffSub_all(:,3)), ...
                    mean(qeshuffSub_all(:,2))];
semqeShuffsubA1_monkey  =  [std(qeshuffSub_all(:,3))/sqrt(size(qeshuffSub_all,1)), ...
                     std(qeshuffSub_all(:,2))/sqrt(size(qeshuffSub_all,1))];


% Adaptive weighted blend (replaces qeshuffSub)
max_info = log2(4);  % same as in your first function
[weighted_all, wA1M] = apply_adaptive_weighting(qe_all, shuffSub_all, max_info, get_weight);
[~, p_A1_monkey_weighted] = ttest(weighted_all(:,3), weighted_all(:,2));
meanJointA1_monkey_weighted = mean(weighted_all(:,1),'omitnan');
semJointA1_monkey_weighted  = std(weighted_all(:,1),'omitnan')/sqrt(size(weighted_all,1));
meanWeightedA1_monkey = [mean(weighted_all(:,3),'omitnan'), mean(weighted_all(:,2),'omitnan')];
semWeightedA1_monkey  = [std(weighted_all(:,3),'omitnan')/sqrt(size(weighted_all,1)), ...
                         std(weighted_all(:,2),'omitnan')/sqrt(size(weighted_all,1))];


%% Prepare Gaussian data
filename = ['Results/Real_data_analysis/Gaussian_data.mat'];
load(filename);
plugin_all      = [];
resample_all    = [];
resample_shuff_all   = [];
shuffSub_all    = [];

plugin_all = reshape(squeeze(results(:, :, 20, 1, :)), [], 4);
resample_all = reshape(squeeze(results(:, :, 20, 2, :)), [], 4);
shuffSub_all = reshape(squeeze(results(:, :, 20, 3, :)), [], 4);
resample_shuff_all = cat(3, resample_all, shuffSub_all);
resample_shuff_all = squeeze(mean(resample_shuff_all,3));

[h, p_gaussian_plugin, ci, stats] = ttest(plugin_all(:,1), plugin_all(:,2));
[h, p_gaussian_resample, ci, stats] = ttest(resample_all(:,1), resample_all(:,2));
[h, p_gaussian_shuffle, ci, stats] = ttest(shuffSub_all(:,1), shuffSub_all(:,2));
[h, p_gaussian_resample_shuff, ci, stats] = ttest(resample_shuff_all(:,1), resample_shuff_all(:,2));


% Calculate joint means and SEM (column 5)
meanJoint_gaussian = [mean(plugin_all(:,4), 'omitnan'),  mean(resample_all(:,4), 'omitnan'), ...
               mean(shuffSub_all(:,4), 'omitnan'), mean(resample_shuff_all(:,4), 'omitnan')];
semJoint_gaussian  = [std(plugin_all(:,4), 'omitnan')/sqrt(sum(~isnan(plugin_all(:,4)))), ...
                      std(resample_all(:,4), 'omitnan')/sqrt(sum(~isnan(resample_all(:,4)))), ...
                      std(shuffSub_all(:,4), 'omitnan')/sqrt(sum(~isnan(shuffSub_all(:,4)))), ...
                      std(resample_shuff_all(:,4), 'omitnan')/sqrt(sum(~isnan(resample_shuff_all(:,4))))];

meanplugin_gaussian = [mean(plugin_all(:,2), 'omitnan'), ...
                       mean(plugin_all(:,1), 'omitnan')];

meanresample_gaussian = [mean(resample_all(:,2), 'omitnan'), ...
                   mean(resample_all(:,1), 'omitnan')];

meanshuffsub_gaussian = [mean(shuffSub_all(:,2), 'omitnan'), ...
                         mean(shuffSub_all(:,1), 'omitnan')];

meanresample_shuff_gaussian = [mean(resample_shuff_all(:,2), 'omitnan'), ...
                           mean(resample_shuff_all(:,1), 'omitnan')];


semplugin_gaussian  = [std(plugin_all(:,2), 'omitnan') / sqrt(sum(~isnan(plugin_all(:,2)))), ...
                       std(plugin_all(:,1), 'omitnan') / sqrt(sum(~isnan(plugin_all(:,1))))];

semresample_gaussian  = [std(resample_all(:,2), 'omitnan') / sqrt(sum(~isnan(resample_all(:,2)))), ...
                         std(resample_all(:,1), 'omitnan') / sqrt(sum(~isnan(resample_all(:,1))))];

semshuffsub_gaussian  = [std(shuffSub_all(:,3), 'omitnan') / sqrt(sum(~isnan(shuffSub_all(:,2)))), ...
                         std(shuffSub_all(:,2), 'omitnan') / sqrt(sum(~isnan(shuffSub_all(:,1))))];

semresample_shuff_gaussian  = [std(resample_shuff_all(:,2), 'omitnan') / sqrt(sum(~isnan(resample_shuff_all(:,2)))), ...
                          std(resample_shuff_all(:,1), 'omitnan') / sqrt(sum(~isnan(resample_shuff_all(:,1))))];

%% Plot
figure('Position', [100, 100, 800, 600]); 
tiledlayout(4,5);
subplotIdx = 1;
maxV = 0.11;
minV = 0;
plotJointInfo(meanJointA1, semJointA1, subplotIdx);
subplotIdx = subplotIdx + 1;
plotPIDterms(meanpluginA1, sempluginA1, subplotIdx, 'plugin', minV, maxV, false);
hold on;
add_pvalue(p_A1_plugin, 0.75, 2.25, 0.002, 0.1);
subplotIdx = subplotIdx + 1;
plotPIDterms(meanqeA1, semqeA1, subplotIdx, 'qe', minV,maxV, false);
hold on;
add_pvalue(p_A1_qe, 0.75, 2.25, 0.002, 0.1);
subplotIdx = subplotIdx + 1;
plotPIDterms(meanshuffsubA1, semshuffsubA1, subplotIdx, 'shuffSub',minV,maxV, false);
hold on;
add_pvalue(p_A1_shuffSub, 0.75, 2.25, 0.002, 0.1);
subplotIdx = subplotIdx + 1;
plotPIDterms(meanWeightedA1, semWeightedA1, subplotIdx, 'merged', minV, maxV, false);
add_pvalue(p_A1_weighted, 0.75, 2.25, 0.002, 0.1)



maxV = 0.025;
subplotIdx = subplotIdx + 1;
plotJointInfo(meanJointCA1, semJointCA1, subplotIdx);
subplotIdx = subplotIdx + 1;
plotPIDterms(meanpluginCA1, sempluginCA1, subplotIdx, 'plugin', minV,maxV, false);
hold on;
add_pvalue(p_CA1_plugin, 0.75, 2.25, 0.001, 0.022);
subplotIdx = subplotIdx + 1;
plotPIDterms(meanqeCA1, semqeCA1, subplotIdx, 'qe',minV,maxV, false);
hold on;
add_pvalue(p_CA1_qe, 0.75, 2.25, 0.001, 0.015);
subplotIdx = subplotIdx + 1;
plotPIDterms(meanshuffsubCA1, semshuffsubCA1, subplotIdx, 'shuffSub',minV,maxV, false);
hold on;
add_pvalue(p_CA1_shuffSub, 0.75, 2.25, 0.001, 0.015);
subplotIdx = subplotIdx + 1;
plotPIDterms(meanWeightedCA1, semWeightedCA1, subplotIdx, 'merged', minV, maxV, false);
add_pvalue(p_CA1_weighted, 0.75, 2.25, 0.001, 0.015);



maxV = 0.012;
minV =0;
subplotIdx = subplotIdx + 1;
plotJointInfo(meanJointA1_monkey, semJointA1_monkey, subplotIdx, true);
subplotIdx = subplotIdx + 1;
plotPIDterms(meanpluginA1_monkey, sempluginA1_monkey, subplotIdx, 'plugin', minV,maxV, false, true);
hold on;
add_pvalue(p_A1_monkey_plugin, 0.75, 2.25, 0.0007, 0.01);
subplotIdx = subplotIdx + 1;
plotPIDterms(meanqeA1_monkey, semqeA1_monkey, subplotIdx, 'qe', minV,maxV, false, true);
hold on;
add_pvalue(p_A1_monkey_qe, 0.75, 2.25, 0.001, 0.005);
subplotIdx = subplotIdx + 1;
plotPIDterms(meanshuffsubA1_monkey, semshuffsubA1_monkey, subplotIdx, 'shuffSub',minV,maxV, false, true);
hold on;
add_pvalue(p_A1_monkey_shuffSub, 0.75, 2.25,0.001, 0.005);
subplotIdx = subplotIdx + 1;
plotPIDterms(meanWeightedA1_monkey, semWeightedA1_monkey, subplotIdx, 'merged', minV, maxV, false, true);
add_pvalue(p_A1_monkey_weighted, 0.75, 2.25, 0.001, 0.005);

maxV = 9.5;
minV =0;
subplotIdx = subplotIdx + 1;
plotJointInfo(meanJoint_gaussian, semJoint_gaussian, subplotIdx, true);
subplotIdx = subplotIdx + 1;
plotPIDterms(meanplugin_gaussian, semplugin_gaussian, subplotIdx, 'plugin', minV,maxV, false, true);
hold on;
add_pvalue(p_gaussian_plugin, 0.75, 2.25, 0.1, 9);
subplotIdx = subplotIdx + 1;
plotPIDterms(meanresample_gaussian, semresample_gaussian, subplotIdx, 'resample', minV,maxV, false, true);
hold on;
add_pvalue(p_gaussian_resample, 0.75, 2.25, 0.1,   9);
subplotIdx = subplotIdx + 1;
plotPIDterms(meanshuffsub_gaussian, semshuffsub_gaussian, subplotIdx, 'shuffSub',minV,maxV, false, true);
hold on;
add_pvalue(p_gaussian_shuffle, 0.75, 2.25,0.1,  9);
subplotIdx = subplotIdx + 1;
plotPIDterms(meanresample_shuff_gaussian, semresample_shuff_gaussian, subplotIdx, 'merged',minV,maxV, false, true);
hold on;
add_pvalue(p_gaussian_resample_shuff, 0.75, 2.25, 0.1,  9);

outDir = fullfile(fileparts(pwd), 'Figures_mat');   % ein Verzeichnis über pwd
set(gcf, 'Renderer', 'painters');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end
exportgraphics(gcf, fullfile(outDir, 'Figure_6.svg'), 'ContentType','vector');



function ax = plotJointInfo(meanJoint, semJoint, subplotIdx,pltxlabel)
    if nargin<4
        pltxlabel = false;
    end
    if nargin<5
        pltylabel = false;
    end
    if nargin<6
        plttitle = false;
    end
    barLabels = {'Plugin', 'qe', 'shuffSub', 'weighted'};
    colorIJoint = [0.9290, 0.6940, 0.1250];
    ax = nexttile(subplotIdx);
    b = bar(meanJoint, 'FaceColor', colorIJoint, 'EdgeColor', 'k');
    hold on;
    errorbar(1:length(meanJoint), meanJoint, 3*semJoint, 'k', 'LineStyle', 'none', 'LineWidth', 1.5);
    xticks(1:length(barLabels));
    if true
        xticklabels(barLabels);
    else
        xticklabels([]);
    end
    if true
        ylabel('Information [bits]');
    else
        ylabel([]);
    end
    grid on;
    box on;
    hold off;
end
function ax = plotPIDterms(meanJoint, semJoint, subplotIdx, plotTitle, minV,maxLim, isPercent,pltxlabel,pltylabel,plttitle)
    if nargin<8
        pltxlabel = false;
    end
    if nargin<9
        pltylabel = false;
    end
    if nargin<10
        plttitle = false;
    end
    barLabels = {'Red', 'Syn'};
    colors = [0, 0.4470, 0.7410;  
              0.4660, 0.6740, 0.1880]; 
    ax =nexttile(subplotIdx);
    b = bar(meanJoint, 'FaceColor', 'flat', 'EdgeColor', 'k');
    for i = 1:length(meanJoint)
        b.CData(i, :) = colors(i, :);
    end
    hold on;
    errorbar(1:length(meanJoint), meanJoint, 3*semJoint, 'k', 'LineStyle', 'none', 'LineWidth', 1.5);
    xticks(1:length(barLabels));
    if true
        xticklabels(barLabels);
    else
        xticklabels([]);
    end
    if true
        if isPercent
            ylabel('% of Joint Info');
        else
            ylabel('Information [bits]');
        end
    else
        ylabel([]);
    end
    if true
        title(plotTitle, 'Interpreter', 'none');
    else
        title([]);
    end
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


% Apply weighting row-by-row to match the real-data table layout
%   qe_mat, shuff_mat : [nDatasets x nCols] (your columns: 1=Joint, 3=Red, 2=Syn, etc.)
%   max_info          : cap for Joint (same 2 bits you used: log2(4))
% Returns:
%   weighted_mat      : same size as inputs
%   w_used            : per-row weights for reference
function [weighted_mat, w_used] = apply_adaptive_weighting(qe_mat, shuff_mat, max_info, get_weight)
    if isempty(qe_mat)
        weighted_mat = qe_mat; w_used = [];
        return;
    end
    % Joint per dataset (column 1 in your real-data arrays)
    joint_qe    = qe_mat(:,1);
    joint_shuff = shuff_mat(:,1);
    joint_avg   = (joint_qe + joint_shuff)/2;

    % Convert to info fraction in [0,1]
    info_frac = joint_avg ./ max_info;

    % Row-wise weights
    w_used = arrayfun(get_weight, info_frac);

    % Blend all columns component-wise (same weight per row across columns)
    weighted_mat = w_used .* qe_mat + (1 - w_used) .* shuff_mat;
end
