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
%% Prepare A1 data
% filename = ['compare_bins/Real_data_analysis/A1_3bins', '.mat'];
filename = ['Data/DATA_A1/Real_data_analysis/A1_3bins', '.mat'];

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


%% Prepare CA1 data
filename = ['compare_bins/Real_data_analysis/A1_2bins','.mat'];
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


%% Prepare A1_monkey data
filename = ['compare_bins/Real_data_analysis/KayserA1_PID_' chosenAtom '.mat'];
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


%% Prepare Gaussian data
filename = ['compare_bins/Real_data_analysis/Gaussian_data.mat'];
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
maxV = 0.2;
minV = 0;
plotJointInfo(meanJointA1, semJointA1, subplotIdx)
subplotIdx = subplotIdx + 1;
plotPIDterms(meanpluginA1, sempluginA1, subplotIdx, 'plugin', minV, maxV, false)
hold on;
add_pvalue(p_A1_plugin, 0.75, 2.25, 0.002, 0.1)
subplotIdx = subplotIdx + 1;
plotPIDterms(meanqeA1, semqeA1, subplotIdx, 'qe', minV,maxV, false)
hold on;
add_pvalue(p_A1_qe, 0.75, 2.25, 0.002, 0.1)
subplotIdx = subplotIdx + 1;
plotPIDterms(meanshuffsubA1, semshuffsubA1, subplotIdx, 'shuffSub',minV,maxV, false)
hold on;
add_pvalue(p_A1_shuffSub, 0.75, 2.25, 0.002, 0.1)
subplotIdx = subplotIdx + 1;
plotPIDterms(meanqeShuffsubA1, semqeShuffsubA1, subplotIdx, 'qe shuffSub',minV,maxV, false)
hold on;
add_pvalue(p_A1_qeshuffSub, 0.75, 2.25, 0.002, 0.1)

maxV = 0.2;
subplotIdx = subplotIdx + 1;
plotJointInfo(meanJointCA1, semJointCA1, subplotIdx)
subplotIdx = subplotIdx + 1;
plotPIDterms(meanpluginCA1, sempluginCA1, subplotIdx, 'plugin', minV,maxV, false)
hold on;
add_pvalue(p_CA1_plugin, 0.75, 2.25, 0.001, 0.022)
subplotIdx = subplotIdx + 1;
plotPIDterms(meanqeCA1, semqeCA1, subplotIdx, 'qe',minV,maxV, false)
hold on;
add_pvalue(p_CA1_qe, 0.75, 2.25, 0.001, 0.015)
subplotIdx = subplotIdx + 1;
plotPIDterms(meanshuffsubCA1, semshuffsubCA1, subplotIdx, 'shuffSub',minV,maxV, false)
hold on;
add_pvalue(p_CA1_shuffSub, 0.75, 2.25, 0.001, 0.015)
subplotIdx = subplotIdx + 1;
plotPIDterms(meanqeShuffsubCA1, semqeShuffsubCA1, subplotIdx, 'qe shuffSub',minV,maxV, false)
hold on;
add_pvalue(p_CA1_qeshuffSub, 0.75, 2.25, 0.001, 0.015)

maxV = 0.012;
minV =0;
subplotIdx = subplotIdx + 1;
plotJointInfo(meanJointA1_monkey, semJointA1_monkey, subplotIdx)
subplotIdx = subplotIdx + 1;
plotPIDterms(meanpluginA1_monkey, sempluginA1_monkey, subplotIdx, 'plugin', minV,maxV, false)
hold on;
add_pvalue(p_A1_monkey_plugin, 0.75, 2.25, 0.0007, 0.01)
subplotIdx = subplotIdx + 1;
plotPIDterms(meanqeA1_monkey, semqeA1_monkey, subplotIdx, 'qe', minV,maxV, false)
hold on;
add_pvalue(p_A1_monkey_qe, 0.75, 2.25, 0.001, 0.005)
subplotIdx = subplotIdx + 1;
plotPIDterms(meanshuffsubA1_monkey, semshuffsubA1_monkey, subplotIdx, 'shuffSub',minV,maxV, false)
hold on;
add_pvalue(p_A1_monkey_shuffSub, 0.75, 2.25,0.001, 0.005)
subplotIdx = subplotIdx + 1;
plotPIDterms(meanqeShuffsubA1_monkey, semqeShuffsubA1_monkey, subplotIdx, 'qe shuffSub',minV,maxV, false)
hold on;
add_pvalue(p_A1_monkey_qeshuffSub, 0.75, 2.25, 0.001, 0.005)

maxV = 9.5;
minV =0;
subplotIdx = subplotIdx + 1;
plotJointInfo(meanJoint_gaussian, semJoint_gaussian, subplotIdx)
subplotIdx = subplotIdx + 1;
plotPIDterms(meanplugin_gaussian, semplugin_gaussian, subplotIdx, 'plugin', minV,maxV, false)
hold on;
add_pvalue(p_gaussian_plugin, 0.75, 2.25, 0.1, 9)
subplotIdx = subplotIdx + 1;
plotPIDterms(meanresample_gaussian, semresample_gaussian, subplotIdx, 'resample', minV,maxV, false)
hold on;
add_pvalue(p_gaussian_resample, 0.75, 2.25, 0.1,   9)
subplotIdx = subplotIdx + 1;
plotPIDterms(meanshuffsub_gaussian, semshuffsub_gaussian, subplotIdx, 'shuffSub',minV,maxV, false)
hold on;
add_pvalue(p_gaussian_shuffle, 0.75, 2.25,0.1,  9)
subplotIdx = subplotIdx + 1;
plotPIDterms(meanresample_shuff_gaussian, semresample_shuff_gaussian, subplotIdx, 'resample_shuff',minV,maxV, false)
hold on;
add_pvalue(p_gaussian_resample_shuff, 0.75, 2.25, 0.1,  9)



filename = ['Figure_5_' chosenAtom '.svg'];
saveas(gcf, filename);
filename = ['Figure_5_' chosenAtom '.png'];
print(filename, '-dpng', '-r300');

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