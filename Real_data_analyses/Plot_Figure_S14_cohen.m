clc, clear, close all;


dataset = '256tp'; % TO BE CHANGED TO "256tp" PLOT 256 tp results
if strcmp(dataset, 'Full')
    filename = '../Results/Real_data_analysis/Gaussian_data_fullTimeseries.mat';
elseif strcmp(dataset, '256tp')
    filename = '../Results/Real_data_analysis/Gaussian_data.mat';
end

load(filename);
plugin_all      = [];
resample_all    = [];
resample_shuff_all   = [];
shuffSub_all    = [];

% Order: Syn, Red, Unq, Joint
plugin_all = squeeze(results(:, :, :, 1, :));
resample_all = squeeze(results(:, :, :, 2, :));
shuffSub_all = squeeze(results(:, :, :, 3, :));
Venkatesh_all = squeeze(results(:, :, :, 4, :));
resample_shuff_all = cat(6, resample_all, shuffSub_all);
resample_shuff_all = squeeze(mean(resample_shuff_all,6));

s = size(plugin_all);
plugin_all = reshape(plugin_all, s(1)*s(2), s(3), s(4));
resample_all = reshape(resample_all, s(1)*s(2), s(3), s(4));
shuffSub_all = reshape(shuffSub_all, s(1)*s(2), s(3), s(4));
Venkatesh_all = reshape(Venkatesh_all, s(1)*s(2), s(3), s(4));
resample_shuff_all = reshape(resample_shuff_all, s(1)*s(2), s(3), s(4));

plugin_all(:,:,3)         = plugin_all(:,:,3) ./ 2;
resample_all(:,:,3)       = resample_all(:,:,3) / 2;
shuffSub_all(:,:,3)       = shuffSub_all(:,:,3) / 2;
Venkatesh_all(:,:,3)      = Venkatesh_all(:,:,3) / 2;
resample_shuff_all(:,:,3) = resample_shuff_all(:,:,3) / 2;


n = size(plugin_all, 1);  % =1500

mean_plugin         = squeeze(mean(plugin_all, 1, 'omitnan'));
sem_plugin          = squeeze(2 * std(plugin_all, 0, 1, 'omitnan') ./ sqrt(n));

n_resample          = sum(~isnan(resample_all(:,1,1)));
mean_resample       = squeeze(mean(resample_all, 1, 'omitnan'));
sem_resample        = squeeze(2 * std(resample_all, 0, 1, 'omitnan') ./ sqrt(n));

n_shuffSub          = sum(~isnan(shuffSub_all(:,1,1)));
mean_shuffSub       = squeeze(mean(shuffSub_all, 1, 'omitnan'));
sem_shuffSub        = squeeze(2 * std(shuffSub_all, 0, 1, 'omitnan') ./ sqrt(n));

n_Venkatesh         = sum(~isnan(Venkatesh_all(:,1,1)));
mean_Venkatesh      = squeeze(mean(Venkatesh_all, 1, 'omitnan'));
sem_Venkatesh       = squeeze(2 * std(Venkatesh_all, 0, 1, 'omitnan') ./ sqrt(n));

n_resample_shuff    = sum(~isnan(resample_shuff_all(:,1,1)));
mean_resample_shuff = squeeze(mean(resample_shuff_all, 1, 'omitnan'));
sem_resample_shuff  = squeeze(2 * std(resample_shuff_all, 0, 1, 'omitnan') ./ sqrt(n));


colors = [...
    51, 51, 172;  
    75, 148, 56;  
    171, 67, 63;  
    187, 188, 89;  
    102, 102, 102]./255; 
labels = {'plugin', 'resample', 'shuffSub', 'Venkatesh', 'resample\_shuff'};
labels2 = {'Joint', 'Synergy', 'Redundancy', 'Unique'};

means = {mean_plugin, mean_resample, mean_shuffSub, mean_Venkatesh, mean_resample_shuff};
sems  = {sem_plugin,  sem_resample,  sem_shuffSub,  sem_Venkatesh,  sem_resample_shuff};

x = 1:20;

dIdx = [4, 1, 2, 3];

figure;
for d = 1:4
    dVal = dIdx(d);
    subplot(2,2,d); hold on;    
    for i = 1:numel(means)
        if d == 4
            m = means{i}(:, dVal)/2;
        else
            m = means{i}(:, dVal);
        end
        s = sems{i}(:, dVal);
       fill([x, fliplr(x)], [m-s; flipud(m+s)]', ...
     colors(i,:), 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        plot(x, m, 'Color', colors(i,:), 'LineWidth', 1.5);
    end
    title(labels2{d});
    xlabel('dimension');
    ylabel('Information [bits]');
    xlim([1 20]);
   xticks([1, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20]);
end
legend(labels)

dir_to_save = '../Figures/';
fileName = 'RealDataGaussFullTimeseries.svg';
filePath = fullfile(dir_to_save, fileName);
saveas(gca, filePath, 'svg');




%%
% Define method labels
method_labels = {'Plugin', 'Resample', 'Shuffle', 'Merged', 'Venkatesh'};
means = {mean_plugin, mean_resample, mean_shuffSub, mean_resample_shuff, mean_Venkatesh};
sems  = {sem_plugin,  sem_resample,  sem_shuffSub ,  sem_resample_shuff,   sem_Venkatesh};
dim_idx=20;
% Initialize arrays to store results
synergy_vals = zeros(5,2);  % rows = methods, cols = [mean, 2SEM]
redund_vals  = zeros(5,2);

% Compute mean and 2SEM for synergy (column 2) and redundancy (column 3)
for i = 1:5
    synergy_mean = mean(means{i}(dim_idx,1), 'omitnan');
    synergy_sem  = mean(sems{i}(dim_idx,1), 'omitnan');

    redund_mean = mean(means{i}(dim_idx,2), 'omitnan');
    redund_sem  = mean(sems{i}(dim_idx,2), 'omitnan');

    synergy_vals(i,:) = [synergy_mean, synergy_sem];
    redund_vals(i,:)  = [redund_mean,  redund_sem];
end

% Print Markdown table header
fprintf('|Measure|%s|\n', strjoin(method_labels, '|'));
fprintf('|-|%s|\n', repmat('-|', 1, numel(method_labels)));

% Print Synergy row
fprintf('|Synergy');
for i = 1:5
    fprintf('|%.3f ± %.3f', synergy_vals(i,1), synergy_vals(i,2));
end
fprintf('|\n');

% Print Redundancy row
fprintf('|Redundancy');
for i = 1:5
    fprintf('|%.3f ± %.3f', redund_vals(i,1), redund_vals(i,2));
end
fprintf('|\n');
