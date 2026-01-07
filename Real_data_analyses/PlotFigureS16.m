function PlotFigureS16(dataset)
% PlotFigureS16  Reproduces the Gaussian real-data plot (Figure S16) and saves as SVG.
%
% Usage:
%   PlotFigureS16('256tp')
%   PlotFigureS16('Full')
%
% Input:
%   dataset : '256tp' or 'Full'
%
% Output:
%   Creates ../Figures_mat/Figure_S16.svg (one level above current pwd)

if nargin < 1 || isempty(dataset)
    dataset = '256tp';
end

% --- Select file based on dataset ---
switch lower(dataset)
    case 'full'
        filename = 'Results/Real_data_analysis/Gaussian_data_fullTimeseries.mat';
    case '256tp'
        filename = 'Results/Real_data_analysis/Gaussian_data.mat';
    otherwise
        error('Unknown dataset "%s". Use "256tp" or "Full".', dataset);
end

% --- Load ---
D = load(filename);
if ~isfield(D, 'results')
    error('Expected variable "results" not found in %s.', filename);
end
results = D.results;

% --- Extract methods (Order in results: method index 1..4) ---
% Comment from your script: Order within "results": Syn, Red, Unq, Joint
plugin_all    = squeeze(results(:, :, :, 1, :));
resample_all  = squeeze(results(:, :, :, 2, :));
shuffSub_all  = squeeze(results(:, :, :, 3, :));
Venkatesh_all = squeeze(results(:, :, :, 4, :));

% merged: resample + shuffSub average
resample_shuff_all = cat(6, resample_all, shuffSub_all);
resample_shuff_all = squeeze(mean(resample_shuff_all, 6));

% --- Reshape to [N x dims x components] ---
s = size(plugin_all); % [A x B x dims x comps] after squeeze
plugin_all        = reshape(plugin_all,        s(1)*s(2), s(3), s(4));
resample_all      = reshape(resample_all,      s(1)*s(2), s(3), s(4));
shuffSub_all      = reshape(shuffSub_all,      s(1)*s(2), s(3), s(4));
Venkatesh_all     = reshape(Venkatesh_all,     s(1)*s(2), s(3), s(4));
resample_shuff_all= reshape(resample_shuff_all,s(1)*s(2), s(3), s(4));

% --- Unique term handling (your original: divide component 3 by 2) ---
plugin_all(:,:,3)          = plugin_all(:,:,3) ./ 2;
resample_all(:,:,3)        = resample_all(:,:,3) ./ 2;
shuffSub_all(:,:,3)        = shuffSub_all(:,:,3) ./ 2;
Venkatesh_all(:,:,3)       = Venkatesh_all(:,:,3) ./ 2;
resample_shuff_all(:,:,3)  = resample_shuff_all(:,:,3) ./ 2;

% --- Means & SEMs (same logic as your script: 2*SEM) ---
n = size(plugin_all, 1);

mean_plugin         = squeeze(mean(plugin_all,        1, 'omitnan'));
sem_plugin          = squeeze(2 * std(plugin_all,     0, 1, 'omitnan') ./ sqrt(n));

mean_resample       = squeeze(mean(resample_all,      1, 'omitnan'));
sem_resample        = squeeze(2 * std(resample_all,   0, 1, 'omitnan') ./ sqrt(n));

mean_shuffSub       = squeeze(mean(shuffSub_all,      1, 'omitnan'));
sem_shuffSub        = squeeze(2 * std(shuffSub_all,   0, 1, 'omitnan') ./ sqrt(n));

mean_Venkatesh      = squeeze(mean(Venkatesh_all,     1, 'omitnan'));
sem_Venkatesh       = squeeze(2 * std(Venkatesh_all,  0, 1, 'omitnan') ./ sqrt(n));

mean_resample_shuff = squeeze(mean(resample_shuff_all,    1, 'omitnan'));
sem_resample_shuff  = squeeze(2 * std(resample_shuff_all, 0, 1, 'omitnan') ./ sqrt(n));

% --- Colors/labels (unchanged from your script) ---
colors = [...
    0.60,  0.60,  0.60;   % plugin
    0.929, 0.694, 0.125;    % resample
    0.850, 0.325, 0.098;    % shuffSub
    0.466, 0.674, 0.188;    % Venkatesh
    0.494, 0.184, 0.556]; % resample_shuff

labels  = {'plugin', 'resample', 'shuffle', 'Venkatesh', 'merged'};
labels2 = {'Joint', 'Synergy', 'Redundancy', 'Unique'};

means = {mean_plugin, mean_resample, mean_shuffSub, mean_Venkatesh, mean_resample_shuff};
sems  = {sem_plugin,  sem_resample,  sem_shuffSub,  sem_Venkatesh,  sem_resample_shuff};

% x-axis / mapping (unchanged)
x = 1:20;
dIdx = [4, 1, 2, 3]; % (Joint, Syn, Red, Unq) mapped from [Syn Red Unq Joint]

% --- Plot ---
figure('Color','w');
for d = 1:4
    dVal = dIdx(d);
    subplot(2,2,d); hold on;

    for i = 1:numel(means)
        if d == 4
            m = means{i}(:, dVal) / 2; % keep your original special-case
        else
            m = means{i}(:, dVal);
        end
        s_ = sems{i}(:, dVal);

        fill([x, fliplr(x)], [m - s_; flipud(m + s_)]', ...
            colors(i,:), 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        plot(x, m, 'Color', colors(i,:), 'LineWidth', 1.5);
    end

    title(labels2{d});
    xlabel('dimension');
    ylabel('Information [bits]');
    xlim([1 20]);
    xticks([1, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20]);
    box off;
end
legend(labels, 'Interpreter','none');

% --- Save (parent of pwd / Figures_mat) ---
outDir = fullfile(pwd, 'Figures_mat');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end


switch lower(dataset)
    case 'full'
        outFile = fullfile(outDir, 'Figure_S16_full.svg');
    case '256tp'
        outFile = fullfile(outDir, 'Figure_S16_256tp.svg');
end
exportgraphics(gcf, outFile, 'ContentType','vector');

end
