function PlotFigure2gauss()
% Three-panel Gaussian (plugin) bias figure:
% 1) bias vs trials  (d=20, alpha=3)
% 2) bias vs dims M  (trials=256, alpha=3)
% 3) bias vs alpha   (trials=256, d=20)
%
% Lines: Joint (comp=1), Syn (comp=6), Red (comp=5)
% Bias = mean(measured) - ground truth (GT)
%
% Output: PNG + SVG saved next to the script.

clc;

%% ---------------------- Locate and load data ----------------------
resultFolder = 'Results/Gauss';
file_main    = fullfile(resultFolder, 'Finalresults_across_M_and_ntrials.mat');
file_d80     = fullfile(resultFolder, 'Finalresults_across_M_and_ntrials_d80.mat');

if isfile(file_main)
    D = load(file_main);
elseif isfile(file_d80)
    D = load(file_d80);
else
    error(['Could not find Gaussian results file:\n  %s\nor\n  %s'], file_main, file_d80);
end

% Expect:
%   sampled_results: comps x nTrials x nM x nBias x nAlpha x nIters
%   GT_results     : comps x nTrials x nM x nBias x nAlpha x nIters OR comps x nTrials x nM x nBias x nAlpha
%   ntrials_vals   : [1 x nTrials]
%   M_vals         : [1 x nM]
%
Samp = D.sampled_results;
GT   = D.GT_results;
ntrials_vals = D.ntrials_vals;
M_vals       = D.M_vals;

% Indices/constants
comp_joint = 1; comp_syn = 6; comp_red = 5;
bias_plugin = 1;                % plugin correction index in your arrays
alpha_idx   = 3;                % fixed alpha for panels 1 & 2
M_fixed     = 20;               % fixed M for panel 1 & 3
trials_fixed = 128;             % fixed trials for panel 2 & 3

% Find index for d=20 and trials=256 (closest match if exact not present)
[~, M_idx]      = min(abs(M_vals - M_fixed));
[~, trials_idx] = min(abs(ntrials_vals - trials_fixed));

% Iteration count (last dim)
nIters = size(Samp, 6);

% Helper to get GT slice consistently (average across GT iterations if needed)
getGT = @(A) squeeze(mean(A, ndims(A), 'omitnan'));  % average over iter dim if present

%% ---------------------- Panel 1: bias vs trials (d=20, alpha=3) ----------------------
% Slice: comps x nTrials x 1 x 1 x 1 x nIters  (then squeeze to comps x nTrials x nIters)
S1  = squeeze(Samp(:, :, M_idx, bias_plugin, alpha_idx, :));  % comps x nTrials x nIters
GT1 = squeeze(GT(:,  :, M_idx, bias_plugin, alpha_idx));   % comps x nTrials x (maybe nIters)
%GT1 = getGT(GT1);                                             % comps x nTrials

[joint_m1, joint_s1] = mean_sem_bias(S1(comp_joint,:,:), GT1(comp_joint,:));
[syn_m1,   syn_s1  ] = mean_sem_bias(S1(comp_syn,  :,:), GT1(comp_syn,  :));
[red_m1,   red_s1  ] = mean_sem_bias(S1(comp_red,  :,:), GT1(comp_red,  :));

%% ---------------------- Panel 2: bias vs dimensions (trials=256, alpha=3) ----------------------
% Slice: comps x 1 x nM x 1 x 1 x nIters  -> comps x nM x nIters
S2  = squeeze(Samp(:, trials_idx, :, bias_plugin, alpha_idx, :));  % comps x nM x nIters
GT2 = squeeze(GT(:,  trials_idx, :, bias_plugin, alpha_idx));   % comps x nM x (maybe nIters)
%GT2 = getGT(GT2);                                                  % comps x nM

[joint_m2, joint_s2] = mean_sem_bias(S2(comp_joint,:,:), GT2(comp_joint,:));
[syn_m2,   syn_s2  ] = mean_sem_bias(S2(comp_syn,  :,:), GT2(comp_syn,  :));
[red_m2,   red_s2  ] = mean_sem_bias(S2(comp_red,  :,:), GT2(comp_red,  :));

%% ---------------------- Panel 3: bias vs alpha (trials=256, d=20) ----------------------
nAlpha = size(Samp, 5);
alphas = [1., 1.1, 1.2, 1.5, 2., 10000.]; % 1:nAlpha;  % if you have actual alpha values vector, replace this

% Slice: comps x 1 x 1 x 1 x nAlpha x nIters -> comps x nAlpha x nIters
S3  = squeeze(Samp(:, trials_idx, M_idx, bias_plugin, :, :));   % comps x nAlpha x nIters
GT3 = squeeze(GT(:,  trials_idx, M_idx, bias_plugin, :));    % comps x nAlpha x (maybe nIters)
%GT3 = getGT(GT3);                                               % comps x nAlpha

[joint_m3, joint_s3] = mean_sem_bias(S3(comp_joint,:,:), GT3(comp_joint,:));
[syn_m3,   syn_s3  ] = mean_sem_bias(S3(comp_syn,  :,:), GT3(comp_syn,  :));
[red_m3,   red_s3  ] = mean_sem_bias(S3(comp_red,  :,:), GT3(comp_red,  :));

%% ---------------------- Plotting ----------------------
set(0,'DefaultTextFontName','Arial'); 
set(0,'DefaultAxesFontName','Arial'); 
set(0,'DefaultLegendFontName','Arial');

fig = figure('Units','centimeters','Position',[1 1 18 6]);
tl = tiledlayout(fig, 1, 3, 'TileSpacing','compact', 'Padding','compact');

lw  = 1;
clr_joint = [0.9290 0.6940 0.12580];  % green-ish
clr_red   = [0       0.4470 0.7410]; % blue-ish
clr_syn   = [0.4660 0.6740 0.1880];
alphaFill = 0.18;  % yellow-ish
FontSize = 9;
% -------- Panel 1: trials (log-x) --------
nexttile; hold on;
plot_with_ci(ntrials_vals, joint_m1, joint_s1, clr_joint, lw, alphaFill);
plot_with_ci(ntrials_vals, syn_m1,   syn_s1,   clr_syn,   lw, alphaFill);
plot_with_ci(ntrials_vals, red_m1,   red_s1,   clr_red,   lw, alphaFill);
set(gca,'XScale','log','LineWidth',lw,'FontSize',FontSize);
xticks([100, 1000])
xlabel('Trials (N)'); ylabel('Bias [bits]');
title('Trials (d=20)');
yline(0,'k-','LineWidth',lw,'HandleVisibility','off');
legend({'Joint','Syn','Red'},'Location','best','Box','off');

% -------- Panel 2: dimensions M --------
nexttile; hold on;
plot_with_ci(M_vals, joint_m2, joint_s2, clr_joint, lw, alphaFill);
plot_with_ci(M_vals, syn_m2,   syn_s2,   clr_syn,   lw, alphaFill);
plot_with_ci(M_vals, red_m2,   red_s2,   clr_red,   lw, alphaFill);
set(gca,'LineWidth',lw,'FontSize', FontSize);
xlabel('Number of dimensions (d)'); %ylabel('Bias [bits]');
title(sprintf('Dimensions (N=%d)', ntrials_vals(trials_idx)));
yline(0,'k-','LineWidth',lw,'HandleVisibility','off');

% -------- Panel 3: alpha sweep --------
nexttile; hold on;
plot_with_ci(GT3(comp_joint,:), joint_m3, joint_s3, clr_joint, lw, alphaFill);
plot_with_ci(GT3(comp_joint,  :), syn_m3,   syn_s3,   clr_syn,   lw, alphaFill);
plot_with_ci(GT3(comp_joint,  :), red_m3,   red_s3,   clr_red,   lw, alphaFill);
xlim([0, max(GT3(comp_joint,:))])
set(gca,'LineWidth',lw,'FontSize', FontSize);
xlabel('Information level [bits]'); % ylabel('Bias [bits]');
title(sprintf('\\alpha (N=%d, d=%d)', ntrials_vals(trials_idx), M_vals(M_idx)));
yline(0,'k-','LineWidth',lw,'HandleVisibility','off');

% Save
set(gcf,'Renderer','painters');
%print(gcf,'-dpng','gauss_plugin_bias_panels.png','-r300');outDir = 'Figures_mat';
outDir = 'Figures_mat';
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

exportgraphics(gcf, fullfile(outDir, 'Figure_2_gauss.svg'), 'ContentType','vector');

end % main function


%% ======================= Helpers =======================
function [mBias, sSEM] = mean_sem_bias(comp_x_iter, GT_vec)
% comp_x_iter: [nX x nIters] (after squeeze)
% GT_vec     : [1 x nX] or [nX x 1]
X = squeeze(comp_x_iter);           % nX x nIters
if size(X,1) == 1 && size(X,2) > 1
    X = X.'; % make it nX x nIters if needed
end
GT_vec = squeeze(GT_vec(:)).';      % row 1 x nX
GT_rep = repmat(GT_vec, size(X,2), 1).';  % nX x nIters

biasMat = X - GT_rep;               % nX x nIters
mBias = mean(biasMat, 2, 'omitnan');         % nX x 1
nIters = size(biasMat, 2);
sSEM  = 3 * std(biasMat, 0, 2, 'omitnan') / sqrt(max(nIters,1));  % 3*SEM band
end

function plot_with_ci(x, m, s, color, lw, a)
x = x(:).';
m = m(:).';
s = s(:).';
xx = [x, fliplr(x)];
yy = [m+s, fliplr(m-s)];
fill(xx, yy, color, 'FaceAlpha', a, 'EdgeColor','none', 'HandleVisibility','off');
plot(x, m, 'LineWidth', lw, 'Color', color);
end
