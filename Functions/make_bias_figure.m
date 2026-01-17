function make_bias_figure()
% make_bias_figure
% Builds a 1x3 figure with bias (measured - GT) curves for Joint, Synergy, Redundancy
% using the PLUGIN correction method in the DISCRETE case.
%
% PANEL 1 (left):    bias vs number of trials,  bins = 4, alpha = 3
% PANEL 2 (middle):  bias vs number of bins,    trials = 512, alpha = 3
% PANEL 3 (right):   bias vs alpha values,      trials = 512, bins = 4
%                    GT for this panel = SHUFFLE at the HIGHEST trial count.
%
% This script assumes the following files (as in your previous codebase):
%   Results/Broja/Simuldata_bit_of_all_high_I_BROJA.mat
%   Results/Broja/GroundTruth_bit_of_all_high.mat
%   Results/Bin_sweep/Simuldata_bit_of_all_low_I_BROJA.mat
%   Results/Bin_sweep/GroundTruth_bit_of_all_low.mat
%   Results/Alpha_sweep/Simuldata_bit_of_all_alpha_sweep_I_BROJA.mat
%
% If your filenames differ, adjust the `cfg` section below.

% ----------------------- CONFIG -----------------------
cfg.simul_case          = 'bit_of_all';
cfg.redundancy_measure  = 'I_BROJA';

% --- Panel 1: trials sweep (discrete, "Broja" folder in your old code) ---
cfg.p1.folder           = 'Results/Broja';
cfg.p1.info_amount      = 'high';
cfg.p1.trial_list       = [16 32 64 128 256 512 1024 2048];
cfg.p1.n_bins_fixed     = 4;     % as requested
cfg.p1.alpha_fixed      = 3;     % as requested

% --- Panel 2: bins sweep (discrete, Bin_sweep) ---
cfg.p2.folder           = 'Results/Bin_sweep';
cfg.p2.info_amount      = 'low'; % file naming in your original code
cfg.p2.bin_list         = 2:9;
cfg.p2.trials_fixed     = 512;   % as requested
cfg.p2.alpha_fixed      = 3;     % as requested

% --- Panel 3: alpha sweep (discrete, Alpha_sweep) ---
cfg.p3.folder           = 'Results/Alpha_sweep';
cfg.p3.info_amount      = 'alpha_sweep'; % naming used previously
cfg.p3.alpha_list       = 1:15;
cfg.p3.trials_fixed     = 512;   % as requested (for plugin curve)
cfg.p3.n_bins_fixed     = 4;     % as requested
% IMPORTANT: GT for panel 3 will be PID_v_shuff at the MAX trial count present.

% Line/plot settings
FontSize = 9;
lw       = 1.2;
cols     = struct('Joint',[0.9290,0.6940,0.1250], ...
                  'Syn',  [0.4660,0.6740,0.1880], ...
                  'Red',  [0.0000,0.4470,0.7410]);

% ----------------------- BUILD FIGURE -----------------------
figure('Units','centimeters','Position',[1 1 15 5]);
tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

% ==== PANEL 1: bias vs trials (bins=4, alpha=3) ====
nexttile; hold on;
[x_trials, bias_trials] = load_bias_vs_trials(cfg);
plot(x_trials, bias_trials.Joint, 'Color', cols.Joint, 'LineWidth', lw);
plot(x_trials, bias_trials.Syn,   'Color', cols.Syn,   'LineWidth', lw);
plot(x_trials, bias_trials.Red,   'Color', cols.Red,   'LineWidth', lw);
yline(0,'k-','LineWidth',0.75,'HandleVisibility','off');
set(gca,'XScale','log','FontSize',FontSize,'LineWidth',0.75);
xlabel('trials'); ylabel('bias [bits]');
title('Trials (bins=4, \alpha=3)');
legend({'Joint','Synergy','Redundancy'},'Box','off','Location','best');

% ==== PANEL 2: bias vs bins (trials=512, alpha=3) ====
nexttile; hold on;
[x_bins, bias_bins] = load_bias_vs_bins(cfg);
plot(x_bins, bias_bins.Joint, 'Color', cols.Joint, 'LineWidth', lw);
plot(x_bins, bias_bins.Syn,   'Color', cols.Syn,   'LineWidth', lw);
plot(x_bins, bias_bins.Red,   'Color', cols.Red,   'LineWidth', lw);
yline(0,'k-','LineWidth',0.75,'HandleVisibility','off');
set(gca,'FontSize',FontSize,'LineWidth',0.75);
xlabel('number of bins'); yticks([]); % left panel already has y-label
title('Bins (trials=512, \alpha=3)');

% ==== PANEL 3: bias vs alpha (trials=512, bins=4; GT = shuff @ max trials) ====
nexttile; hold on;
[x_alpha, bias_alpha] = load_bias_vs_alpha_gtIsShuffMaxTrial(cfg);
plot(x_alpha, bias_alpha.Joint, 'Color', cols.Joint, 'LineWidth', lw);
plot(x_alpha, bias_alpha.Syn,   'Color', cols.Syn,   'LineWidth', lw);
plot(x_alpha, bias_alpha.Red,   'Color', cols.Red,   'LineWidth', lw);
yline(0,'k-','LineWidth',0.75,'HandleVisibility','off');
set(gca,'FontSize',FontSize,'LineWidth',0.75);
xlabel('\alpha'); yticks([]);
title('\alpha (trials=512, bins=4) — GT: shuff @ max trials');

set(gcf,'Renderer','painters');
% Uncomment if you want files:
% saveas(gcf, 'discrete_bias_panels.png');
% print(gcf, '-dsvg', 'discrete_bias_panels');
end % make_bias_figure


% ======================= HELPERS =======================

function [x_trials, bias] = load_bias_vs_trials(cfg)
% Loads PLUGIN estimates across trials (Broja) and subtracts scalar GT from GT file.

data_f = fullfile(cfg.p1.folder, ...
    sprintf('Simuldata_%s_%s_%s.mat', cfg.simul_case, cfg.p1.info_amount, cfg.redundancy_measure));
gt_f   = fullfile(cfg.p1.folder, ...
    sprintf('GroundTruth_%s_%s.mat',  cfg.simul_case, cfg.p1.info_amount));
needvar = 'PID_v_plugin';

assert_exist(data_f); assert_exist(gt_f);

S  = load(data_f, needvar);
if ~isfield(S, needvar), error('Variable "%s" not found in %s', needvar, data_f); end
PID = S.(needvar);   % nTrials x comps x iters

G  = load(gt_f, 'GroundTruth_value');
if ~isfield(G, 'GroundTruth_value'), error('GroundTruth_value not found in %s', gt_f); end

GTv = G.GroundTruth_value;
if iscell(GTv)
    gt_cell = pick_cell_for_discrete(GTv, cfg.p1.n_bins_fixed, cfg.p1.alpha_fixed);
    GT = gt_cell(:);
else
    GT = GTv(:);
end
% GT layout: [Red, U1, U2, Syn, IInd]
GT_joint = sum(GT(1:4));  GT_syn = GT(4);  GT_red = GT(1);

m_joint = squeeze(mean(PID(:,1,:), 3, 'omitnan'));
m_syn   = squeeze(mean(PID(:,2,:), 3, 'omitnan'));
m_red   = squeeze(mean(PID(:,3,:), 3, 'omitnan'));

bias.Joint = m_joint - GT_joint;
bias.Syn   = m_syn   - GT_syn;
bias.Red   = m_red   - GT_red;

x_trials = cfg.p1.trial_list(:);
end


function [x_bins, bias] = load_bias_vs_bins(cfg)
% Loads PLUGIN estimates across #bins at fixed trials (=512) (Bin_sweep)
% and subtracts per-bin GT from GT file.

data_f = fullfile(cfg.p2.folder, ...
    sprintf('Simuldata_%s_%s_%s.mat', cfg.simul_case, cfg.p2.info_amount, cfg.redundancy_measure));
gt_f   = fullfile(cfg.p2.folder, ...
    sprintf('GroundTruth_%s_%s.mat',  cfg.simul_case, cfg.p2.info_amount));
needvar = 'PID_v_plugin';

assert_exist(data_f); assert_exist(gt_f);

S = load(data_f, needvar);
if ~isfield(S, needvar), error('Variable "%s" not found in %s', needvar, data_f); end
PID_all = S.(needvar);  % nTrials x comps x iters x nBins

trial_axis = infer_trial_axis(PID_all, cfg.p2.trials_fixed);
trialIdx   = trial_axis.best_idx;

PID_slice = squeeze(PID_all(trialIdx, :, :, :)); % comps x iters x bins
PID_slice = permute(PID_slice, [3 1 2]);         % bins x comps x iters

G  = load(gt_f, 'GroundTruth_value');
if ~isfield(G, 'GroundTruth_value'), error('GroundTruth_value not found in %s', gt_f); end

[GT_bins, bins_vec] = gt_bins_from_cell(G.GroundTruth_value, cfg.p2.bin_list);
GT_joint = sum(GT_bins(:,1:4), 2); GT_syn = GT_bins(:,4); GT_red = GT_bins(:,1);

m_joint = squeeze(mean(PID_slice(:,1,:), 3, 'omitnan'));
m_syn   = squeeze(mean(PID_slice(:,2,:), 3, 'omitnan'));
m_red   = squeeze(mean(PID_slice(:,3,:), 3, 'omitnan'));

bias.Joint = m_joint - GT_joint;
bias.Syn   = m_syn   - GT_syn;
bias.Red   = m_red   - GT_red;

x_bins = bins_vec(:);
end


function [x_alpha, bias] = load_bias_vs_alpha_gtIsShuffMaxTrial(cfg)
% Panel 3 requirement: Use PLUGIN at fixed trials (e.g., 512) for measured,
% but use SHUFFLE at the HIGHEST available trial count as the "ground truth".
%
% Data assumed in: Results/Alpha_sweep/Simuldata_<case>_alpha_sweep_<measure>.mat
% Variables: PID_v_plugin, PID_v_shuff with dims (nTrials x comps x iters x nAlpha)

data_f   = fullfile(cfg.p3.folder, ...
            sprintf('Simuldata_%s_%s_%s200.mat', cfg.simul_case, cfg.p3.info_amount, cfg.redundancy_measure));
var_plug = 'PID_v_plugin';
var_shuf = 'PID_v_shuff';

assert_exist(data_f);

L = load(data_f);
if ~isfield(L, var_plug), error('Variable "%s" not found in %s', var_plug, data_f); end
if ~isfield(L, var_shuf), error('Variable "%s" not found in %s', var_shuf, data_f); end

PID_plug  = L.(var_plug); % nTrials x comps x iters x nAlpha
PID_shuff = L.(var_shuf); % nTrials x comps x iters x nAlpha

% ----- choose measured (plugin) at desired trials -----
trial_axis = infer_trial_axis(PID_plug, cfg.p3.trials_fixed);
trialIdx_meas = trial_axis.best_idx;

% Slice plugin at trialIdx_meas -> alpha x comps x iters
PL = squeeze(PID_plug(trialIdx_meas, :, :, :));   % comps x iters x alpha
PL = permute(PL, [3 1 2]);                        % alpha x comps x iters

% ----- choose GT (shuffle) at MAX trials available -----
nTrialsAvail = size(PID_shuff, 1);
trialIdx_GT  = nTrialsAvail;                      % highest available trial index
SH = squeeze(PID_shuff(trialIdx_GT, :, :, :));    % comps x iters x alpha
SH = permute(SH, [3 1 2]);                        % alpha x comps x iters

% Average over iterations
m_joint_PL = squeeze(mean(PL(:,1,:), 3, 'omitnan'));  % alpha x 1
m_syn_PL   = squeeze(mean(PL(:,2,:), 3, 'omitnan'));
m_red_PL   = squeeze(mean(PL(:,3,:), 3, 'omitnan'));

gt_joint_SH = squeeze(mean(SH(:,1,:), 3, 'omitnan')); % alpha x 1
gt_syn_SH   = squeeze(mean(SH(:,2,:), 3, 'omitnan'));
gt_red_SH   = squeeze(mean(SH(:,3,:), 3, 'omitnan'));

% Bias (plugin − shuffle@maxtrials) per alpha
bias.Joint = m_joint_PL - gt_joint_SH;
bias.Syn   = m_syn_PL   - gt_syn_SH;
bias.Red   = m_red_PL   - gt_red_SH;

% x-axis = alpha list (assumed to match 1:nAlpha or provided cfg list)
nAlpha = size(PL,1);
if ~isempty(cfg.p3.alpha_list) && numel(cfg.p3.alpha_list) == nAlpha
    x_alpha = cfg.p3.alpha_list(:);
else
    x_alpha = (1:nAlpha).';
end
end


% ------------------- GT utilities (used by panels 1 & 2) -------------------

function cell_entry = pick_cell_for_discrete(GTv, n_bins, alpha_val)
% Try to find a GT cell that matches bins and/or alpha if encoded; fallback to last numeric.
cell_entry = [];
if iscell(GTv)
    for k = 1:numel(GTv)
        v = GTv{k};
        if isstruct(v)
            ok_bins  = ~isfield(v,'bins')  || v.bins  == n_bins;
            ok_alpha = ~isfield(v,'alpha') || v.alpha == alpha_val;
            if ok_bins && ok_alpha && isfield(v,'value')
                cell_entry = v.value; break;
            end
        elseif isnumeric(v)
            if numel(v) >= 4, cell_entry = v; break; end
        end
    end
    if isempty(cell_entry)
        for k = numel(GTv):-1:1
            v = GTv{k};
            if isnumeric(v) && numel(v) >= 4
                cell_entry = v; break;
            end
        end
    end
else
    cell_entry = GTv;
end
if isempty(cell_entry), error('Could not parse GroundTruth_value cell array.'); end
end

function [GT_bins, bins_vec] = gt_bins_from_cell(GTv, bins_wanted)
% Convert GT cell → [bins x comps]
if ~iscell(GTv), error('Expected GroundTruth_value to be a cell array for bin sweep.'); end
K = numel(GTv);
tmp = nan(K, 5);
for k = 1:K
    v = GTv{k};
    if isstruct(v) && isfield(v,'value'), v = v.value; end
    if isnumeric(v)
        v = v(:).'; tmp(k,1:min(5,numel(v))) = v(1:min(5,numel(v)));
    end
end
bins_vec = (2:K+1).'; % heuristic default
if nargin>=2 && ~isempty(bins_wanted)
    map = max(1, min(K, bins_wanted - 1));
    GT_bins = tmp(map, :); bins_vec = bins_wanted(:);
else
    GT_bins = tmp;
end
end


% ------------------- misc utilities -------------------

function out = infer_trial_axis(PID, desired_trials)
% Infer trial counts; return nearest index to desired_trials.
trial_candidates = [16 32 64 128 256 512 1024 2048 4096];
nT = size(PID, 1);
if nT <= numel(trial_candidates)
    tvals = trial_candidates(1:nT);
else
    tvals = 2.^(4:(4+nT-1));
end
[~, idx] = min(abs(tvals - desired_trials));
out.values   = tvals;
out.best_idx = idx;
end

function assert_exist(fname)
if ~exist(fname, 'file')
    error('File not found: %s', fname);
end
end
