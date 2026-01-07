function pid_plot_grid2(cfg)
%PID_PLOT_GRID  Generic grid plotter for PID simulations (discrete/Gauss/SVM).
%
% Required cfg fields (common):
%   cfg.ResultFolder        e.g. 'Broja' | 'IMMI' | 'Gauss' | 'SVM_results'
%   cfg.mode                'discrete' | 'gauss' | 'svm'      (auto if missing)
%   cfg.simul_cases         cellstr, rows (e.g., {'bit_of_all','uncorr_unique',...})
%   cfg.info_types          cellstr columns (subset of {'Joint','Syn','Red','Unq'})
%   cfg.bias_correction     cellstr (order = legend order)
%   cfg.colors              [N x 3] RGB aligned with bias_correction
%   cfg.row_names           cellstr pretty names per simul case (same length)
%
% Data/plot settings:
%   cfg.redundancy_measure  e.g. 'I_BROJA' or 'I_MMI' (discrete only)
%   cfg.info_amount         'low'|'high' (discrete, SVM) or numeric index (gauss)
%   cfg.trial_categories    numeric vector (e.g., [16 32 ... 2048])
%   cfg.trials_multiplier   default 1 (use 4 for SVM to match your originals)
%   cfg.row_ylim            [nRows x 2] [ymin ymax] per row (optional)
%   cfg.thumb_rule_x        scalar (optional; draws xline if provided)
%   cfg.save_path           'Figures/Figure_X.svg' (optional)
%
% Notes:
% - Discrete/SVM expect files: Results/<Folder>/Simuldata_<case>_<infoAmt>_<redund>.mat
%   containing variables PID_v_<bias>.
% - Discrete expects GT in Results/<Folder>/GroundTruth_<case>_<infoAmt>.mat with
%   GroundTruth_value (length 4: [Red, U1, U2, Syn]).
% - Gauss expects a row-level file:
%     Results/Gauss/Finalresults_across_M_and_ntrials[_<suffix>].mat
%   with sampled_results (comp x trials x M x bias x info x iter),
%   GT_results (comp x trials x M x bias x info), M_vals, ntrials_vals.
% -------------------------------------------------------------------------

% ---------- Defaults & cosmetics ----------
if ~isfield(cfg,'mode') || isempty(cfg.mode)
    if strcmpi(cfg.ResultFolder,'Gauss'), cfg.mode = 'gauss';
    elseif strcmpi(cfg.ResultFolder,'SVM_results'), cfg.mode = 'svm';
    else, cfg.mode = 'discrete';
    end
end
if ~isfield(cfg,'trials_multiplier') || isempty(cfg.trials_multiplier)
    cfg.trials_multiplier = strcmpi(cfg.mode,'svm') * 4 + ~strcmpi(cfg.mode,'svm') * 1;
end
if ~isfield(cfg,'trial_categories'), cfg.trial_categories = [16 32 64 128 256 512 1024 2048]; end
if ~isfield(cfg,'info_types') || isempty(cfg.info_types), cfg.info_types = {'Joint','Syn','Red'}; end
if ~isfield(cfg,'row_names') || isempty(cfg.row_names),  cfg.row_names  = cfg.simul_cases;        end

% Fonts
set(0,'DefaultTextFontName','Arial');
set(0,'DefaultAxesFontName','Arial');
set(0,'DefaultLegendFontName','Arial');
set(0,'DefaultTextFontSize',8);
set(0,'DefaultAxesFontSize',8);
set(0,'DefaultLegendFontSize',8);

% Figure / layout
nrows = numel(cfg.simul_cases);
ncols = numel(cfg.info_types);
figure_handle = figure('Units','centimeters','Position',[1,1,14,2.5*(nrows)]);
tl = tiledlayout(figure_handle, nrows, ncols, 'TileSpacing','compact','Padding','compact');

% ---------- Main loops ----------
for r = 1:nrows
    simul_case = cfg.simul_cases{r};

    switch lower(cfg.mode)
        case 'gauss'
            RowData = load_row_gauss(cfg, simul_case);
        otherwise % 'discrete' or 'svm'
            RowData = load_row_discrete(cfg, simul_case);
    end

    switch lower(cfg.mode)
        case 'gauss'
            x = RowData.trials(:);
        otherwise % 'discrete' or 'svm'
            x = 4*RowData.trials(:);
    end
    % row-specific X vector
    
    if isempty(x), continue; end

    for c = 1:ncols
        info_type = cfg.info_types{c};
        ax = nexttile(tl, (r-1)*ncols + c); hold(ax,'on');

        % Shaded SEM + mean for each bias
        Lh = gobjects(1, numel(cfg.bias_correction));
        for b = 1:numel(cfg.bias_correction)
            bc = cfg.bias_correction{b};
            bc_key = bias_key(bc);  % sanitized field name

            if ~isfield(RowData.means, bc_key), continue; end
            if ~isfield(RowData.means.(bc_key), info_type), continue; end

            mu  = RowData.means.(bc_key).(info_type)(:);
            sem = RowData.sems.(bc_key).(info_type)(:);

            % ---- align lengths to avoid fill errors ----
            n = min([numel(x), numel(mu), numel(sem)]);
            if n == 0, continue; end
            xx  = x(1:n);
            muu = mu(1:n);
            see = sem(1:n);

            % shaded band (~95% CI as 2*SEM)
            Xfill = [xx; flipud(xx)];
            Yfill = [muu + 2*see; flipud(muu - 2*see)];
            fill(Xfill, Yfill, cfg.colors(b,:), ...
                 'FaceAlpha',0.15, 'EdgeColor','none', 'HandleVisibility','off');

            % mean line
            Lh(b) = plot(xx, muu, 'LineWidth',1.5, 'Color', cfg.colors(b,:), ...
                         'DisplayName', pretty_label(bc));
        end

        % Ground truth: line or vector
        if isfield(RowData,'GT') && isfield(RowData.GT, info_type)
            gt = RowData.GT.(info_type);
            if isscalar(gt)
                plot([x(1) x(end)], [gt gt], 'k--', 'LineWidth',1.2);
            else
                % align GT length if needed
                ngt = min(numel(x), numel(gt));
                plot(x(1:ngt), gt(1:ngt), 'k--', 'LineWidth',1.2);
            end
        end

        % Optional vertical thumb rule
        if isfield(cfg,'thumb_rule_x') && ~isempty(cfg.thumb_rule_x)
            xline(cfg.thumb_rule_x, ':', 'Color',[0 0 0], 'LineWidth',1.0, 'HandleVisibility','off');
        end

        % Axes
        set(ax, 'XScale','log', 'LineWidth',1.2);
        xlim([x(1) x(end)]); xticks([100, 1000]);
        if r == nrows
            xlabel('Trials');
            % xticks_like(ax, x);
        else
            set(ax,'XTickLabel',[]);
        end
        if c == 1
            ylabel(sprintf('%s\nInfo [bits]', cfg.row_names{r}));
        else
            set(ax,'YTickLabel',[]);
        end
        yline(0,'k-','LineWidth',1.0,'HandleVisibility','off');

        % Row-specific y-lims if provided
        if isfield(cfg,'row_ylim') && size(cfg.row_ylim,1) >= r
            ylim(cfg.row_ylim(r,:));
        end

        % Titles on top row
        if r == 1, title(info_type, 'FontWeight','bold'); end

        % Legend once in top-right, using only valid handles
        if r == 1 && c == ncols
            validMask = isgraphics(Lh);
            H = Lh(validMask);
            if ~isempty(H)
                labels = get(H,'DisplayName');
                if ischar(labels), labels = {labels}; end
                lgd = legend(ax, H, labels, 'Location','northeast', 'Box','off');
                lgd.ItemTokenSize = [12,6];
            end
        end

        hold(ax,'off');
    end
end

% Save if requested
if isfield(cfg,'save_path') && ~isempty(cfg.save_path)
    outdir = fileparts(cfg.save_path);
    if ~exist(outdir,'dir'), mkdir(outdir); end
    [~,~,ext] = fileparts(cfg.save_path);
    if isempty(ext) || strcmpi(ext,'.svg')
        saveas(gcf, cfg.save_path, 'svg');
    else
        exportgraphics(tl, cfg.save_path, 'ContentType','vector');
    end
end

end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END MAIN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% ---------------------- Helpers: Data loading ---------------------------

function Row = load_row_discrete(cfg, simul_case)
% Loads per-bias PID_v_<bias> cubes and optional ground truth for a row.
% PID_v dims expected: [nTrials x nComponents x nIters]
Row.means = struct(); Row.sems = struct();

% X vector for this row (discrete/SVM use cfg)
Row.trials = (cfg.trial_categories(:)' * cfg.trials_multiplier).';

% Ground truth (if present)
gt_path = fullfile('Results', cfg.ResultFolder, ...
    sprintf('GroundTruth_%s_%s.mat', simul_case, cfg.info_amount));
if exist(gt_path,'file')
    S = load(gt_path, 'GroundTruth_value');
    if isfield(S,'GroundTruth_value')
        GT = S.GroundTruth_value(:);
        Row.GT = struct('Joint', sum(GT(1:4)), ...
                        'Syn',   GT(4), ...
                        'Red',   GT(1), ...
                        'Unq',   mean(GT(2:3)));
    end
end

% Data per bias
data_file = fullfile('Results', cfg.ResultFolder, ...
    sprintf('Simuldata_%s_%s_%s.mat', simul_case, cfg.info_amount, cfg.redundancy_measure));

for b = 1:numel(cfg.bias_correction)
    bc = cfg.bias_correction{b};
    var_name = sprintf('PID_v_%s', bc);
    if ~isfile(data_file), continue; end
    D = load(data_file, var_name);
    if ~isfield(D, var_name), continue; end
    PID_v = D.(var_name); % [T x comp x K]

    k = bias_key(bc);  % sanitize field
    Row.means.(k) = struct();
    Row.sems.(k)  = struct();

    Row.means.(k).Joint = mean_over_iter(PID_v(:,1,:));
    Row.sems.(k).Joint  = sem_over_iter(PID_v(:,1,:));

    Row.means.(k).Syn   = mean_over_iter(PID_v(:,2,:));
    Row.sems.(k).Syn    = sem_over_iter(PID_v(:,2,:));

    Row.means.(k).Red   = mean_over_iter(PID_v(:,3,:));
    Row.sems.(k).Red    = sem_over_iter(PID_v(:,3,:));

    if any(strcmpi(cfg.info_types,'Unq'))
        UNQ = squeeze(mean(PID_v(:,4:5,:), 2)); % T x K
        Row.means.(k).Unq = mean(UNQ, 2, 'omitnan')';
        Row.sems.(k).Unq  = sem_matrix(UNQ);
    end
end
end


function Row = load_row_gauss(cfg, simul_case)
% Loads the Gauss file for this row and aggregates biases.
Row.means = struct(); Row.sems = struct();

% Choose file by case
switch lower(strrep(simul_case,' ','_'))
    case 'high_syn'
        filename = fullfile('Results','Gauss','Finalresults_across_M_and_ntrials_high_synergy.mat');
    case 'zero_syn'
        filename = fullfile('Results','Gauss','Finalresults_across_M_and_ntrials_zero_synergy.mat');
    case 'both_unq'
        filename = fullfile('Results','Gauss','Finalresults_across_M_and_ntrials_both_unique.mat');
    case {'bit_of_all','bitofall'}
        filename = fullfile('Results','Gauss','Finalresults_across_M_and_ntrials_bit_of_all.mat');
    otherwise
        filename = fullfile('Results','Gauss','Finalresults_across_M_and_ntrials.mat');
end

L = load(filename); % expects sampled_results, GT_results, M_vals, ntrials_vals
PID_all = L.sampled_results;  % [comp x trials x M x bias x info x iter]
GT_all  = L.GT_results;       % [comp x trials x M x bias x info]
M_vals  = L.M_vals;

% Row-specific X vector from file
Row.trials = L.ntrials_vals(:);

% Dimension selection
if ~isfield(cfg,'gauss_dim') || isempty(cfg.gauss_dim), cfg.gauss_dim = 20; end
dimidx = find(M_vals == cfg.gauss_dim, 1, 'first');
if isempty(dimidx), error('Requested gauss_dim=%d not found.', cfg.gauss_dim); end

% Ground truth (plugin slice; GT identical across bias)
GT_sel = squeeze(GT_all(:, :, dimidx, 1, cfg.info_amount));  % [comp x trials]

% Component index mapping (two variants exist)
if size(PID_all,1) == 13
    map.Joint = 1; map.Syn = 6; map.Red = 5; map.Unq = 2; use_half_for_unq = true;
else
    map.Joint = 1; map.Syn = 7; map.Red = 6; map.Unq = [2 3]; use_half_for_unq = false;
end

% GT map (vector per trials)
GT.Joint = GT_sel(map.Joint, :);
GT.Syn   = GT_sel(map.Syn,   :);
GT.Red   = GT_sel(map.Red,   :);
if use_half_for_unq, GT.Unq = GT_sel(map.Unq,:)/2; else, GT.Unq = mean(GT_sel(map.Unq,:),1); end
Row.GT = GT;

% ---- 1) compute base biases first (no hyphens in field names) ----
all_b = lower(cfg.bias_correction);
base_list = setdiff(all_b, {'shuff-resample'});  % compute these directly

for b = 1:numel(base_list)
    bc = base_list{b};
    PID_b = select_bias_slice_gauss(PID_all, dimidx, cfg.info_amount, bc); % [comp x trials x iter]
    k = bias_key(bc);

    Row.means.(k) = struct(); Row.sems.(k) = struct();

    Row.means.(k).Joint = mean_over_iter(PID_b(map.Joint,:,:));
    Row.sems.(k).Joint  = sem_over_iter(PID_b(map.Joint,:,:));

    Row.means.(k).Syn   = mean_over_iter(PID_b(map.Syn,:,:));
    Row.sems.(k).Syn    = sem_over_iter(PID_b(map.Syn,:,:));

    Row.means.(k).Red   = mean_over_iter(PID_b(map.Red,:,:));
    Row.sems.(k).Red    = sem_over_iter(PID_b(map.Red,:,:));

    if any(strcmpi(cfg.info_types,'Unq'))
        if use_half_for_unq
            A = squeeze(PID_b(map.Unq,:,:))/2; % [trials x iter]
        else
            A = squeeze(mean(PID_b(map.Unq,:,:),1)); % avg U1,U2 -> [trials x iter]
        end
        Row.means.(k).Unq = mean(A, 2, 'omitnan')';
        Row.sems.(k).Unq  = sem_matrix(A);
    end
end

% ---- 2) if user asked for 'shuff-resample', combine corrected values ----
if any(strcmpi(all_b, 'shuff-resample'))
    if isfield(Row.means,'resample') && isfield(Row.means,'shuff')
        k = 'shuff_resample';   % sanitized field
        Row.means.(k) = struct(); Row.sems.(k) = struct();

        fields = fieldnames(Row.means.resample);
        for f = 1:numel(fields)
            nm  = fields{f};
            muR = Row.means.resample.(nm);
            muS = Row.means.shuff.(nm);
            seR = Row.sems.resample.(nm);
            seS = Row.sems.shuff.(nm);

            % Mean of corrected curves
            Row.means.(k).(nm) = 0.5*(muR + muS);

            % SEM for average of two independent estimates:
            % var(avg) = (varR + varS)/4  -> sem = sqrt(seR^2 + seS^2)/2
            Row.sems.(k).(nm)  = 0.5*sqrt(seR.^2 + seS.^2);
        end
    else
        warning('Requested ''shuff-resample'' but base ''shuff''/''resample'' missing.');
    end
end
end

%% ---------------------- Helpers: math & slices ---------------------------

function mu = mean_over_iter(A) % A: [T x 1 x K] or [1 x T x K] or [1 x 1 x K]
    B = squeeze(A);              % -> [T x K] or [T]
    if isvector(B), mu = B(:)'; else, mu = mean(B, 2, 'omitnan')'; end
end

function s = sem_over_iter(A) % like above, returns [1 x T]
    B = squeeze(A);                        % [T x K] or [T]
    if isvector(B), s = zeros(1,numel(B)); return; end
    s = std(B, 0, 2, 'omitnan')' ./ sqrt(max(1,sum(~isnan(B),2))');
end

function s = sem_matrix(B) % B: [T x K]
    s = std(B, 0, 2, 'omitnan')' ./ sqrt(max(1,sum(~isnan(B),2))');
end

function xticks_like(ax, xvals)
%XTICKS_LIKE Set decade ticks/labels on a log-scaled x-axis.
%   XTICKS_LIKE(AX, XVALS) chooses ticks at powers of 10 spanning XVALS
%   (positive values only), sets them on axes AX, and labels them as 10^k.

    if nargin < 2 || ~ishandle(ax) || ~isa(ax, 'matlab.graphics.axis.Axes')
        error('First input must be an axes handle, second a numeric vector of x-values.');
    end
    xvals = xvals(:);
    xvals = xvals(isfinite(xvals) & xvals > 0); % log scale requires positive values
    if isempty(xvals), return; end

    xmin = min(xvals); xmax = max(xvals);
    eMin = floor(log10(xmin)); eMax = ceil(log10(xmax));
    exps = eMin:eMax;
    ticks = 10.^exps;

    inRange = ticks >= xmin & ticks <= xmax;
    ticks = ticks(inRange);
    exps  = exps(inRange);
    if isempty(ticks), return; end

    set(ax, 'XScale', 'log');
    xticks(ax, ticks);
    lbl = arrayfun(@(k) sprintf('10^{%d}', k), exps, 'UniformOutput', false);
    xticklabels(ax, lbl);
end

function name = pretty_label(bc)
    switch lower(bc)
        case 'qe',               name = 'QE';
        case 'shuff',            name = 'shuff-sub';
        case 'qeshuff',          name = 'merged';
        case 'infocorr',         name = 'resampling';
        case 'shuff-resample',   name = 'shuff-resample';
        otherwise,               name = bc;
    end
end

function PID_b = select_bias_slice_gauss(PID_all, dimidx, info_amount, bias_name)
% Returns [comp x trials x iter] for a single base bias.
    switch lower(bias_name)
        case 'plugin',    bidx = 1;
        case 'resample',  bidx = 2;
        case 'shuff',     bidx = 3;
        case 'venkatesh', bidx = 4;
        otherwise
            error('Unknown base bias for Gauss: %s', bias_name);
    end
    PID_b = squeeze(PID_all(:,:,dimidx,bidx,info_amount,:));
end

function key = bias_key(bc)
% Lowercase and make a valid struct field (replaces '-' with '_', etc.)
    key = matlab.lang.makeValidName(lower(bc));
end