function pid_plot_grid(cfg)

methodColors = containers.Map( ...
    {'plugin','qe','resample','resampling','shuff','shuffle','qe_shuff','qeshuff', ...
     'merged','shuff-resamp','shuff-resample','infocorr','weighted','venkatesh','Venkatesh'}, ...
    { [0.60 0.60 0.60], ...
      [0.00 0.447 0.741], ...
      [0.929 0.694 0.125], ...
      [0.929 0.694 0.125], ...
      [0.85  0.325 0.098], ...
      [0.85  0.325 0.098], ...
      [0.494 0.184 0.556], ...
      [0.494 0.184 0.556], ...
      [0.494 0.184 0.556], ...
      [0.494 0.184 0.556], ...
      [0.494 0.184 0.556], ...
      [0.929 0.694 0.125], ...
      [0.60  0.60  0.60], ...
      [0.466 0.674 0.188], ...
      [0.466 0.674 0.188] } );

if ~isfield(cfg,'mode') || isempty(cfg.mode)
    cfg.mode = 'discrete';
end
if ~isfield(cfg,'trial_categories') || isempty(cfg.trial_categories)
    cfg.trial_categories = [16 32 64 128 256 512 1024 2048];
end

cfg.info_types = {'Joint','Synergy','Redundancy'};

if ~isfield(cfg,'row_names') || isempty(cfg.row_names)
    cfg.row_names = cfg.simul_cases;
end

rm = lower(string(cfg.redundancy_measure));
use_asymptotic_gt = strcmpi(cfg.mode,'discrete') && ismember(rm, ["imin","i_min","immi","i_mmi", 'i_broja']);

set(0,'DefaultTextFontName','Arial');
set(0,'DefaultAxesFontName','Arial');
set(0,'DefaultLegendFontName','Arial');
set(0,'DefaultAxesFontSize',8);

nrows = numel(cfg.simul_cases);
ncols = numel(cfg.info_types);

figure('Units','centimeters','Position',[1 1 21 2.6*nrows*1.5]);
tl = tiledlayout(nrows,ncols,'TileSpacing','compact','Padding','compact');

for r = 1:nrows
    simul_case = cfg.simul_cases{r};

    if strcmpi(cfg.mode,'gauss')
        Row = load_row_gauss(cfg, simul_case);
    else
        Row = load_row_discrete(cfg, simul_case);
    end

    x = Row.trials(:)';


    for c = 1:ncols
        ax = nexttile; hold(ax,'on');

        info_label = cfg.info_types{c};
        key = switch_key(info_label);

        L = gobjects(numel(cfg.bias_correction),1);

        for b = 1:numel(cfg.bias_correction)
            bc = cfg.bias_correction{b};
            k  = bias_key(bc);

            if ~isfield(Row.means,k), continue; end
            if ~isfield(Row.means.(k),key), continue; end

            mu = Row.means.(k).(key);
            se = Row.sems.(k).(key);

            n = min([numel(x),numel(mu),numel(se)]);
            if n==0, continue; end

            xx = x(1:n);
            mu = mu(1:n);
            se = se(1:n);

            clr = methodColors(lower(bc));

            fill([xx fliplr(xx)], ...
                 [mu+2*se fliplr(mu-2*se)], ...
                 clr,'FaceAlpha',0.15,'EdgeColor','none', ...
                 'HandleVisibility','off');

            L(b) = plot(xx,mu,'LineWidth',1.5, ...
                'Color',clr,'DisplayName',pretty_label(bc));
        end

        if strcmpi(cfg.mode,'gauss')
            if isfield(Row,'GT') && isfield(Row.GT,key)
                gt = Row.GT.(key);
                ng = min(numel(x),numel(gt));
                plot(x(1:ng),gt(1:ng),'k--','LineWidth',1.2,'HandleVisibility','off');
            end
        else
            if use_asymptotic_gt
                if isfield(Row.means,'qeshuff') && isfield(Row.means.qeshuff,key)
                    gt_curve = Row.means.qeshuff.(key)(:);
                    gt_val = mean(gt_curve(max(end-1,1):end),'omitnan');
                    plot([x(1) x(end)],[gt_val gt_val],'k--','LineWidth',1.2,'HandleVisibility','off');
                end
            elseif isfield(Row,'GT') && isfield(Row.GT,key)
                plot([x(1) x(end)],[Row.GT.(key) Row.GT.(key)],'k--','LineWidth',1.2,'HandleVisibility','off');
            end
        end

        if strcmpi(cfg.mode,'gauss') && isfield(cfg,'gauss_dim')
            xline(12*cfg.gauss_dim,'k--','LineWidth',1.2,'HandleVisibility','off');
        end

        yline(0,'k','LineWidth',1,'HandleVisibility','off');
        set(ax,'XScale','log','LineWidth',1.1);
        xlim([x(1) x(end)]);

        % Rule-of-thumb vertical line
        if isfield(cfg,'thumb_rule_x') && ~isempty(cfg.thumb_rule_x) && ~strcmpi(cfg.mode,'gauss')
            xr = cfg.thumb_rule_x;
            if xr >= x(1) && xr <= x(end)  % avoid drawing outside x-limits
                xline(ax, xr, 'k--', 'LineWidth', 1.2, 'HandleVisibility', 'off');
            end
        end


        if r==nrows
            xlabel('Trials (N)');
        else
            ax.XTickLabel = [];
        end

        if c==1
            ylabel(sprintf('%sInformation [bits]',cfg.row_names{r}));
        else
            ax.YTickLabel = [];
        end

        if isfield(cfg,'row_ylim') && ~isempty(cfg.row_ylim)
            ylim(cfg.row_ylim(r,:));
        end

        if r==1
            title(info_label,'FontWeight','bold');
        end

        if r==1 && c==ncols
            legend(ax,L(isgraphics(L)),'Location','northwest','Box','off');
        end

        hold(ax,'off');
    end
end

if isfield(cfg,'save_path') && ~isempty(cfg.save_path)
    outdir = fileparts(cfg.save_path);
    if ~exist(outdir,'dir'), mkdir(outdir); end
    exportgraphics(tl,cfg.save_path,'ContentType','vector');
end
end

% =======================================================================
% ========================= LOADERS =====================================
% =======================================================================

function Row = load_row_discrete(cfg,simul_case)

Row.means = struct(); Row.sems = struct();
Row.trials = cfg.trial_categories(:)' * 4;

datafile = fullfile('Results',cfg.ResultFolder, ...
    sprintf('Simuldata_%s_%s_%s.mat', simul_case, cfg.info_amount, cfg.redundancy_measure));
D = load(datafile);

for b = 1:numel(cfg.bias_correction)
    bc = cfg.bias_correction{b};
    v  = ['PID_v_' bc];
    if ~isfield(D,v), continue; end
    V = D.(v);
    k = bias_key(bc);

    Row.means.(k).Joint = mean(squeeze(V(:,1,:)),2,'omitnan')';
    Row.means.(k).Syn   = mean(squeeze(V(:,2,:)),2,'omitnan')';
    Row.means.(k).Red   = mean(squeeze(V(:,3,:)),2,'omitnan')';

    Row.sems.(k).Joint = sem(squeeze(V(:,1,:)));
    Row.sems.(k).Syn   = sem(squeeze(V(:,2,:)));
    Row.sems.(k).Red   = sem(squeeze(V(:,3,:)));
end
end

function Row = load_row_gauss(cfg, simul_case)

Row.means = struct(); Row.sems = struct();

switch lower(strrep(simul_case,' ','_'))
    case 'high_syn'
        fname = 'Finalresults_across_M_and_ntrials_high_synergy.mat';
    case 'zero_syn'
        fname = 'Finalresults_across_M_and_ntrials_zero_synergy.mat';
    case 'both_unq'
        fname = 'Finalresults_across_M_and_ntrials_both_unique.mat';
    case {'bit_of_all','bitofall'}
        fname = 'Finalresults_across_M_and_ntrials_bit_of_all.mat';
    otherwise
        fname = 'Finalresults_across_M_and_ntrials.mat';
end

S = load(fullfile('Results','Gauss',fname));
PID = S.sampled_results;
Row.trials = S.ntrials_vals(:)';

mid = find(S.M_vals == cfg.gauss_dim,1);

if size(PID,1)==13
    map.Joint=1; map.Red=5; map.Syn=6;
else
    map.Joint=1; map.Red=6; map.Syn=7;
end

% ---------- GROUND TRUTH (no repetitions) ----------
if isfield(S,'GT_results') && ~isempty(S.GT_results)
    GT = S.GT_results;   % expected: PID dims without repetitions
    Row.GT = struct();
    Row.GT.Joint = gt_curve_from_GT(GT, map.Joint, mid, cfg.info_amount, PID);
    Row.GT.Red   = gt_curve_from_GT(GT, map.Red,   mid, cfg.info_amount, PID);
    Row.GT.Syn   = gt_curve_from_GT(GT, map.Syn,   mid, cfg.info_amount, PID);
end

for b = 1:numel(cfg.bias_correction)
    bc = lower(cfg.bias_correction{b});
    switch bc
        case 'plugin',    bidx=1;
        case 'resample',  bidx=2;
        case 'shuff',     bidx=3;
        case 'venkatesh', bidx=4;
        otherwise, continue;
    end

    X = squeeze(PID(:,:,mid,bidx,cfg.info_amount,:));
    k = bias_key(bc);

    Row.means.(k).Joint = mean(squeeze(X(map.Joint,:,:)),2,'omitnan')';
    Row.means.(k).Red   = mean(squeeze(X(map.Red,:,:)),2,'omitnan')';
    Row.means.(k).Syn   = mean(squeeze(X(map.Syn,:,:)),2,'omitnan')';

    Row.sems.(k).Joint = sem(squeeze(X(map.Joint,:,:)));
    Row.sems.(k).Red   = sem(squeeze(X(map.Red,:,:)));
    Row.sems.(k).Syn   = sem(squeeze(X(map.Syn,:,:)));
end

% ---------- MERGED (resample + shuffle) ----------
if isfield(Row.means,'resample') && isfield(Row.means,'shuff')

    km = bias_key('merged');

    Row.means.(km).Joint = 0.5*(Row.means.resample.Joint + Row.means.shuff.Joint);
    Row.means.(km).Red   = 0.5*(Row.means.resample.Red   + Row.means.shuff.Red);
    Row.means.(km).Syn   = 0.5*(Row.means.resample.Syn   + Row.means.shuff.Syn);

    Row.sems.(km).Joint = 0.5*sqrt(Row.sems.resample.Joint.^2 + Row.sems.shuff.Joint.^2);
    Row.sems.(km).Red   = 0.5*sqrt(Row.sems.resample.Red.^2   + Row.sems.shuff.Red.^2);
    Row.sems.(km).Syn   = 0.5*sqrt(Row.sems.resample.Syn.^2   + Row.sems.shuff.Syn.^2);

    % alias so shuff-resample works
    Row.means.(bias_key('shuff-resample')) = Row.means.(km);
    Row.sems.(bias_key('shuff-resample'))  = Row.sems.(km);
end

end

function v = gt_curve_from_GT(GT, compIdx, mid, ia, PID)
% Extract a 1×ntrials ground-truth curve for one PID component.
% GT is expected to be PID without the repetitions dimension (5D vs 6D),
% but we handle common dimension-order variants robustly.

nd = ndims(GT);
sz = size(GT);

% Default empty
v = [];

if nd == 5
    % Most common: (comp, ntrials, M, bidx, info_amount)
    if numel(sz) >= 5 && sz(4) == size(PID,4) && sz(5) == size(PID,5)
        v = squeeze(GT(compIdx,:,mid,1,ia)); % bidx not meaningful for GT -> pick 1
    % Swapped last dims: (comp, ntrials, M, info_amount, bidx)
    elseif numel(sz) >= 5 && sz(5) == size(PID,4) && sz(4) == size(PID,5)
        v = squeeze(GT(compIdx,:,mid,ia,1));
    % Fallback: average over remaining trailing dims
    else
        v = squeeze(mean(GT(compIdx,:,mid,:,:), [4 5], 'omitnan'));
    end

elseif nd == 4
    % Possible: (comp, ntrials, M, info_amount)
    if sz(4) >= ia
        v = squeeze(GT(compIdx,:,mid,ia));
    else
        v = squeeze(GT(compIdx,:,mid,1));
    end

else
    % Last resort: squeeze whatever matches
    v = squeeze(GT(compIdx,:,mid));
end

v = v(:)'; % enforce row vector
end


function key = switch_key(label)
switch label
    case 'Synergy',    key='Syn';
    case 'Redundancy', key='Red';
    otherwise,         key='Joint';
end
end

function s = sem(X)
X = squeeze(X);
s = std(X,0,2,'omitnan')' ./ sqrt(size(X,2));
end

function k = bias_key(bc)
k = matlab.lang.makeValidName(lower(bc));
end

function name = pretty_label(bc)
switch lower(bc)
    case {'merged','shuff-resample','shuff-resamp','qeshuff'}
        name = 'merged';
    case 'shuff'
        name = 'shuffle';
    case 'infocorr'
        name = 'resample';
    otherwise
        name = bc;
end
end
