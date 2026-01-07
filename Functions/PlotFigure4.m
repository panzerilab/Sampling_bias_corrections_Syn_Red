function PlotFigure4(trial_categories, simul_case, info_amount, redundancy_measure, ...
    ResultFolder, maxvals, datatype)
% PlotFigure4
% Panels: Joint | Syn | Red
% X-axis uses Trials (N).
%
% Fixes:
%  - Gaussian case uses exact trial matching (no nearest-neighbor duplicates)
%  - Robust Gaussian component mapping and GT extraction
%  - Colors match reference; resample/resampling forced to YELLOW
%  - Correct rule-of-thumb vertical line:
%      * Discrete:  N = 4*nRbins^3, nRbins=4 -> 256
%      * Gaussian:  N = 4*3*d, d=20         -> 240

load('adaptive_weight_matrix_allcases.mat', 'weight_matrix', 'info_levels_fine');

% -------------------- Visual defaults --------------------
set(0,'DefaultTextFontName','Arial');
set(0,'DefaultAxesFontName','Arial');
set(0,'DefaultLegendFontName','Arial');

FontSize = 9;
legendFontSize = 8;
lw = 1;
n_iter = 100;

% Titles exactly as you want (plot titles)
comp_titles = {'Joint','Synergy','Redundancy'};

% -------------------- X AXIS (Trials N) --------------------
% Discrete convention is trial_categories*4. Use the same for Gaussian so
% both figures are directly comparable and match your reference.
xvals = trial_categories(:)' * 4;

% -------------------- Method colors (reference) --------------------
% Match the reference figure, except resample/resampling is YELLOW.
methodColors = containers.Map( ...
    {'plugin','qe','resample','resampling','shuff','shuffle','qe_shuff','merged','shuff-resamp','weighted','venkatesh','Venkatesh'}, ...
    { [0.60 0.60 0.60], ...         % plugin (gray) [not used here but safe]
      [0.00 0.447 0.741], ...       % QE (blue)
      [0.929 0.694 0.125], ...      % resample (YELLOW)  <-- requested change
      [0.929 0.694 0.125], ...      % resampling (YELLOW)
      [0.85 0.325 0.098], ...       % shuff (orange)
      [0.85 0.325 0.098], ...       % shuffle (orange)
      [0.494 0.184 0.556], ...      % qe_shuff (purple)
      [0.494 0.184 0.556], ...      % merged (purple)
      [0.494 0.184 0.556], ...      % shuff-resamp (purple)
      [0.60 0.60 0.60], ...         % weighted (gray)
      [0.466 0.674 0.188], ...      % venkatesh (green)
      [0.466 0.674 0.188] } );      % Venkatesh (green)

series = struct('name',{},'color',{},'mean',{},'sem',{});

% ======================= DISCRETE ========================
if strcmpi(datatype,'discr')

    data_fname = sprintf('Results/%s/Simuldata_%s_%s_%s.mat', ...
        ResultFolder, simul_case, info_amount, redundancy_measure);

    gt_fname = sprintf('Results/%s/GroundTruth_%s_%s.mat', ...
        ResultFolder, simul_case, info_amount);

    S = load(gt_fname,'GroundTruth_value');
    GT = S.GroundTruth_value;

    GT_joint = sum(GT(1:4));
    GT_syn   = GT(4);
    GT_red   = GT(1);
    gt_scalar = [GT_joint, GT_syn, GT_red];

    % Rule-of-thumb: N = 4*nRbins^3 with nRbins=4 -> 256
    thumbRule = 4 * 4^3;

    % Methods to try (as in your code)
    try_bias = {'qe','shuff','qe_shuff','weighted','venkatesh'};
    ser_ix = 0;

    for bi = 1:numel(try_bias)
        bc = try_bias{bi};
        try
            if strcmp(bc,'weighted')
                L1 = load(data_fname,'PID_v_qe');
                L2 = load(data_fname,'PID_v_shuff');

                PID_qe    = L1.PID_v_qe;
                PID_shuff = L2.PID_v_shuff;

                max_info = log2(4);
                PID_joint = (PID_qe(:,1,:) + PID_shuff(:,1,:))/2;
                info_frac = squeeze(mean(PID_joint,3)) / max_info;

                PID_weighted = zeros(size(PID_qe));
                for i = 1:size(PID_qe,1)
                    [~,idx] = min(abs(info_levels_fine - info_frac(i)));
                    w = weight_matrix(3,i,idx);
                    PID_weighted(i,:,:) = w * PID_qe(i,:,:) + (1-w) * PID_shuff(i,:,:);
                end
                PID_v = PID_weighted;
            else
                varname = sprintf('PID_v_%s', bc);
                L = load(data_fname,varname);
                PID_v = L.(varname);
            end

            % Means/SEM (3*SEM) across iters
            m = zeros(numel(xvals),3);
            s = zeros(numel(xvals),3);
            for c = 1:3
                vals = squeeze(PID_v(:,c,:));
                m(:,c) = mean(vals,2,'omitnan');
                s(:,c) = 3 * std(vals,0,2,'omitnan') / sqrt(n_iter);
            end

            ser_ix = ser_ix + 1;
            series(ser_ix).name  = canonical_bias_name(bc);
            series(ser_ix).color = methodColors(color_key_for_method(bc));
            series(ser_ix).mean  = m;
            series(ser_ix).sem   = s;
        catch
            % silently skip missing methods, as your original code intended
        end
    end

% ======================= GAUSSIAN ========================
else
    filename = sprintf('Results/%s/Finalresults_across_M_and_ntrials.mat', ResultFolder);
    D = load(filename);

    PID_all = D.sampled_results;   % expected: [comp x trials x M x bias x info x iter] (common)
    GT_all  = D.GT_results;        % expected: [comp x trials x M x bias x info] (common)
    M_vals  = D.M_vals;
    ntrials_vals = D.ntrials_vals(:)';

    d = 20;
    dimidx = find(M_vals==d,1);
    if isempty(dimidx)
        error('M=20 not found.');
    end

    % Rule-of-thumb: N = 4*3*d with d=20 -> 240
    thumbRule = 4 * 3 * d;

    % ---- Exact x matching (NO nearest neighbor) ----
    [tf, tidx] = ismember(xvals, ntrials_vals);
    if ~all(tf)
        missing = xvals(~tf);
        warning('Dropping %d x-points not present in Gaussian file ntrials_vals (e.g. %s ...).', ...
            numel(missing), mat2str(missing(1:min(5,end))));
    end
    xvals_g = xvals(tf);
    tidx_g  = tidx(tf);

    % Methods in Gaussian results
    bias_types      = {'resample','shuffle','shuff-resamp','venkatesh'};
    bias_type_names = {'resampling','shuffle','merged','Venkatesh'}; % display names
    ser_ix = 0;

    % Determine component mapping once (robust to format)
    compMap = gauss_component_map(size(PID_all,1));

    for bi = 1:numel(bias_types)
        try
            [PID_v, GT_v] = local_pick_gauss_trials( ...
                PID_all, GT_all, tidx_g, dimidx, bias_types{bi}, info_amount);

            % PID_v: [comp x Tsel x iter]
            % GT_v : [comp x Tsel]
            m = [
                squeeze(mean(PID_v(compMap.Joint,:,:),3))', ...
                squeeze(mean(PID_v(compMap.Syn  ,:,:),3))', ...
                squeeze(mean(PID_v(compMap.Red  ,:,:),3))'
            ];

            s = [
                3*std(squeeze(PID_v(compMap.Joint,:,:)),0,2)/sqrt(n_iter), ...
                3*std(squeeze(PID_v(compMap.Syn  ,:,:)),0,2)/sqrt(n_iter), ...
                3*std(squeeze(PID_v(compMap.Red  ,:,:)),0,2)/sqrt(n_iter)
            ];

            ser_ix = ser_ix + 1;
            series(ser_ix).name  = bias_type_names{bi};
            series(ser_ix).color = methodColors(color_key_for_method(bias_types{bi}));
            series(ser_ix).mean  = m;
            series(ser_ix).sem   = s;

            GT_joint_vec = GT_v(compMap.Joint,:);
            GT_syn_vec   = GT_v(compMap.Syn  ,:);
            GT_red_vec   = GT_v(compMap.Red  ,:);
        catch
            % skip missing methods gracefully
        end
    end

    % Use Gaussian-matched x-values for plotting
    xvals = xvals_g;
end

% ======================= PLOTTING ========================
figure('Units','centimeters','Position',[1 1 16.5 6]);
tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

for c = 1:3
    nexttile; hold on;

    for k = 1:numel(series)

        if strcmpi(datatype,'discr')
            gt = gt_scalar(c);
        else
            gtmat = [GT_joint_vec; GT_syn_vec; GT_red_vec];
            gt = gtmat(c,:);
        end

        m = series(k).mean(:,c) - gt(:);
        s = series(k).sem(:,c);

        % robust fill
        xx = xvals(:)'; m = m(:)'; s = s(:)';
        Xf = [xx, fliplr(xx)];
        Yf = [m+s, fliplr(m-s)];

        hfill = fill(Xf, Yf, series(k).color, 'FaceAlpha',0.18,'EdgeColor','none');
        hfill.Annotation.LegendInformation.IconDisplayStyle = 'off';

        plot(xvals, m, 'LineWidth', lw, 'Color', series(k).color);
    end

    set(gca,'XScale','log','FontSize',FontSize,'LineWidth',lw);
    xlim([min(xvals) max(xvals)]);

    xticks([1e2 1e3]);
    xticklabels({'10^2','10^3'});

    yline(0,'k-','LineWidth',lw);
    xline(thumbRule,':k','LineWidth',lw);

    title(comp_titles{c},'FontWeight','bold');
    xlabel('Trials (N)');

    if c==1
        ylabel('Residual bias [bits]');
        legend({series.name},'Box','off','FontSize',legendFontSize,'Location','northwest');
    end

    hold off;
end

end % PlotFigure4


% ====================== HELPERS ==========================
function [PID_v, GT_v] = local_pick_gauss_trials(PID_all, GT_all, trial_idx_vec, dimidx, bias_name, info_idx)
% Returns:
%   PID_v: [comp x Tsel x iter]
%   GT_v : [comp x Tsel]
%
% Assumes common Gaussian layout:
%   PID_all: [comp x trials x M x bias x info x iter]
%   GT_all : [comp x trials x M x bias x info]
%
% This matches your working compute_gauss_means_for_method_weighted() logic.

switch lower(bias_name)
    case 'plugin',    b = 1;
    case 'resample',  b = 2;
    case 'shuffle',   b = 3;
    case 'venkatesh', b = 4;
    case 'shuff-resamp'
        PID_r = squeeze(PID_all(:, trial_idx_vec, dimidx, 2, info_idx, :));
        PID_s = squeeze(PID_all(:, trial_idx_vec, dimidx, 3, info_idx, :));
        PID_v = 0.5*PID_r + 0.5*PID_s;

        GT_v = squeeze(GT_all(:, trial_idx_vec, dimidx, 1, info_idx));
        return
    otherwise
        error('Unknown bias: %s', bias_name);
end

PID_v = squeeze(PID_all(:, trial_idx_vec, dimidx, b, info_idx, :));
GT_v  = squeeze(GT_all(:,  trial_idx_vec, dimidx, 1, info_idx));
end


function map = gauss_component_map(nComp)
% Robust mapping for the two common Gaussian result formats.
% nComp==13 is the "older" indexing seen in some pipelines.
if nComp == 13
    map.Joint = 1;
    map.Syn   = 6;
    map.Red   = 5;
else
    map.Joint = 1;
    map.Syn   = 7;
    map.Red   = 6;
end
end


function name = canonical_bias_name(b)
b = lower(string(b));
switch b
    case 'qe'
        name = 'QE';
    case {'shuff','shuffle'}
        name = 'shuffle';
    case {'qe_shuff','shuff-resamp','shuff_resamp'}
        name = 'merged';
    case {'resample','resampling'}
        name = 'resampling';
    case 'venkatesh'
        name = 'Venkatesh';
    otherwise
        name = char(b);
end
end


function key = color_key_for_method(b)
% Normalize method keys to those used in methodColors.
b = lower(string(b));
switch b
    case 'qe'
        key = 'qe';
    case {'shuff','shuffle'}
        key = 'shuffle';
    case {'qe_shuff','shuff-resamp','shuff_resamp'}
        key = 'merged';
    case {'resample','resampling'}
        key = 'resampling';
    case 'weighted'
        key = 'weighted';
    case 'venkatesh'
        key = 'venkatesh';
    otherwise
        key = char(b);
end
end
