function [worst_all, worst_suff_all] = PlotFigure5( ...
    trialIdx, categories, simul_case, bias_types, info_amount, redundancy_measure, ...
    subplot_idx, fig_handle, rows, cols, ResultFolder, min_val, max_val, isGauss, isUnion, addInset, insetRange)

if nargin < 16 || isempty(addInset)
    addInset = false;
end
if nargin < 17
    insetRange = []; % auto per-backend
end

% Visual defaults (unchanged)
set(0,'DefaultTextFontName','Arial');
set(0,'DefaultAxesFontName','Arial');
set(0,'DefaultLegendFontName','Arial');

FontSize = 8;
legendFontSize = 8;
titleFontWeight = 'bold';

comp_names = {'Joint','Syn','Red'};
n_iter = 100;
lw = 1;

% ===================== COLOR MAP (USER-SPECIFIED) ======================
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

normalize_name = @(s) lower(strrep(strrep(strrep(s,'-','_'),' ','_'), '__','_'));

get_method_color = @(rawName) local_get_method_color(methodColors, rawName, normalize_name);

% ======================================================================
% ======================= DISCRETE (Alpha_Bin_sweep) ===================
% ======================================================================
if ~isGauss && ~isUnion

    data_fname = sprintf('Results/%s/Simuldata_%s_%s_%s.mat', ResultFolder, simul_case, info_amount, redundancy_measure);
    gt_fname   = sprintf('Results/%s/GroundTruth_%s_%s.mat', ResultFolder, simul_case, info_amount);

    % Ground truth from qeShuff at last trial (unchanged)
    Lgt = load(gt_fname, 'GroundTruth_value'); % trials x comps x iters x cats
    PID_gt_all = Lgt.GroundTruth_value;
    GT_all = vertcat(PID_gt_all{:});           % comps x iters x cats

    GT_joint = sum(GT_all(:,1:4),2);
    GT_red   = GT_all(:,1);
    GT_syn   = GT_all(:,4);
    GT_UNQ   = GT_all(:,2:3);
    GT_unq   = mean(GT_UNQ,2);
    GT_vec   = [GT_joint, GT_syn, GT_red];

    % thumb rule (bins) (unchanged)
    if trialIdx == 1
        N=256;
    elseif trialIdx==2
        N=512;
    elseif trialIdx==3
        N=1024;
    else
        N=512;
    end
    thumbRule = floor(nthroot(N/(4*4), 3));

    % Load all bias series
    series = struct('name',{},'color',{},'mean',{},'sem',{});

    for bi = 1:numel(bias_types)
        bname_raw = bias_types{bi};              % keep original label
        bname_norm = normalize_name(bname_raw);  % for matching/loading

        try
            % ---- Weighted (your existing logic kept) ----
            if strcmpi(bname_norm,'weighted')
                PID_qe    = local_load_discrete_slice(data_fname, 'qe', trialIdx);
                PID_shuff = local_load_discrete_slice(data_fname, 'shuff', trialIdx);

                max_info = log2(2:9);
                PID_qeshuff = (PID_qe(:,1,:)+PID_shuff(:,1,:))/2;
                infoval_Joint = squeeze(mean(PID_qeshuff,3));
                info_frac_matrix = infoval_Joint' ./ max_info;

                S = load('adaptive_weight_matrix_allcases.mat', 'weight_matrix', 'info_levels_fine');
                weight_matrix = S.weight_matrix;
                info_levels_fine = S.info_levels_fine;

                PID_weighted = zeros(size(PID_shuff));
                for i = 1:size(PID_shuff, 1)
                    val = info_frac_matrix(i);
                    [~, idx] = min(abs(info_levels_fine - val));
                    w = weight_matrix(i,4,idx);
                    QE_slice = squeeze(PID_qe(i, :, :));
                    Shuff_slice = squeeze(PID_shuff(i, :, :));
                    PID_weighted(i, :, :) = w * QE_slice + (1 - w) * Shuff_slice;
                end
                PID_v = PID_weighted;

            else
                % ---- Standard methods ----
                PID_v = local_load_discrete_slice(data_fname, bname_norm, trialIdx);
            end

            UNQ     = squeeze(mean(PID_v(:,4:5,:), 2)); % cats x iters
            m_joint = (squeeze(mean(PID_v(:,1,:), 3, 'omitnan')) - GT_joint);
            m_red   = (squeeze(mean(PID_v(:,3,:), 3, 'omitnan')) - GT_red);
            m_unq   = (squeeze(mean(UNQ,          3, 'omitnan')) - GT_unq);
            m_syn   = (squeeze(mean(PID_v(:,2,:), 3, 'omitnan')) - GT_syn);

            s_joint = 3 * std(squeeze(PID_v(:,1,:)), 0, 2, 'omitnan') / sqrt(n_iter);
            s_red   = 3 * std(squeeze(PID_v(:,3,:)), 0, 2, 'omitnan') / sqrt(n_iter);
            s_unq   = 3 * std(UNQ,                   0, 2, 'omitnan') / sqrt(n_iter);
            s_syn   = 3 * std(squeeze(PID_v(:,2,:)), 0, 2, 'omitnan') / sqrt(n_iter);

            series(end+1).name  = bname_raw; %#ok<AGROW>
            series(end).color   = get_method_color(bname_raw); % <<< ONLY COLOR SOURCE
            series(end).mean    = [m_joint, m_syn, m_red];
            series(end).sem     = [s_joint, s_syn, s_red];

        catch ME
            warning('Skipping bias "%s": %s', bname_raw, ME.message);
        end
    end

    x_vals = categories;

    if isempty(insetRange)
        insetRange = [2 5]; % bins 5..9 for discrete (unchanged behavior)
    end

    x_label = 'Number of bins';
    cat_sufficient = 5;
    y_label = 'Bias [bits]';

% ======================================================================
% ======================== GAUSSIAN (dimension sweep) ==================
% ======================================================================
elseif isGauss && ~isUnion

    if trialIdx == 256
        fname = sprintf('Results/%s/Finalresults_across_M_and_ntrials_d80.mat', ResultFolder);
    else
        fname = sprintf('Results/%s/Finalresults_across_M_and_ntrials.mat', ResultFolder);
    end

    D = load(fname);
    GT_results      = D.GT_results;
    sampled_results = D.sampled_results;
    trial_vals      = D.ntrials_vals;
    trialidx        = find(trial_vals == trialIdx, 1);

    if isempty(trialidx)
        error('trialIdx=%d not found in ntrials_vals. Available: %s', trialIdx, mat2str(trial_vals));
    end

    PID_all = squeeze(sampled_results(:, trialidx, :, :, :, :));
    GT_all  = squeeze(GT_results(:,    trialidx, :, :, :, :));
    info_idx = info_amount;

    pick_bias = @(bias_name) local_pick_gauss_bias(PID_all, GT_all, bias_name, info_idx);

    series = struct('name',{},'color',{},'mean',{},'sem',{});
    thumbRule = floor(trialIdx/(4*3));

    for bi = 1:numel(bias_types)
        bname_raw = bias_types{bi};
        bname_norm = normalize_name(bname_raw);

        try
            [PID_v, GT] = pick_bias(bname_norm); % comps x cats x iters

            m_joint = (squeeze(mean(PID_v(1,:,:), 3, 'omitnan')) - squeeze(GT(1,:)));
            m_red   = (squeeze(mean(PID_v(5,:,:), 3, 'omitnan')) - squeeze(GT(5,:)));
            m_unq_t = (squeeze(mean(PID_v(2,:,:), 3, 'omitnan')) - squeeze(GT(2,:)));
            m_unq   = m_unq_t / 2;
            m_syn   = (squeeze(mean(PID_v(6,:,:), 3, 'omitnan')) - squeeze(GT(6,:)));

            s_joint = 3 * squeeze(std(PID_v(1,:,:), 0, 3, 'omitnan')) / sqrt(n_iter);
            s_red   = 3 * squeeze(std(PID_v(5,:,:), 0, 3, 'omitnan')) / sqrt(n_iter);
            s_unq   = 3 * squeeze(std(PID_v(2,:,:), 0, 3, 'omitnan')) / sqrt(n_iter) / 2;
            s_syn   = 3 * squeeze(std(PID_v(6,:,:), 0, 3, 'omitnan')) / sqrt(n_iter);

            series(end+1).name  = bname_raw; %#ok<AGROW>
            series(end).color   = get_method_color(bname_raw); % <<< ONLY COLOR SOURCE
            series(end).mean    = [m_joint(:), m_syn(:), m_red(:)];
            series(end).sem     = [s_joint(:), s_syn(:), s_red(:)];

        catch ME
            warning('Skipping bias "%s": %s', bname_raw, ME.message);
        end
    end

    if trialIdx == 256
        x_vals = D.M_vals;
    else
        x_vals = categories;
    end

    if isempty(insetRange)
        insetRange = [28 36];
    end

    x_label = 'Dimensions';
    cat_sufficient = 32;
    y_label = 'Bias [bits]';

else
    error('Union path not implemented in this function.');
end

if isempty(series)
    error('No bias series found.');
end

% ================= Plot: three adjacent tiles (Joint | Syn | Red) =======
for c = 1:3
    nexttile(fig_handle, subplot_idx + (c-1)); hold on;
    xline(thumbRule, ':', 'Color','k','LineWidth',lw, 'HandleVisibility','off');

    for k = 1:numel(series)
        m = series(k).mean(:,c);
        s = series(k).sem(:,c);
        col = series(k).color;

        xv = [x_vals(:); flipud(x_vals(:))];
        yv = [m + s; flipud(m - s)];

        % NOTE: kept identical; fill uses same color as line (as in your original approach)
        fill(xv, yv, col, 'FaceAlpha', 0.18, 'EdgeColor','none', 'HandleVisibility','off');
        plot(x_vals, m, 'LineWidth',lw, 'Color', col);
    end

    yline(0, 'k-', 'LineWidth', lw, 'HandleVisibility','off');
    xlim([min(x_vals) max(x_vals)]);

    if isGauss
        ylim([-2, 8]);
    else
        ylim([-.4, .2]);
    end

    if c~=1
        yticklabels({});
    end

    set(gca, 'FontSize', FontSize, 'LineWidth',lw);
    title(comp_names{c}, 'FontSize', FontSize, 'FontWeight', titleFontWeight);
    if c==1
        ylabel(y_label, 'FontSize', FontSize, 'FontWeight','normal');
    end
    xlabel(x_label, 'FontSize', FontSize);

    % ================= Legend text fixes (ONLY TEXT CHANGED) ==============
    if c==1
        if subplot_idx<=3
            % was: {'QE', 'shuff-sub', 'merged'}
            lgd = legend({'QE', 'shuffle', 'merged'}, ...
                'Box','off', 'FontSize', legendFontSize, 'Location','northwest');
        else
            % was: {'resampling', 'shuff-sub', 'resamp, shuff-sub','Venkatesh'}
            lgd = legend({'resampling', 'shuffle', 'merged', 'Venkatesh'}, ...
                'Box','off', 'FontSize', legendFontSize, 'Location','northwest');
        end
        lgd.ItemTokenSize = [12,6];
    end

    % ================= Inset (unchanged) =================================
    if addInset && c ~= 1
        parentAx = gca;
        parentPos = get(parentAx, 'Position');

        insetW = 0.40; insetH = 0.40;
        offX   = 0.1;  offY   = 0.58;

        left   = parentPos(1) + offX * parentPos(3);
        bottom = parentPos(2) + offY * parentPos(4);
        width  = insetW * parentPos(3);
        height = insetH * parentPos(4);

        axIn = axes('Position', [left, bottom, width, height]);
        hold(axIn,'on');

        xr = insetRange;
        idxInset = (x_vals >= xr(1)) & (x_vals <= xr(2));
        if ~any(idxInset)
            mid = round(numel(x_vals)/2);
            win = max(1, mid-2):min(numel(x_vals), mid+2);
            idxInset = false(size(x_vals)); idxInset(win) = true;
        end
        xv = x_vals(idxInset);

        for kk = 1:numel(series)
            m = series(kk).mean(idxInset, c);
            s = series(kk).sem(idxInset, c);
            col = series(kk).color;

            xfill = [xv(:); flipud(xv(:))];
            yfill = [m + s; flipud(m - s)];
            fill(xfill, yfill, col, 'FaceAlpha', 0.18, 'EdgeColor','none', 'HandleVisibility','off');
            plot(axIn, xv, m, 'LineWidth', lw, 'Color', col);
        end

        xline(thumbRule, ':', 'Color','k','LineWidth',lw, 'HandleVisibility','off');
        axis(axIn, 'tight');
        yl = ylim(axIn);
        pad = 0.06 * max(eps, diff(yl));

        if ~isGauss
            if c==2
                ylim(axIn, [-0.11, 0.05]);
            elseif c==3
                ylim(axIn, [-0.005, 0.012]);
            else
                ylim(axIn, [yl(1)-pad, yl(2)+pad]);
            end
        else
            if c==2
                ylim(axIn, [-1, 3.5]);
            elseif c==3
                ylim(axIn, [-.5, 3.5]);
            else
                ylim(axIn, [yl(1)-pad, yl(2)+pad]);
            end
        end

        if yl(1) < 0 && yl(2) > 0
            yline(axIn, 0, '-', 'Color', [0 0 0], 'LineWidth', lw);
        end

        xlim(axIn, [min(xv) max(xv)]);
        if ~isGauss
            xticks(axIn, 3:2:7);
        else
            x0 = min(xv); x1 = max(xv);
            t0 = ceil(x0/4)*4;
            t1 = floor(x1/4)*4;
            if t1 >= t0
                xticks(axIn, t0:4:t1);
            else
                xticks(axIn, []);
            end
        end

        set(axIn, 'Box','on', 'LineWidth',lw, 'FontSize', FontSize);
        hold(axIn,'off');
        axes(parentAx);
    end

    hold off;
end

% ================= Return "worst" metrics (unchanged) ====================
cat_idx = find(x_vals == cat_sufficient, 1);
if isempty(cat_idx), cat_idx = round(median(1:numel(x_vals))); end

if ~isGauss
    syn_gt_vec = GT_vec(:,2);
    red_gt_vec = GT_vec(:,3);
else
    if trialIdx == 256
        fname = sprintf('Results/%s/Finalresults_across_M_and_ntrials_d80.mat', ResultFolder);
    else
        fname = sprintf('Results/%s/Finalresults_across_M_and_ntrials.mat', ResultFolder);
    end
    D = load(fname);
    GT_results = D.GT_results;
    trial_vals = D.ntrials_vals;
    trialidx   = find(trial_vals == trialIdx, 1);
    GT_all2 = squeeze(GT_results(:, trialidx, :, :, :, :));
    [~, GT_base] = local_pick_gauss_bias(squeeze(D.sampled_results(:,trialidx,:,:,:,:)), GT_all2, 'plugin', info_amount);
    syn_gt_vec = squeeze(GT_base(6,:)).';
    red_gt_vec = squeeze(GT_base(5,:)).';
end

worst_all = struct();
worst_suff_all = struct();

for k = 1:numel(series)
    nm = matlab.lang.makeValidName(series(k).name);

    syn_end = series(k).mean(end,2);   syn_se_end = series(k).sem(end,2);
    red_end = series(k).mean(end,3);   red_se_end = series(k).sem(end,3);

    syn_s   = series(k).mean(cat_idx,2); syn_se_s = series(k).sem(cat_idx,2);
    red_s   = series(k).mean(cat_idx,3); red_se_s = series(k).sem(cat_idx,3);

    syn_gt_end = syn_gt_vec(end); red_gt_end = red_gt_vec(end);
    syn_gt_s   = syn_gt_vec(cat_idx); red_gt_s = red_gt_vec(cat_idx);

    worst_all.(nm).syndif       = syn_end;
    worst_all.(nm).syndif_sem   = syn_se_end;
    worst_all.(nm).synratio     = 100 * syn_end / syn_gt_end;
    worst_all.(nm).synratio_sem = 100 * syn_se_end / syn_gt_end;

    worst_all.(nm).reddif       = red_end;
    worst_all.(nm).reddif_sem   = red_se_end;
    worst_all.(nm).redratio     = 100 * red_end / red_gt_end;
    worst_all.(nm).redratio_sem = 100 * red_se_end / red_gt_end;

    worst_suff_all.(nm).syndif       = syn_s;
    worst_suff_all.(nm).syndif_sem   = syn_se_s;
    worst_suff_all.(nm).synratio     = 100 * syn_s / syn_gt_s;
    worst_suff_all.(nm).synratio_sem = 100 * syn_se_s / syn_gt_s;

    worst_suff_all.(nm).reddif       = red_s;
    worst_suff_all.(nm).reddif_sem   = red_se_s;
    worst_suff_all.(nm).redratio     = 100 * red_s / red_gt_s;
    worst_suff_all.(nm).redratio_sem = 100 * red_se_s / red_gt_s;
end

end % ======================= END MAIN FUNCTION ===========================


% ================================ SUBFUNCTIONS ===========================

function col = local_get_method_color(methodColors, rawName, normalize_name)
% Robust color lookup without changing plot behavior.
% Uses exact key if present; else tries normalized variants; else falls back to black.
fallback = [0 0 0];

if methodColors.isKey(rawName)
    col = methodColors(rawName);
    return;
end

nm = normalize_name(rawName);

% common aliases -> map to your keys (NO VISUAL CHANGE EXCEPT COLOR)
switch nm
    case {'shuff_sub','shuffsub','shuff'}
        nm2 = 'shuffle';
    case {'shuff_resamp','shuffresamp','shuff_resample','shuffresample','resamp_shuff','resamp_shuffle','shuff-resamp','shuff-resample'}
        nm2 = 'merged';
    case {'qe_shuff','qeshuff','qe_shuff_sub'}
        nm2 = 'merged';
    otherwise
        nm2 = nm;
end

% try direct normalized, then alias
if methodColors.isKey(nm2)
    col = methodColors(nm2);
elseif methodColors.isKey(nm)
    col = methodColors(nm);
else
    col = fallback;
end
end


function PID_v = local_load_discrete_slice(fname, bias_name, trialIdx)
% File layout (DISCRETE): trials x comps x iters x cats
% Returns K x C x I (cats x comps x iters)
bn = lower(bias_name);

switch bn
    case {'plugin','qe','shuff','shuffle','qeshuff','qe_shuff','merged'}

        % Map requested bias to variable name in file
        % (keeps your storage conventions; does not change plot results)
        switch bn
            case 'shuffle'
                vn = 'PID_v_shuff';
            case {'qeshuff','qe_shuff','merged'}
                % your files often use qeShuff camel-case
                vn = local_find_existing_var(fname, {'PID_v_qeShuff','PID_v_qeshuff','PID_v_qe_shuff'});
            otherwise
                vn = ['PID_v_' bn];
        end

        L = load(fname, vn);
        PID_v_all = L.(vn);                         % T x C x I x K
        PID_v = squeeze(PID_v_all(trialIdx,:,:,:)); % C x I x K
        PID_v = permute(PID_v, [3 1 2]);            % K x C x I

    case 'venkatesh'
        % Uses plugin + qeShuff (as in your original helper)
        Lp = load(fname, 'PID_v_plugin');
        vnGT = local_find_existing_var(fname, {'PID_v_qeShuff','PID_v_qeshuff','PID_v_qe_shuff'});
        Lm = load(fname, vnGT);

        P  = permute(squeeze(Lp.PID_v_plugin(trialIdx,:,:,:)), [3 1 2]);   % K x C x I
        GT = permute(squeeze(Lm.(vnGT)(end,:,:,:)), [3 1 2]);              % K x C x I

        joint_gt = squeeze(mean(GT(:,1,:), 3, 'omitnan'));
        joint_pl = squeeze(mean(P(:,1,:),  3, 'omitnan'));
        debias_factor = joint_gt ./ joint_pl;

        i1_gt = squeeze(mean(GT(:,3,:),3,'omitnan')) + squeeze(mean(GT(:,4,:),3,'omitnan'));
        i2_gt = squeeze(mean(GT(:,3,:),3,'omitnan')) + squeeze(mean(GT(:,5,:),3,'omitnan'));

        union_est = squeeze(sum(P(:,3:5,:), 2)); % K x I
        for it = 1:size(union_est,2)
            union_est(:,it) = union_est(:,it) .* debias_factor;
        end

        union_min = max(i1_gt, i2_gt);
        union_max = min(i1_gt + i2_gt, joint_gt);

        for it = 1:size(union_est,2)
            union_est(:,it) = max(union_est(:,it), union_min);
            union_est(:,it) = min(union_est(:,it), union_max);
        end

        syn = joint_gt - union_est;
        red = i1_gt + i2_gt - union_est;
        un1 = union_est - i2_gt;
        un2 = union_est - i1_gt;

        PID_v = zeros(size(P));
        PID_v(:,1,:) = repmat(joint_gt, 1, 1, size(P,3));
        PID_v(:,2,:) = syn;
        PID_v(:,3,:) = red;
        PID_v(:,4,:) = un1;
        PID_v(:,5,:) = un2;

    otherwise
        error('Unknown discrete bias type: %s', bias_name);
end
end


function vn = local_find_existing_var(fname, candidates)
% Returns the first candidate that exists in the MAT-file; errors if none exist.
vars = who('-file', fname);
for i = 1:numel(candidates)
    if any(strcmp(vars, candidates{i}))
        vn = candidates{i};
        return;
    end
end
error('None of the variables exist in %s. Tried: %s', fname, strjoin(candidates, ', '));
end


function [PID_v, GT] = local_pick_gauss_bias(PID_all, GT_all, bias_name, info_idx)
% GAUSS layout: comps x cats x bias x info x iters
switch lower(bias_name)
    case 'plugin',       b = 1;
    case {'resample','resampling'}, b = 2;
    case {'shuffle','shuff'},       b = 3;
    case {'shuff-resamp','shuffresamp','shuff_resamp','shuff-resample','shuff_resample'}
        PID2 = squeeze(PID_all(:, :, 2, info_idx, :));
        PID3 = squeeze(PID_all(:, :, 3, info_idx, :));
        PID_v = (PID2 + PID3) / 2;

        GT2  = squeeze(GT_all(:, :, 2, info_idx, :));
        GT3  = squeeze(GT_all(:, :, 3, info_idx, :));
        GT   = squeeze(mean((GT2 + GT3)/2, 3, 'omitnan'));
        return

    case 'venkatesh',    b = 4;
    otherwise
        error('Unknown Gaussian bias type: %s', bias_name);
end

PID_v = squeeze(PID_all(:, :, b, info_idx, :)); % comps x cats x iters
GT    = squeeze(GT_all(:,  :, b, info_idx));    % comps x cats x iters
end
