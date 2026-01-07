function PlotFigureS2discr()
% MAKE_BIN_SWEEP_FIGURE_COLUMNS_BY_INFO
% Layout: 3 rows (trials) × 4 columns (information types).
% Output: Figures/S2discr.svg

%% Figure & layout
figure_handle = figure('Units', 'centimeters', 'Position', [1, 1, 21, 12]);
tiled = tiledlayout(figure_handle, 3, 4, ...
    'TileSpacing', 'compact', 'Padding', 'compact');

%% Config
trials = [1, 2, 3];                                % rows
info_cols = {'Joint','Syn','Red','Unq'};           % columns
methods = {'plugin','shuff','qe','qeShuff'};       % methods

method_labels = containers.Map( ...
    {'plugin','shuff','qe','qeShuff'}, ...
    {'plugin','shuffle','QE','merged'} );

% Colors matched to attached figure
methodColors = containers.Map( ...
    {'plugin','qe','resample','resampling','shuff','shuffle','qe_shuff','qeShuff', ...
     'merged','shuff-resamp','weighted','venkatesh','Venkatesh'}, ...
    { [0.60 0.60 0.60], ...         % plugin (gray)
      [0.00 0.447 0.741], ...       % QE (blue)
      [0.929 0.694 0.125], ...      % resample (YELLOW)
      [0.929 0.694 0.125], ...      % resampling (YELLOW)
      [0.85  0.325 0.098], ...      % shuff (orange)
      [0.85  0.325 0.098], ...      % shuffle (orange)
      [0.494 0.184 0.556], ...      % qe_shuff (purple)
      [0.494 0.184 0.556], ...      % qeShuff (purple)
      [0.494 0.184 0.556], ...      % merged (purple)
      [0.494 0.184 0.556], ...      % shuff-resamp (purple)
      [0.60  0.60  0.60], ...       % weighted (gray)
      [0.466 0.674 0.188], ...      % venkatesh (green)
      [0.466 0.674 0.188] } );      % Venkatesh (green)

% Data params
ResultFolder = 'Bin_sweep';
categories = 2:9;
simul_case = 'bit_of_all';
info_amount = 'low';
redundancy_measure = 'I_Broja';
isGauss = false;
isUnion = false;

% Fonts
FontSize = 8;
legendFontSize = 8;

%% Plot grid
subplot_idx = 1;

for r = 1:numel(trials)
    trialIdx = trials(r);

    % Per-row y-limits
    if r == 1
        y_min = -0.55; y_max = 1.0;
    elseif r == 2
        y_min = -0.40; y_max = 0.6;
    else
        y_min = -0.22; y_max = 0.4;
    end

    % Precompute curves
    curves = struct();
    for m = 1:numel(methods)
        meth = methods{m};
        curves.(meth) = compute_means_for_method( ...
            trialIdx, categories, simul_case, meth, ...
            info_amount, redundancy_measure, ...
            ResultFolder, isGauss, isUnion);
    end

    % Columns
    for c = 1:numel(info_cols)
        nexttile(tiled, subplot_idx); hold on;

        switch c
            case 1, comp = 'joint'; titleTxt = 'Joint';
            case 2, comp = 'syn';   titleTxt = 'Synergy';
            case 3, comp = 'red';   titleTxt = 'Redundancy';
            case 4, comp = 'unq';   titleTxt = 'Unique';
        end

        % Plot methods
        for m = 1:numel(methods)
            meth = methods{m};
            clr  = methodColors(meth);
            cats = curves.(meth).categories;

            y   = curves.(meth).(comp).mean;
            sem = curves.(meth).(comp).sem;

            % SEM band
            fill([cats(:); flipud(cats(:))], ...
                 [y(:)+sem(:); flipud(y(:)-sem(:))], ...
                 clr, 'FaceAlpha', 0.18, ...
                 'EdgeColor', 'none', ...
                 'HandleVisibility', 'off');

            % Mean line
            plot(cats, y, 'Color', clr, ...
                'LineWidth', 2.0, ...
                'DisplayName', method_labels(meth));
        end

        % Thumb-rule vertical line (DASHED)
        N = (trialIdx==1)*256 + (trialIdx==2)*512 + (trialIdx==3)*1024;
        thumbRule = floor(nthroot(N/4, 3));
        xline(thumbRule, '--', 'Color', 'k', ...
              'LineWidth', 1.5, 'HandleVisibility', 'off');

        % Axes
        xlim([categories(1) categories(end)]);
        xticks([2 4 6 8]);
        ylim([y_min y_max]);
        yline(0, 'k-', 'LineWidth', 1.2, ...
              'HandleVisibility', 'off');

        set(gca, 'FontSize', FontSize, 'LineWidth', 1.2);
        xlabel('Number of bins');

        if c == 1
            ylabel('Residual bias [bits]');
        end

        if r == 1
            title(titleTxt, 'FontWeight', 'bold');
        end

        % Legend inside plot (top-left panel)
        if r == 1 && c == 1
            lgd = legend('Location','northwest');
            lgd.Box = 'off';
            lgd.FontSize = legendFontSize;
        end

        % Cleanup shared axes
        if r ~= numel(trials)
            xticklabels([]); xlabel([]);
        end
        if c ~= 1
            yticklabels([]); ylabel([]);
        end

        hold off;
        subplot_idx = subplot_idx + 1;
    end
end

%% Save
set(gcf, 'Renderer', 'painters');
outDir = 'Figures_mat';
if ~exist(outDir, 'dir')
    mkdir(outDir);
end
exportgraphics(gcf, fullfile(outDir, 'Figure_S2_discrete.svg'), 'ContentType','vector');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = compute_means_for_method(trialIdx, categories, simul_case, bias_correction, info_amount, redundancy_measure, ResultFolder, isGauss, isUnion)
% Returns struct with fields:
% out.categories
% out.joint.mean/sem, out.syn.mean/sem, out.red.mean/sem, out.unq.mean/sem
%
% This mirrors your discrete (non-Union) branch calculations and keeps
% your ground-truth handling through the qeShuff final slice.

% Fonts (match your defaults; harmless to set here too)
set(0, 'DefaultTextFontName', 'Arial');
set(0, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultUicontrolFontName', 'Arial');
set(0, 'DefaultUitableFontName', 'Arial');
set(0, 'DefaultUipanelFontName', 'Arial');
set(0, 'DefaultLegendFontName', 'Arial');

out = struct();
out.categories = categories(:)';

assert(~isUnion, 'This helper is configured for non-Union discrete case.');

filename = sprintf('Results/%s/Simuldata_%s_%s_%s.mat', ResultFolder, simul_case, info_amount, redundancy_measure);
n_iter = 100;

% Load method and ground-truth (taken from qeShuff end, as in your code)
if strcmpi(bias_correction,'venkatesh')
    % Kept for compatibility; not used in main flow
    nameLoad =  sprintf('PID_v_%s', 'plugin');
    PID_v_big = load(filename, nameLoad);
    PID_v_plug = PID_v_big.(nameLoad);

    PID_v_merg = load(filename, sprintf('PID_v_%s', 'qeShuff'));
    PID_v_merg=PID_v_merg.(sprintf('PID_v_%s', 'qeShuff'));
    GroundTruth = squeeze(PID_v_merg(end,:, :, :));
    PID_v_plug = squeeze(PID_v_plug(trialIdx,:, :, :));
    PID_v = PID_v_plug;

    debias_factor =  squeeze(mean(GroundTruth(1,:,:), 2, 'omitnan'))./squeeze(mean(PID_v_plug(1,:,:), 2, 'omitnan'));

    joint_gt = squeeze(mean(GroundTruth(1,:,:), 2, 'omitnan'));
    i1_gt    = squeeze(mean(GroundTruth(3,:,:), 2, 'omitnan'))+squeeze(mean(GroundTruth(4,:,:), 2, 'omitnan'));
    i2_gt    = squeeze(mean(GroundTruth(3,:,:), 2, 'omitnan'))+squeeze(mean(GroundTruth(5,:,:), 2, 'omitnan'));
    union = squeeze(sum(PID_v_plug(3:5,:,:),1))' .* debias_factor;
    union = max(union, max(i1_gt, i2_gt));
    union = min(union, min(i1_gt + i2_gt, joint_gt));
    syn = joint_gt-union;
    red = i1_gt + i2_gt - union;
    un1 = union - i2_gt;
    un2 = union - i1_gt;
    PID_v(1,:,:) = repmat(joint_gt', 100, 1);
    PID_v(2,:,:) = syn';
    PID_v(3,:,:) = red';
    PID_v(4,:,:) = un1';
    PID_v(5,:,:) = un2';

elseif strcmpi(bias_correction,'diff')
    % Kept for completeness; not used in final layout
    nameLoad =  sprintf('PID_v_%s', 'plugin');
    PID_v_big = load(filename, nameLoad);
    PID_v_plug = PID_v_big.(nameLoad);

    PID_v_merg = load(filename, sprintf('PID_v_%s', 'qeShuff'));
    PID_v_merg=PID_v_merg.(sprintf('PID_v_%s', 'qeShuff'));
    GroundTruth = squeeze(PID_v_merg(end,:, :, :));

    nameLoad =  sprintf('PID_v_%s', 'qe');
    PID_v_qe = load(filename, nameLoad);
    PID_v_qe = squeeze(PID_v_qe.(nameLoad));
    PID_v_qe = squeeze(PID_v_qe(trialIdx,:, :, :));

    nameLoad =  sprintf('PID_v_%s', 'shuff');
    PID_v_sh = load(filename, nameLoad);
    PID_v_sh = squeeze(PID_v_sh.(nameLoad));
    PID_v_sh = squeeze(PID_v_sh(trialIdx,:, :, :));

    PID_v = abs(PID_v_qe-PID_v_sh);

else
    nameLoad =  sprintf('PID_v_%s', bias_correction);
    PID_v_big = load(filename, nameLoad);
    PID_v = PID_v_big.(nameLoad);

    PID_v_merg = load(filename, sprintf('PID_v_%s', 'qeShuff'));
    PID_v_merg=PID_v_merg.(sprintf('PID_v_%s', 'qeShuff'));
    GroundTruth = squeeze(PID_v_merg(end,:, :, :));
    PID_v = squeeze(PID_v(trialIdx,:, :, :));
end

% Compute means vs. GT (discrete, non-Union branch from your function)
GroundTruth_UNQ = mean(GroundTruth(4:5,:,:),1);
mean_Joint_GT = squeeze(mean(GroundTruth(1,:,:), 2, 'omitnan'));
mean_RED_GT   = squeeze(mean(GroundTruth(3,:,:), 2, 'omitnan'));
mean_UNQ_GT   = squeeze(mean(GroundTruth_UNQ,     2, 'omitnan'));
mean_SYN_GT   = squeeze(mean(GroundTruth(2,:,:),  2, 'omitnan'));

UNQ = mean(PID_v(4:5,:,:),1);
mean_Joint = squeeze(mean(PID_v(1,:,:), 2, 'omitnan')) - mean_Joint_GT;
mean_RED   = squeeze(mean(PID_v(3,:,:), 2, 'omitnan')) - mean_RED_GT;
mean_UNQ   = squeeze(mean(UNQ,          2, 'omitnan')) - mean_UNQ_GT;
mean_SYN   = squeeze(mean(PID_v(2,:,:), 2, 'omitnan')) - mean_SYN_GT;

SEM_Joint = squeeze(3 * std(PID_v(1,:,:),1,2,'omitnan') / sqrt(n_iter));
SEM_RED   = squeeze(3 * std(PID_v(3,:,:),1,2,'omitnan') / sqrt(n_iter));
SEM_UNQ   = squeeze(3 * std(UNQ,        1,2,'omitnan') / sqrt(n_iter));
SEM_SYN   = squeeze(3 * std(PID_v(2,:,:),1,2,'omitnan') / sqrt(n_iter));

% Pack results
out.joint.mean = mean_Joint(:)'; out.joint.sem = SEM_Joint(:)';
out.red.mean   = mean_RED(:)';   out.red.sem   = SEM_RED(:)';
out.unq.mean   = mean_UNQ(:)';   out.unq.sem   = SEM_UNQ(:)';
out.syn.mean   = mean_SYN(:)';   out.syn.sem   = SEM_SYN(:)';


end
