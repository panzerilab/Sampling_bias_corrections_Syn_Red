function PlotFigureS2gauss()
% Layout: 3 rows (trials 64,128,256) × 4 columns
% Output: Figures/S2gauss.svg

%% ---------- Figure & layout ----------
figure('Units','centimeters','Position',[1 1 21 12]);
tiledlayout(3,4,'TileSpacing','compact','Padding','compact');

%% ---------- Config ----------
trials = [64, 128, 256];
methods = {'plugin','shuffle','resample','merged','Venkatesh'};

method_labels = containers.Map( ...
    {'plugin','shuffle','resample','merged','Venkatesh'}, ...
    {'plugin','shuffle','resample','merged','Venkatesh'} );

% --- Colors matched to discrete figure ---
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

ResultFolder = 'Gauss';
info_amount = 4;
isGauss = true;
isUnion = false;

FontSize = 8;
legendFontSize = 8;

%% ---------- Plot ----------
for r = 1:numel(trials)
    trialIdx = trials(r);

    if r == 1
        categories = 4:2:10;
    elseif r == 2
        categories = 4:2:20;
    else
        categories = 4:4:40;
    end

    % Precompute curves
    curves = struct();
    for m = 1:numel(methods)
        curves.(methods{m}) = compute_gauss_means_for_method_weighted( ...
            trialIdx, categories, methods{m}, ...
            info_amount, ResultFolder, isGauss, isUnion);
    end

    comps  = {'joint','syn','red','unq'};
    titles = {'Joint','Synergy','Redundancy','Unique'};

    for c = 1:4
        nexttile; hold on;

        % ---- Y limits (unchanged) ----
        if r == 1
            ylims = (c==1)*[-1 5] + (c~=1)*[-1 2];
        elseif r == 2
            ylims = (c==1)*[-3 40]/2 + (c~=1)*[-3 20]/2;
        else
            ylims = (c==1)*[-5 100]/2 + (c~=1)*[-5 30]/2;
        end

        % ---- Plot methods ----
        for m = 1:numel(methods)
            meth = methods{m};
            clr  = methodColors(meth);

            y   = curves.(meth).(comps{c}).mean;
            sem = curves.(meth).(comps{c}).sem;

            fill([categories fliplr(categories)], ...
                 [y+sem fliplr(y-sem)], ...
                 clr,'FaceAlpha',0.18,'EdgeColor','none', ...
                 'HandleVisibility','off');

            plot(categories,y,'Color',clr,'LineWidth',2, ...
                 'DisplayName',method_labels(meth));
        end

        % ---- Gaussian rule-of-thumb: N > 12 d ----
        thumbRule = floor(trialIdx / 12);
        thumbRule = min(max(thumbRule, categories(1)), categories(end));
        xline(thumbRule,'--k','LineWidth',1.5,'HandleVisibility','off');

        % ---- Axes formatting ----
        yline(0,'k','LineWidth',1.2,'HandleVisibility','off');
        xlim([categories(1) categories(end)]);
        ylim(ylims);

        % ---- Titles ----
        if r == 1
            title(titles{c},'FontWeight','bold');
        end

        % ---- Y axis: columns 1 & 2 only ----
        if c <= 2
            ylabel('Residual bias [bits]');
        else
            yticklabels([]);
        end

        % ---- X axis: ALL rows ----
        xlabel('Number of dim. (d)');

        set(gca,'FontSize',FontSize,'LineWidth',1.2);

        % ---- Legend (once) ----
        if r == 1 && c == 1
            lgd = legend('Location','northwest');
            lgd.Box = 'off';
            lgd.FontSize = legendFontSize;
        end

        hold off;
    end
end

%% ---------- Save ----------
set(gcf, 'Renderer', 'painters');
outDir = 'Figures_mat';
if ~exist(outDir, 'dir')
    mkdir(outDir);
end
exportgraphics(gcf, fullfile(outDir, 'Figure_S2_gauss.svg'), 'ContentType','vector');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = compute_gauss_means_for_method_weighted( ...
    trialIdx, categories, bias_correction, ...
    info_amount, ResultFolder, isGauss, isUnion)

% Gaussian non-Union helper

set(0,'DefaultTextFontName','Arial');
set(0,'DefaultAxesFontName','Arial');
set(0,'DefaultLegendFontName','Arial');

out = struct();
out.categories = categories(:)';

assert(isGauss && ~isUnion, ...
    'This helper is configured for Gaussian non-Union.');

% Load Gaussian results
if trialIdx == 256
    filename = sprintf('Results/%s/Finalresults_across_M_and_ntrials_d80.mat', ResultFolder);
else
    filename = sprintf('Results/%s/Finalresults_across_M_and_ntrials.mat', ResultFolder);
end
D = load(filename);

GT_all   = D.GT_results;
PID_all  = D.sampled_results;
trials_v = D.ntrials_vals;

trialidx = find(trials_v == trialIdx);
PID = squeeze(PID_all(:, trialidx, :, :, :, :));
GT  = squeeze(GT_all(:,  trialidx, :, :, :, :));
info_idx = info_amount;

get_slice = @(b) deal( ...
    squeeze(PID(:,1:numel(categories),b,info_idx,:)), ...
    squeeze(GT(:, 1:numel(categories),b,info_idx)) );

switch lower(bias_correction)
    case 'plugin'
        [PID_v, GT_v] = get_slice(1);
    case 'resample'
        [PID_v, GT_v] = get_slice(2);
    case 'shuffle'
        [PID_v, GT_v] = get_slice(3);
    case 'venkatesh'
        [PID_v, GT_v] = get_slice(4);
    case 'merged'
        [PID_res, ~]  = get_slice(2);
        [PID_shf, ~]  = get_slice(3);
        [~, GT_v]     = get_slice(1);
        PID_v = 0.5*PID_res + 0.5*PID_shf;
    otherwise
        error('Unknown method');
end

n_iter = 100;

mean_Joint = mean(PID_v(1,:,:),3) - GT_v(1,:);
mean_RED   = mean(PID_v(5,:,:),3) - GT_v(5,:);
mean_UNQ   = (mean(PID_v(2,:,:),3) - GT_v(2,:)) / 2;
mean_SYN   = mean(PID_v(6,:,:),3) - GT_v(6,:);

SEM_Joint = 3*std(PID_v(1,:,:),[],3)/sqrt(n_iter);
SEM_RED   = 3*std(PID_v(5,:,:),[],3)/sqrt(n_iter);
SEM_UNQ   = (3*std(PID_v(2,:,:),[],3)/sqrt(n_iter))/2;
SEM_SYN   = 3*std(PID_v(6,:,:),[],3)/sqrt(n_iter);

out.joint.mean = mean_Joint(:)'; out.joint.sem = SEM_Joint(:)';
out.red.mean   = mean_RED(:)';   out.red.sem   = SEM_RED(:)';
out.unq.mean   = mean_UNQ(:)';   out.unq.sem   = SEM_UNQ(:)';
out.syn.mean   = mean_SYN(:)';   out.syn.sem   = SEM_SYN(:)';
end
