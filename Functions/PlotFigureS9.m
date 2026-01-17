function PlotFigureS9(trial_categories, category, atom, bias_corrections, ...
    subplot_idx, fig_handle, ResultFolder, max_v, min_v)

% ==================== GLOBAL STYLE ====================
set(0,'DefaultTextFontName','Arial');
set(0,'DefaultAxesFontName','Arial');

n_iter = 100;
numColors = 256;

% ==================== COLORMAPS ====================
greenColors = [229 245 224;199 233 192;161 217 155;116 196 118;65 171 93;35 139 69]/255;
SynColormap = interp1(linspace(0,1,size(greenColors,1)),greenColors,linspace(0,1,numColors));

yellowColors = [254 227 145;254 196 79;254 153 41;217 95 14;153 52 4]/255;
JointColormap = interp1(linspace(0,1,size(yellowColors,1)),yellowColors,linspace(0,1,numColors));

blueColors = [198 219 239;158 202 225;107 174 214;66 146 198;33 113 181;8 69 148]/255;
RedColormap = interp1(linspace(0,1,size(blueColors,1)),blueColors,linspace(0,1,numColors));

% ==================== LOAD DATA ====================
data   = load(sprintf('Results/%s/Finalresults_across_M_and_ntrials.mat', ResultFolder));
PID_all = data.sampled_results;
M_vals  = data.M_vals;

dimidx  = find(M_vals == category,1);
PID_all = squeeze(PID_all(:,:,dimidx,:,:,:));

thumbRule = 4 *3* category;

% ==================== MAIN LOOP ====================
for biasIdx = 1:numel(bias_corrections)

    bias = bias_corrections{biasIdx};

    switch bias
        case 'plugin'
            PID_v = squeeze(PID_all(:,:,1,:,:));
        case 'resample'
            PID_v = squeeze(PID_all(:,:,2,:,:));
        case {'shuff','shuffle','shuff-sub'}
            PID_v = squeeze(PID_all(:,:,3,:,:));
        case 'Venkatesh'
            PID_v = squeeze(PID_all(:,:,4,:,:));
        case 'shuff-resample'
            PID_v = squeeze(mean(PID_all(:,:,2:3,:,:),3));
        otherwise
            error('PlotFigureS9:UnknownBias','Unknown bias correction "%s".',bias);
    end

    % ==================== ATOM SELECTION ====================
    switch atom
        case 'Joint'
            mean_v = squeeze(mean(PID_v(1,:,:,:),4,'omitnan'));
            SEM    = squeeze(3*std(PID_v(1,:,:,:),1,4,'omitnan')/sqrt(n_iter));
            cmap   = JointColormap;
            hasTitle = true;

        case 'Synergy'
            mean_v = squeeze(mean(PID_v(6,:,:,:),4,'omitnan'));
            SEM    = squeeze(3*std(PID_v(6,:,:,:),1,4,'omitnan')/sqrt(n_iter));
            cmap   = SynColormap;
            hasTitle = false;

        case 'Redundancy'
            mean_v = squeeze(mean(PID_v(5,:,:,:),4,'omitnan'));
            SEM    = squeeze(3*std(PID_v(5,:,:,:),1,4,'omitnan')/sqrt(n_iter));
            cmap   = RedColormap;
            hasTitle = false;

        otherwise
            error('PlotFigureS9:UnknownAtom', ...
                  'Unknown atom "%s". Use Joint, Synergy, or Redundancy.', atom);
    end

    % ==================== COLOR ASSIGNMENT ====================
    nLines = size(mean_v,2);
    linecolor = flip(cmap(round(linspace(1,numColors,nLines)),:),1);

    min_plot = min_v;
    max_plot = max_v;

    if biasIdx > 1
        switch atom
            case 'Joint'
                max_plot = 10;
            otherwise
                max_plot = 4;
        end
    end

    % ==================== PLOT ====================
    plotSubplotS9(mean_v, SEM, trial_categories, linecolor, ...
        fig_handle, subplot_idx, bias, atom, ...
        min_plot, max_plot, hasTitle, thumbRule);

    subplot_idx = subplot_idx + 1;
end
end


% =======================================================================
function plotSubplotS9(mean_v, SEM, x, color, fig_handle, subplot_idx, ...
    titletxt, rowtxt, min_val, max_val, hasTitle, thumbRule)

nexttile(fig_handle, subplot_idx);
hold on;

set(gca, ...
    'FontSize',14, ...
    'LineWidth',1.2, ...
    'LooseInset',max(get(gca,'TightInset'),0.02));

switch titletxt
    case {'shuff','shuffle','shuff-sub'}
        titletxt = 'shuffle';
    case 'shuff-resample'
        titletxt = 'merged';
end

for i = 1:size(mean_v,2)
    plot(x, mean_v(:,i), 'Color', color(i,:), 'LineWidth', 1.2);
end

xline(thumbRule, '--', 'LineWidth', 2, 'HandleVisibility','off');

set(gca,'XScale','log');
xlim([min(x) max(x)]);
ylim([min_val max_val]);

if hasTitle
    title(sprintf('\\bf%s\\rm', titletxt));
end

nCols = 5;
colIdx = mod(subplot_idx-1,nCols)+1;
rowIdx = ceil(subplot_idx/nCols);

if rowIdx == 3
    xticks([1e2 1e3 1e4]);
    xticklabels({'10^2','10^3','10^4'});
    xlabel('Trials (N)');
else
    xticklabels({});
end

if colIdx <= 2
    ylabel(sprintf('\\bf%s\\rm\nInformation [bits]', rowtxt), ...
        'FontWeight','normal');
else
    yticklabels({});
end

hold off;
end
