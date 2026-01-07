function PlotFigureS3(trial_categories, category, atom, simul_case, ...
    bias_corrections, info_amount, redundancy_measure, ...
    subplot_idx, fig_handle, ResultFolder, max_v, min_v)

set(0,'DefaultTextFontName','Arial');
set(0,'DefaultAxesFontName','Arial');

simul_case  = char(simul_case);
info_amount = char(info_amount);
redundancy_measure = char(redundancy_measure);

n_iter = 100;
folder = ResultFolder;
numColors = 256;

greenColors = [229 245 224;199 233 192;161 217 155;116 196 118;65 171 93;35 139 69]/255;
SynColormap = interp1(linspace(0,1,size(greenColors,1)),greenColors,linspace(0,1,numColors));

yellowColors = [254 227 145;254 196 79;254 153 41;217 95 14;153 52 4]/255;
JointColormap = interp1(linspace(0,1,size(yellowColors,1)),yellowColors,linspace(0,1,numColors));

greyColors = [217 217 217;189 189 189;150 150 150;115 115 115;82 82 82]/255;
UniqueColormap = interp1(linspace(0,1,size(greyColors,1)),greyColors,linspace(0,1,numColors));

blueColors = [198 219 239;158 202 225;107 174 214;66 146 198;33 113 181;8 69 148]/255;
RedColormap = interp1(linspace(0,1,size(blueColors,1)),blueColors,linspace(0,1,numColors));

for biasIdx = 1:numel(bias_corrections)

    bias_correction = bias_corrections{biasIdx};

    filename = sprintf('Results/%s/Simuldata_%s_%s_%s.mat', ...
        folder, simul_case, info_amount, redundancy_measure);

    nameLoad = sprintf('PID_v_%s', bias_correction);
    S = load(filename, nameLoad);
    PID_v = S.(nameLoad);

    categories   = trial_categories;
    isTrialSweep = true;
    thumbRule    = 4^4;

    switch atom
        case 'Joint'
            mean_v = squeeze(mean(PID_v(:,1,:,:),3,'omitnan'));
            SEM = squeeze(3*std(PID_v(:,1,:,:),1,3,'omitnan')/sqrt(n_iter));
            cmap = JointColormap;
            hasTitle = true;

        case 'Synergy'
            mean_v = squeeze(mean(PID_v(:,2,:,:),3,'omitnan'));
            SEM = squeeze(3*std(PID_v(:,2,:,:),1,3,'omitnan')/sqrt(n_iter));
            cmap = SynColormap;
            hasTitle = false;

        case 'Redundancy'
            mean_v = squeeze(mean(PID_v(:,3,:,:),3,'omitnan'));
            SEM = squeeze(3*std(PID_v(:,3,:,:),1,3,'omitnan')/sqrt(n_iter));
            cmap = RedColormap;
            hasTitle = false;

        case 'Unique'
            UNQ = mean(PID_v(:,4:5,:,:),2);
            mean_v = squeeze(mean(UNQ,3,'omitnan'));
            SEM = squeeze(3*std(UNQ,1,3,'omitnan')/sqrt(n_iter));
            cmap = UniqueColormap;
            hasTitle = false;
    end

    if length(categories) ~= size(mean_v,1)
        error('Dimension mismatch')
    end

    nLines = size(mean_v,2);
    linecolor = cmap(round(linspace(1,numColors,nLines)),:);
    rowtxttrue = (biasIdx == 1);

    plotSubplot(mean_v, SEM, categories, linecolor, ...
        fig_handle, subplot_idx, bias_correction, ...
        rowtxttrue, atom, min_v, max_v, hasTitle, ...
        isTrialSweep, thumbRule);

    subplot_idx = subplot_idx + 1;
end
end


function plotSubplot(mean_v, SEM, categories, color, fig_handle, ...
    subplot_idx, titletxt, rowtxttrue, rowtxt, min_val, max_val, ...
    hasTitle, isTrialSweep, thumbRule)

nexttile(fig_handle, subplot_idx);
hold on;
set(gca,'FontSize',18,'LineWidth',1.2);

switch titletxt
    case 'shuff',     titletxt = 'shuffle';
    case 'qeShuff',   titletxt = 'merged';
    case 'infoCorr',  titletxt = 'resample';
    case 'qe',        titletxt = 'QE';
end

x = categories * 4;

for i = 1:size(mean_v,2)
    plot(x, mean_v(:,i), 'Color', color(i,:), 'LineWidth', 1.2);
end

xline(thumbRule, ':k', 'LineWidth', 2, 'HandleVisibility','off');

if isTrialSweep
    set(gca,'XScale','log');
end

xlim([min(x) max(x)]);
ylim([min_val max_val]);

if hasTitle
    title(sprintf('\\bf%s\\rm', titletxt));
end

nCols = 5;
rowIdx = ceil(subplot_idx / nCols);
colIdx = mod(subplot_idx-1,nCols)+1;

if rowIdx == 3
    xticks([1e2 1e3 1e4]);
    xticklabels({'10^2','10^3','10^4'});
    xlabel('Trials (N)');
else
    xticklabels({});
end

if colIdx ~= 1
    yticklabels({});
end

if rowtxttrue && colIdx == 1
    ylabel(sprintf('\\bf%s\\rm\nInformation [bits]', rowtxt), ...
        'FontWeight','normal');
end

hold off;
end
