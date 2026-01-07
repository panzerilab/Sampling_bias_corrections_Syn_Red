function PlotFigureS8(trial_categories, category, atom, simul_case, bias_corrections, info_amount, redundancy_measure, subplot_idx, fig_handle, ResultFolder,max_v, min_v, isGauss)
set(0, 'DefaultTextFontName', 'Arial');
set(0, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultUicontrolFontName', 'Arial');
set(0, 'DefaultUitableFontName', 'Arial');
set(0, 'DefaultUipanelFontName', 'Arial');
set(0, 'DefaultLegendFontName', 'Arial');


% Convert inputs to strings if necessary
simul_case = char(simul_case);
info_amount = char(info_amount);
n_iter = 100;
folder = ResultFolder;

% Custom colors for each component of the PID
customColors_info_PID = [
    0.9290, 0.6940, 0.1250;  % Color for 'Joint'
    0, 0.4470, 0.7410;       % Color for 'Red'
    0.8500, 0.3250, 0.0980;  % Color for 'Union'
    0.4940, 0.1840, 0.5560;  % Color for 'SR'
    0.4660, 0.6740, 0.1880;  % Color for 'Syn'
    0.5, 0.5, 0.5;           % Color for 'U1'
    0.4, 0.26, 0.13;         % Color for 'U2'
    0.4940, 0.1840, 0.5560]; % Color for 'PInd'

greenColors = [ ...
    229, 245, 224;
    199, 233, 192;
    161, 217, 155;
    116, 196, 118;
    65, 171, 93;
    35, 139, 69
    ] / 255;

numColors = 256;
SynColormap = interp1(linspace(0, 1, size(greenColors, 1)), greenColors, linspace(0, 1, numColors));

yellowColors = [ ...
    254, 227, 145;
    254, 196, 79;
    254, 153, 41;
    217, 95, 14;
    153, 52, 4
    ] / 255;

numColors = 256;
JointColormap = interp1(linspace(0, 1, size(yellowColors, 1)), yellowColors, linspace(0, 1, numColors));

greyColors = [ ...
    217, 217, 217;
    189, 189, 189;
    150, 150, 150;
    115, 115, 115;
    82, 82, 82
    ] / 255;

numColors = 256;
UniqueColormap = interp1(linspace(0, 1, size(greyColors, 1)), greyColors, linspace(0, 1, numColors));

blueColors = [ ...
    198, 219, 239;
    158, 202, 225;
    107, 174, 214;
    66, 146, 198;
    33, 113, 181;
    8, 69, 148
    ] / 255;

numColors = 256;
RedColormap = interp1(linspace(0, 1, size(blueColors, 1)), blueColors, linspace(0, 1, numColors));




for biasIdx = 1:length(bias_corrections)
    bias_correction = bias_corrections{biasIdx};

    % Set up the figure and subplot
    nexttile(fig_handle, subplot_idx);
    %subplot_handle = subplot(rows, cols, subplot_idx);
    hold on;


    if length(trial_categories) > 1
        filename = sprintf('Results/%s/Finalresults_across_M_and_ntrials.mat', ResultFolder);
    else
        filename = sprintf('Results/%s/Finalresults_across_M_and_ntrials_d80.mat', ResultFolder);
    end
    data = load(filename);
    GroundTruth =  data.GT_results;
    PID_v =  data.sampled_results;
    M_vals =  data.M_vals;

    if length(trial_categories) > 1
        dim = category;
        dimidx = find(M_vals==dim);
        thumbRule = 4*dim;
        PID_v = squeeze(PID_v(:,:,dimidx,:,:,:));
        GroundTruth = squeeze(GroundTruth(:,:,dimidx,:,:));
        if strcmp(bias_correction, 'plugin')
            PID_v = squeeze(PID_v(:, :, 1, :, :));
        elseif strcmp(bias_correction, 'shuffle')
            PID_v = squeeze(PID_v(:, :, 3, :, :));
        elseif strcmp(bias_correction, 'resample')
            PID_v = squeeze(PID_v(:, :, 2, :, :));
        elseif strcmp(bias_correction, 'Venkatesh')
            PID_v = squeeze(PID_v(:, :, 4, :, :));
        elseif strcmp(bias_correction, 'shuff-resample')
            PID_v = squeeze(mean(PID_v(:, :, 2:3, :, :),3));
        end
        if strcmp(atom, 'Joint')
            mean_v = squeeze(mean(mean(PID_v(1, :,  :, :),4, 'omitnan'),5, 'omitnan'));
            SEM  = squeeze(3 * (std(PID_v(1, :,  :, :),1,4,'omitnan') / sqrt(n_iter)));
            [nX, nLines] = size(mean_v);
            colorIndices = round(linspace(1, numColors, nLines));
            colorIndices = flip(colorIndices);
            linecolor = JointColormap(colorIndices, :);
            hasXlabel = false;
            hasTitle = true;
        elseif strcmp(atom, 'Synergy')
            mean_v = squeeze(mean(mean(PID_v(6, :,  :, :),4, 'omitnan'),5, 'omitnan'));
            SEM  = squeeze(3 * (std(PID_v(6, :,  :, :),1,4,'omitnan') / sqrt(n_iter)));
            [nX, nLines] = size(mean_v);
            colorIndices = round(linspace(1, numColors, nLines));
            colorIndices = flip(colorIndices);
            linecolor = SynColormap(colorIndices, :);
            hasXlabel = false;
            hasTitle = false;
        elseif strcmp(atom, 'Redundancy')
            mean_v = squeeze(mean(mean(PID_v(5, :,  :, :),4, 'omitnan'),5, 'omitnan'));
            SEM  = squeeze(3 * (std(PID_v(5, :,  :, :),1,4,'omitnan') / sqrt(n_iter)));
            [nX, nLines] = size(mean_v);
            colorIndices = round(linspace(1, numColors, nLines));
            colorIndices = flip(colorIndices);
            linecolor = RedColormap(colorIndices, :);
            hasXlabel = true;
            hasTitle = false;
        elseif strcmp(atom, 'Unique')
            PID_v(2, :, :, :) = PID_v(2, :,  :, :)/2;
            mean_v = squeeze(mean(mean(PID_v(2, :,  :, :),4, 'omitnan'),5, 'omitnan'));
            SEM  = squeeze(3 * (std(PID_v(2, :,  :, :),1,4,'omitnan') / sqrt(n_iter)));
            [nX, nLines] = size(mean_v);
            colorIndices = round(linspace(1, numColors, nLines));
            colorIndices = flip(colorIndices);
            linecolor = UniqueColormap(colorIndices, :);
            hasXlabel = true;
            hasTitle = false;
        end
        categories = trial_categories;
        isTrialSweep = true;
    else
        trial_vals =  data.ntrials_vals;
        trialidx = find(trial_vals==trial_categories);
        thumbRule = floor(trial_categories/4);
        PID_v = squeeze(PID_v(:,trialidx,:,:,:,:));
        GroundTruth = squeeze(GroundTruth(:,trialidx,:,:,:));

        if strcmp(bias_correction, 'plugin')
            PID_v = squeeze(PID_v(:, 1:length(category), 1, :, :));
            GroundTruth = squeeze(GroundTruth(:, 1:length(category), 1, :));
        elseif strcmp(bias_correction, 'shuffle')
            PID_v = squeeze(PID_v(:, 1:length(category), 3,  :, :));
            GroundTruth = squeeze(GroundTruth(:, 1:length(category), 3, :));
        elseif strcmp(bias_correction, 'resample')
            PID_v = squeeze(PID_v(:, 1:length(category), 2, :, :));
            GroundTruth = squeeze(GroundTruth(:, 1:length(category), 2, :));
        elseif strcmp(bias_correction, 'Venkatesh')
            PID_v = squeeze(PID_v(:, 1:length(category), 4, :, :));
            GroundTruth = squeeze(GroundTruth(:, 1:length(category), 4, :));
        elseif strcmp(bias_correction, 'shuff-resamp')
            PID_v = squeeze(mean(PID_v(:, 1:length(category), 2:3, :, :),3));
            GroundTruth = squeeze(mean(GroundTruth(:, 1:length(category), 2:3, :),3));
        end
        if strcmp(atom, 'Joint')
            mean_v = squeeze(mean(PID_v(1, :,  :, :),4, 'omitnan')) - squeeze(GroundTruth(1,:,:));
            SEM  = squeeze(3 * (std(PID_v(1, :,  :, :),1,4,'omitnan') / sqrt(n_iter)));
            [nX, nLines] = size(mean_v);
            colorIndices = round(linspace(1, numColors, nLines));
            colorIndices = flip(colorIndices);
            linecolor = JointColormap(colorIndices, :);
            hasXlabel = false;
            hasTitle = true;
        elseif strcmp(atom, 'Synergy')
            mean_v = squeeze(mean(PID_v(6, :,  :, :),4, 'omitnan')) - squeeze(GroundTruth(6,:,:));
            SEM  = squeeze(3 * (std(PID_v(6, :,  :, :),1,4,'omitnan') / sqrt(n_iter)));
            [nX, nLines] = size(mean_v);
            colorIndices = round(linspace(1, numColors, nLines));
            colorIndices = flip(colorIndices);
            linecolor = SynColormap(colorIndices, :);
            hasXlabel = false;
            hasTitle = false;
        elseif strcmp(atom, 'Redundancy')
            mean_v = squeeze(mean(PID_v(5, :,  :, :),4, 'omitnan')) - squeeze(GroundTruth(5,:,:));
            SEM  = squeeze(3 * (std(PID_v(5, :,  :, :),1,4,'omitnan') / sqrt(n_iter)));
            [nX, nLines] = size(mean_v);
            colorIndices = round(linspace(1, numColors, nLines));
            colorIndices = flip(colorIndices);
            linecolor = RedColormap(colorIndices, :);
            hasXlabel = true;
            hasTitle = false;
        elseif strcmp(atom, 'Unique')
            PID_v(2, :, :, :) = PID_v(2, :,  :, :)/2;
            mean_v = squeeze(mean(PID_v(2, :,  :, :),4, 'omitnan')) - (squeeze(GroundTruth(2,:,:))/2);
            SEM  = squeeze(3 * (std(PID_v(2, :,  :, :),1,4,'omitnan') / sqrt(n_iter)));
            [nX, nLines] = size(mean_v);
            colorIndices = round(linspace(1, numColors, nLines));
            colorIndices = flip(colorIndices);
            linecolor = UniqueColormap(colorIndices, :);
            hasXlabel = true;
            hasTitle = false;
        end
        categories = category;
        isTrialSweep = false;
    end


    if biasIdx == 1
        rowtxttrue = true;
    else
        rowtxttrue = false;
        if isGauss && length(trial_categories) == 1
            min_v = -10;
            max_v =10;
        elseif isGauss && length(trial_categories) > 1
            if strcmp(atom, 'Joint')
                if biasIdx > 1
                    min_v = -0.5;
                    max_v = 11;
                end
            elseif strcmp(atom, 'Synergy')
                if biasIdx > 1
                    min_v = -1;
                    max_v = 4;
                end
            elseif strcmp(atom, 'Redundancy')
                if biasIdx > 1
                    min_v = -1;
                    max_v = 4;
                end
            elseif strcmp(atom, 'Unique')
                if biasIdx > 1
                    min_v = -0.2;
                    max_v = 3;
                end
            end
        end
    end

    % Plot the means
    plotSubplot(mean_v, SEM, categories, linecolor, hasXlabel, fig_handle, subplot_idx, bias_correction, rowtxttrue, atom, min_v, max_v, hasTitle, isGauss, isTrialSweep, thumbRule); %-mean_v(end,:)
    subplot_idx = subplot_idx+1;
end
end




function plotSubplot(mean, SEM, categories, color, hasXlabel, fig_handle, subplot_idx, titletxt, rowtxttrue, rowtxt, min_val, max_val, hasTitle, isGauss, isTrialSweep, thumbRule)

nexttile(fig_handle, subplot_idx);
hold on;
FontSize = 18;
SEM = 0 * SEM;
set(gca, 'FontSize', FontSize, 'LineWidth',1.2);
[nX, nLines] = size(mean);
% Plot the means
hold on;
for i = 1:nLines
    plot(categories, mean(:, i), 'Color', color(i, :), 'LineWidth', 1.2);
    fill([categories(:); flipud(categories(:))],[mean(:, i) + SEM(:, i); flipud(mean(:, i) - SEM(:, i))], color(i, :), ...
        'FaceAlpha', 0.25,'EdgeColor', 'none');
end
xline(thumbRule, ':', 'Color',  'k', 'LineWidth', 2, 'HandleVisibility', 'off');
if hasTitle
    titletext = sprintf('\\bf%s\\rm', titletxt);
    title(titletext);
end
if isTrialSweep
    set(gca, 'XScale', 'log', 'LineWidth',1.2);
end
if isGauss
    if isTrialSweep
        xlim([64 2048]);
    else
        xlim([categories(1) categories(end)]);
    end
else
    if isTrialSweep
        xlim([4*16 4*2048]);
    else
        xlim([categories(1) categories(end)]);
    end
end
if hasXlabel
    if isTrialSweep
        xticks([100, 1000, 10000]);
        xticklabels({'10^2', '10^3', '10^4'});
        xlabel('trials', 'FontSize', FontSize);
    else
        if isGauss
            xlabel('dimensions', 'FontSize', FontSize);
        else
            xlabel('trials', 'FontSize', FontSize);
        end
    end
else
    xticklabels({});
end
if rowtxttrue
    subplot_text = rowtxt;
    subplot_text_bold = sprintf('\\bf%s\\rm', subplot_text);
    if isTrialSweep
        ylabelText = sprintf('%s\nInformation [bits]', subplot_text_bold);
        ylabel(ylabelText, 'FontSize', FontSize, 'FontWeight', 'normal');
    else
        ylabelText = sprintf('%s\nBias [bits]', subplot_text_bold);
        ylabel(ylabelText, 'FontSize', FontSize, 'FontWeight', 'normal');
    end
end
ylim([min_val, max_val])
hold off;
end