%% === Config ===
trial_categories   = [16, 32, 64, 128, 256, 512, 1024, 2048];
bias_correction    = {'plugin','qe','shuff','qeShuff','infoCorr'};
simul_cases        = {'uncorr_unique'};   % add more cases to get more rows
info_amount        = 'low';
ResultFolder       = 'SVM_results';
redundancy_measure = 'I_BROJA';

info_types = {'Joint','Syn','Red'};   % columns
nrows = numel(simul_cases);
ncols = numel(info_types);
fill_alpha = 0.22;                    % transparency for SEM bands

%% === Figure / layout ===
figure_handle = figure('Units','inches','Position',[1, 1, 14, 3+2*(nrows-1)]);
tl = tiledlayout(figure_handle, nrows, ncols, 'TileSpacing','compact','Padding','compact');

%% === Colors per bias method (keyed by canonical names) ===
cmap = containers.Map( ...
  {'plugin','qe','resample','resampling','infocorr','shuff','shuffle','qe_shuff','qeshuff', ...
   'merged','shuff-resamp','weighted','venkatesh'}, ...
  num2cell([ ...
    0.60,  0.60,  0.60;   % plugin (gray)
    0.00,  0.447, 0.741;  % qe (blue)
    0.929, 0.694, 0.125;  % resample (yellow)
    0.929, 0.694, 0.125;  % resampling (yellow)
    0.929, 0.694, 0.125;  % infocorr (yellow)  <-- add this
    0.850, 0.325, 0.098;  % shuff (orange)
    0.850, 0.325, 0.098;  % shuffle (orange)
    0.494, 0.184, 0.556;  % qe_shuff (purple)
    0.494, 0.184, 0.556;  % qeshuff (purple)
    0.494, 0.184, 0.556;  % merged (purple)
    0.494, 0.184, 0.556;  % shuff-resamp (purple)
    0.60,  0.60,  0.60;   % weighted (gray)
    0.466, 0.674, 0.188;  % venkatesh (green)
  ],2) ...
);


% Legend labels in the same order as bias_correction:
legendLabels = cellfun(@prettyMethodLabel, bias_correction, 'UniformOutput', false);

%% === Plot ===
for r = 1:nrows
    simul_case = simul_cases{r};

    % Preload stats for all methods once per row
    stats = struct();
    for m = 1:numel(bias_correction)
        bc = bias_correction{m};
        stats.(bc) = getPIDStats(trial_categories, simul_case, bc, info_amount, redundancy_measure, ResultFolder);
    end

    for c = 1:ncols
        nexttile(tl, (r-1)*ncols + c); hold on;

        % --- draw fills first (kept out of legend) ---
        for m = 1:numel(bias_correction)
            bc = bias_correction{m};
            S  = stats.(bc);
            switch info_types{c}
                case 'Joint', y = S.Joint; e = S.SEM_Joint;
                case 'Syn',   y = S.Syn;   e = S.SEM_Syn;
                case 'Red',   y = S.Red;   e = S.SEM_Red;
                otherwise,    error('Unknown info type: %s', info_types{c});
            end
            key = lower(bc);
            shadeSEM(S.trials, y, e, cmap(key), fill_alpha);
        end

        % --- overlay mean lines (these will appear in legend) ---
        Lh = gobjects(1, numel(bias_correction));
        for m = 1:numel(bias_correction)
            bc = bias_correction{m};
            S  = stats.(bc);
            switch info_types{c}
                case 'Joint', y = S.Joint;
                case 'Syn',   y = S.Syn;
                case 'Red',   y = S.Red;
            end
            key = lower(bc);
            Lh(m) = plot(S.trials, y, 'LineWidth', 1.6, 'Color', cmap(key), ...
                'DisplayName', legendLabels{m});
            plot(S.trials, y(end)*ones(size(S.trials)), 'k--','LineWidth',1)
        end

        
        % Axes cosmetics
        set(gca, 'XScale','log', 'LineWidth',1.2);
        xticks([100 1000 10000]); xticklabels({'10^2','10^3','10^4'});
        if strcmpi(ResultFolder,'SVM_results')
            xlim([4*trial_categories(1), 4*trial_categories(end)]);
        else
            xlim([trial_categories(1), trial_categories(end)]);
        end
        % ylim([0, 1.])
        yline(0,'k-','LineWidth',1.0,'HandleVisibility','off');

        if r == 1
            title(info_types{c}, 'FontWeight','bold');
        end
        if c == 1
            ylabel('Information [bits]');
        else
            ylim([-.05, .45])
        end
        if r == nrows
            xlabel('Trials (N)');
        end

        % legend only once (top-right panel)
        if r == 1 && c == ncols
            lgd = legend(Lh, legendLabels, 'Box','off', 'Location','northeast');
            lgd.ItemTokenSize = [12 6];
        end

        hold off;
    end
end

%% === Save ===
outDir = fullfile(fileparts(pwd), 'Figures_mat');  

if ~exist(outDir, 'dir')
    mkdir(outDir);
end

exportgraphics(gcf, fullfile(outDir, 'Figure_S17.svg'), 'ContentType','vector');
%% ------------------ Helpers ------------------
function S = getPIDStats(trial_categories, simul_case, bias_correction, info_amount, redundancy_measure, ResultFolder)
    % Loads PID_v_<bias> and returns mean curves + SEM for components.
    trials_plot = trial_categories;
    if strcmpi(ResultFolder,'SVM_results')
        trials_plot = trial_categories * 4;  % matches original SVM scaling
    end

    filename = sprintf('../Results/%s/SVM_results/Simuldata_%s_%s_%s_10_alpha1.mat', ...
        ResultFolder, simul_case, info_amount, redundancy_measure);
    varname = sprintf('PID_v_%s', bias_correction);

    L = load(filename, varname);
    if ~isfield(L, varname)
        error('Variable %s not found in %s.', varname, filename);
    end
    PID_v = L.(varname);  % dims: trials x components x iterations

    % Means & SEM across iterations (NaNs ignored)
    [S.Joint, S.SEM_Joint] = meanSEM_3D(PID_v(:,1,:));
    [S.Syn,   S.SEM_Syn]   = meanSEM_3D(PID_v(:,2,:));
    [S.Red,   S.SEM_Red]   = meanSEM_3D(PID_v(:,3,:));

    % Unq available if you add that column back later:
    UNQ_per_iter = squeeze(mean(PID_v(:,4:5,:), 2));   % T x K
    [S.Unq,  S.SEM_Unq]    = meanSEM_2D(UNQ_per_iter);

    S.trials = trials_plot(:)';                        % row vector for plotting
end

function [mu, sem] = meanSEM_3D(A)
    % A: T x 1 x K (or T x K x 1); aggregate along iterations.
    B = squeeze(A);                % -> T x K (or T if K==1)
    [mu, sem] = meanSEM_2D(B);
end

function [mu, sem] = meanSEM_2D(B)
    % B: T x K (trials x iterations), NaNs allowed.
    if isvector(B), B = B(:); end
    K_eff = sum(~isnan(B), 2);
    mu    = mean(B, 2, 'omitnan');
    s     = std(B, 0, 2, 'omitnan');   % N-1 normalization
    sem   = s ./ sqrt(max(K_eff,1));
    mu    = mu.';   % return row for plotting
    sem   = sem.';  % row
end

function shadeSEM(x, y, e, colorRGB, alphaVal)
    % Draws a filled SEM band using fill(); hidden from legend.
    yU = y + 3 * e; 
    yL = y - 3 * e;
    bad = isnan(yU) | isnan(yL) | isnan(x);
    yU(bad) = y(bad);
    yL(bad) = y(bad);
    fill([x, fliplr(x)], [yU, fliplr(yL)], colorRGB, ...
        'FaceAlpha', alphaVal, 'EdgeColor','none', 'HandleVisibility','off');
end

function label = prettyMethodLabel(bc)
    switch lower(bc)
        case 'plugin',    label = 'plugin';
        case 'qe',        label = 'QE';
        case 'shuff',     label = 'shuffle';
        case 'qeshuff',   label = 'merged';
        case 'infocorr',  label = 'resample';
        otherwise,        label = bc;
    end
end
