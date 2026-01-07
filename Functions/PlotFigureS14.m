function PlotFigureS13()
% MAKE_PID_BIAS_FIGURE (Weighted Merged)
% One publication-quality figure for Broja + Gauss results.
% Rows (top->bottom):  Broja Weighted Merged, Broja Plugin, Gauss Merged, Gauss Plugin
% Cols (left->right):  Joint, Syn, Red
%
% Saves PNG+PDF at 600 dpi and prints slope table to console.

%% --------------------------- Paths & Options ----------------------------
basePathBroja = 'Results/Broja/';
basePathGauss = 'Results/Gauss/';
savePath      = 'Figures_mat';   % change if you prefer a different folder

opts.fitType        = 'quadratic';   % 'linear0' (bx), 'linear' (bx+c), or 'quadratic' (ax^2+bx)
opts.weightedFit    = false;         % true = weight by avg SEM
opts.exportName     = 'Figure_S14';
opts.fontName       = 'Arial';   % use a journal-approved font if needed
opts.fontSize       = 8;             % axis tick/title font size
opts.labelSize      = 8;            % axis labels
opts.titleSize      = 8;
opts.lineWidth      = 1.4;
opts.markerSize     = 5;
opts.capSize        = 0;
opts.colors.low     = [0.23 0.49 0.77];   % low SNR color
opts.colors.high    = [0.83 0.24 0.31];   % high SNR color
opts.colors.fit     = [0.10 0.10 0.10];   % fit line color

% ---- Weighted merge controls ----
opts.useWeightedMerged = true;               % << Replace merged with weighted merged
% Weight as a function of information fraction in [0,1].
% Higher information -> trust shuffled (weight QE less). Linear ramp by default.
opts.weightFn = @(r) max(0,min(1, 1 - r));   % w=1 uses QE only, w=0 uses Shuff only

%% --------------------------- Graphics defaults --------------------------
set(0,'DefaultAxesFontName',opts.fontName);
set(0,'DefaultAxesFontSize',opts.fontSize);
set(0,'DefaultTextFontName',opts.fontName);
set(0,'DefaultLineLineWidth',opts.lineWidth);

%% ------------------------------ BROJA -----------------------------------
% Files to load
filesBroja = { ...
    'GroundTruth_bit_of_all_high.mat',        'gt_boa_high'; ...
    'GroundTruth_bit_of_all_low.mat',         'gt_boa_low';  ...
    'Simuldata_bit_of_all_high_I_BROJA.mat',  'pid_boa_high';...
    'Simuldata_bit_of_all_low_I_BROJA.mat',   'pid_boa_low'; ...
    'GroundTruth_unique_high_red_high.mat',   'gt_hr_high';  ...
    'GroundTruth_unique_high_red_low.mat',    'gt_hr_low';   ...
    'Simuldata_unique_high_red_high_I_BROJA.mat','pid_hr_high';...
    'Simuldata_unique_high_red_low_I_BROJA.mat', 'pid_hr_low'; ...
    'GroundTruth_unique_high_syn_high.mat',   'gt_hs_high';  ...
    'GroundTruth_unique_high_syn_low.mat',    'gt_hs_low';   ...
    'Simuldata_unique_high_syn_high_I_BROJA.mat','pid_hs_high';...
    'Simuldata_unique_high_syn_low_I_BROJA.mat', 'pid_hs_low'; ...
    'GroundTruth_uncorr_unique_high.mat',     'gt_uu_high';  ...
    'GroundTruth_uncorr_unique_low.mat',      'gt_uu_low';   ...
    'Simuldata_uncorr_unique_high_I_BROJA.mat','pid_uu_high';...
    'Simuldata_uncorr_unique_low_I_BROJA.mat', 'pid_uu_low'  ...
    };

B = load_broja(basePathBroja, filesBroja);

% Compute means/SEM for Broja
metrics = {'Joint','Syn','Red'}; %'Joint',
trialCatsBroja = [64,128,256,512,1024,2048]*4;
allTrialsBroja = [16,32,64,128,256,512,1024,2048]*4;
[broja, xBroja] = process_broja(B, metrics, trialCatsBroja, allTrialsBroja, opts);

%% ------------------------------ GAUSS -----------------------------------
filesGauss = { ...
    'Finalresults_across_M_and_ntrials.mat',              'pid_boa'; ...
    'Finalresults_across_M_and_ntrials_high_synergy.mat', 'pid_hs';  ...
    'Finalresults_across_M_and_ntrials_zero_synergy.mat', 'pid_zs';  ...
    'Finalresults_across_M_and_ntrials_both_unique.mat',  'pid_bu'   ...
    };
dim   = 20;
G     = load_gauss(basePathGauss, filesGauss, dim);

% NOTE: fixed a likely typo "2058" -> "2048"
trialCatsGauss = [256,512,1024,2058];
[gauss, xGauss] = process_gauss(G, metrics, trialCatsGauss, opts);

%% ------------------------------ Figure ----------------------------------
t = tiledlayout(4,3,'TileSpacing','compact','Padding','compact');
fig = gcf; set(fig,'Color','w','Units','centimeters'); fig.Position = [1 1 19 15]; % fits 1-col + legend

if false && opts.useWeightedMerged
    titles = {@(m)[m ' (Weighted Merged)'], @(m)[m ' (Plugin)']};
    rowLabels = {'Broja – Weighted Merged','Broja – Plugin','Gauss – Merged','Gauss – Plugin'};
else
    titles = {@(m)[m ' (Merged)'], @(m)[m ' (Plugin)']};
    rowLabels = {'Broja – Merged','Broja – Plugin','Gauss – Merged','Gauss – Plugin'};
end

% Panel 1–3: Broja Merged; 4–6: Broja Plugin; 7–9: Gauss Merged; 10–12: Gauss Plugin
plot_panel(broja.merged,  xBroja, 1, titles{1}, rowLabels{1}, metrics, opts);
plot_panel(broja.plugin,  xBroja, 2, titles{2}, rowLabels{2}, metrics, opts);
plot_panel(gauss.merged,  xGauss, 3, @(m)[m ' (Merged)'], rowLabels{3}, metrics, opts, 'gauss');
plot_panel(gauss.plugin,  xGauss, 4, titles{2}, rowLabels{4}, metrics, opts, 'gauss');

% Super labels
% xlabel(t, '(bins)^3 / nTrials   (top 2 rows)      |      3·dim / nTrials   (bottom 2 rows)','FontSize',opts.labelSize,'FontWeight','bold');
% ylabel(t, 'Bias','FontSize',opts.labelSize,'FontWeight','bold');

% Optional overall title:
% title(t, 'PID Bias vs. Sampling Ratio','FontSize',opts.titleSize+1,'FontWeight','bold');

%% ----------------------------- Slopes (table) ---------------------------
fprintf('\n### Linear (no-intercept) slopes from mean of Low/High\n');
fprintf('| Dataset | Method  | Syn     | Red     |\n');
fprintf('|---------|---------|---------|---------|\n');
sBrojaM = slope_table(xBroja, broja.merged,  metrics);
sBrojaP = slope_table(xBroja, broja.plugin,  metrics);
sGaussM = slope_table(xGauss, gauss.merged,  metrics);
sGaussP = slope_table(xGauss, gauss.plugin,  metrics);
fprintf('| Broja   | Merged  | %.4fx | %.4fx |\n', sBrojaM.Syn, sBrojaM.Red);
fprintf('| Broja   | Plugin  | %.4fx | %.4fx |\n', sBrojaP.Syn, sBrojaP.Red);
fprintf('| Gauss   | Merged  | %.4fx | %.4fx |\n', sGaussM.Syn, sGaussM.Red);
fprintf('| Gauss   | Plugin  | %.4fx | %.4fx |\n', sGaussP.Syn, sGaussP.Red);

%% ------------------------------ Export ----------------------------------
% PNG + PDF @ 600 dpi
set(gcf, 'Renderer', 'painters');
print(gcf, '-dsvg', opts.exportName);
print(fig, fullfile(savePath,[opts.exportName '.svg']), '-dsvg','-r600');


end

%% ============================== HELPERS =================================
function S = load_broja(basePath, files)
for i = 1:size(files,1)
    tmp = load(fullfile(basePath, files{i,1}));
    S.(files{i,2}) = tmp;
end
end

function [out, xvals] = process_broja(S, metrics, trialCats, allTrials, opts)
% Extract needed arrays
gt_boa_high = S.gt_boa_high.GroundTruth_value;  gt_boa_low  = S.gt_boa_low.GroundTruth_value;
gt_hr_high  = S.gt_hr_high.GroundTruth_value;   gt_hr_low   = S.gt_hr_low.GroundTruth_value;
gt_hs_high  = S.gt_hs_high.GroundTruth_value;   gt_hs_low   = S.gt_hs_low.GroundTruth_value;
gt_uu_high  = S.gt_uu_high.GroundTruth_value;   gt_uu_low   = S.gt_uu_low.GroundTruth_value;

% Available simulated PID (plugin always present)
pid_boa_high_pl = S.pid_boa_high.PID_v_plugin;  pid_boa_high_mg = S.pid_boa_high.PID_v_qeShuff;
pid_boa_low_pl  = S.pid_boa_low.PID_v_plugin;   pid_boa_low_mg  = S.pid_boa_low .PID_v_qeShuff;
pid_hr_high_pl  = S.pid_hr_high.PID_v_plugin;   pid_hr_high_mg  = S.pid_hr_high .PID_v_qeShuff;
pid_hr_low_pl   = S.pid_hr_low.PID_v_plugin;    pid_hr_low_mg   = S.pid_hr_low  .PID_v_qeShuff;
pid_hs_high_pl  = S.pid_hs_high.PID_v_plugin;   pid_hs_high_mg  = S.pid_hs_high .PID_v_qeShuff;
pid_hs_low_pl   = S.pid_hs_low.PID_v_plugin;    pid_hs_low_mg   = S.pid_hs_low  .PID_v_qeShuff;
pid_uu_high_pl  = S.pid_uu_high.PID_v_plugin;   pid_uu_high_mg  = S.pid_uu_high .PID_v_qeShuff;
pid_uu_low_pl   = S.pid_uu_low.PID_v_plugin;    pid_uu_low_mg   = S.pid_uu_low  .PID_v_qeShuff;

% Optional: if QE and Shuff are provided separately, use them for weighted merge
hasQE  = isfield(S.pid_boa_high,'PID_v_qe');
hasSH  = isfield(S.pid_boa_high,'PID_v_shuff');

map = struct('Joint',1:4,'Syn',4,'Red',1);

nT = numel(trialCats); nM = numel(metrics);
out.merged.mean_low  = zeros(nT,nM);
out.merged.mean_high = zeros(nT,nM);
out.merged.sem_low   = zeros(nT,nM);
out.merged.sem_high  = zeros(nT,nM);
out.plugin = out.merged;

xvals = zeros(nT,1);

for it = 1:nT
    trial = trialCats(it);
    xvals(it) = 4^3 / trial;                % (bins)^3 / nTrials
    tidx = find(allTrials==trial,1);

    for m = 1:nM
        idx = map.(metrics{m});

        o_l_boa = sum(gt_boa_low(idx));  o_h_boa = sum(gt_boa_high(idx));
        o_l_hr  = sum(gt_hr_low(idx));   o_h_hr  = sum(gt_hr_high(idx));
        o_l_hs  = sum(gt_hs_low(idx));   o_h_hs  = sum(gt_hs_high(idx));
        o_l_uu  = sum(gt_uu_low(idx));   o_h_uu  = sum(gt_uu_high(idx));

        % ---------------- merged (weighted if available) -----------------
        if false && opts.useWeightedMerged && hasQE && hasSH
            % pull QE/SHUFF for metric m at this trial across seeds
            [qe_l_boa, sh_l_boa] = get_qe_sh(S.pid_boa_low,  m, tidx);
            [qe_h_boa, sh_h_boa] = get_qe_sh(S.pid_boa_high, m, tidx);
            [qe_l_hr,  sh_l_hr ] = get_qe_sh(S.pid_hr_low,   m, tidx);
            [qe_h_hr,  sh_h_hr ] = get_qe_sh(S.pid_hr_high,  m, tidx);
            [qe_l_hs,  sh_l_hs ] = get_qe_sh(S.pid_hs_low,   m, tidx);
            [qe_h_hs,  sh_h_hs ] = get_qe_sh(S.pid_hs_high,  m, tidx);
            [qe_l_uu,  sh_l_uu ] = get_qe_sh(S.pid_uu_low,   m, tidx);
            [qe_h_uu,  sh_h_uu ] = get_qe_sh(S.pid_uu_high,  m, tidx);

            % compute a single weight per trial using JOINT (metric 'Joint')
            jidx = map.Joint;  % 1:4 (sum)
            max_info_low  = sum(gt_boa_low(jidx));
            max_info_high = sum(gt_boa_high(jidx));
            % joint values (averaged QE/SH)
            [qeLJ, shLJ] = get_qe_sh(S.pid_boa_low,  1, tidx);  % Joint
            [qeHJ, shHJ] = get_qe_sh(S.pid_boa_high, 1, tidx);
            wL = opts.weightFn( mean((qeLJ + shLJ)/2, 'omitnan') / max_info_low );
            wH = opts.weightFn( mean((qeHJ + shHJ)/2, 'omitnan') / max_info_high );

            d_l = [ mean(wL*qe_l_boa + (1-wL)*sh_l_boa - o_l_boa,'omitnan'), ...
                    mean(wL*qe_l_hr  + (1-wL)*sh_l_hr  - o_l_hr , 'omitnan'), ...
                    mean(wL*qe_l_hs  + (1-wL)*sh_l_hs  - o_l_hs , 'omitnan'), ...
                    mean(wL*qe_l_uu  + (1-wL)*sh_l_uu  - o_l_uu , 'omitnan') ];
            d_h = [ mean(wH*qe_h_boa + (1-wH)*sh_h_boa - o_h_boa,'omitnan'), ...
                    mean(wH*qe_h_hr  + (1-wH)*sh_h_hr  - o_h_hr , 'omitnan'), ...
                    mean(wH*qe_h_hs  + (1-wH)*sh_h_hs  - o_h_hs , 'omitnan'), ...
                    mean(wH*qe_h_uu  + (1-wH)*sh_h_uu  - o_h_uu , 'omitnan') ];
        else
            % fallback: use provided merged (PID_v_qeShuff)
            d_l = [ mean(squeeze(pid_boa_low_mg (tidx,m,:)) - o_l_boa,'omitnan'), ...
                mean(squeeze(pid_hr_low_mg  (tidx,m,:)) - o_l_hr , 'omitnan'), ...
                mean(squeeze(pid_hs_low_mg  (tidx,m,:)) - o_l_hs , 'omitnan'), ...
                mean(squeeze(pid_uu_low_mg  (tidx,m,:)) - o_l_uu , 'omitnan') ];
            d_h = [ mean(squeeze(pid_boa_high_mg(tidx,m,:)) - o_h_boa,'omitnan'), ...
                mean(squeeze(pid_hr_high_mg (tidx,m,:)) - o_h_hr , 'omitnan'), ...
                mean(squeeze(pid_hs_high_mg (tidx,m,:)) - o_h_hs , 'omitnan'), ...
                mean(squeeze(pid_uu_high_mg (tidx,m,:)) - o_h_uu , 'omitnan') ];
        end

        out.merged.mean_low (it,m) = mean(d_l,'omitnan');
        out.merged.mean_high(it,m) = mean(d_h,'omitnan');
        out.merged.sem_low  (it,m) = 2*std(d_l,'omitnan')/sqrt(numel(d_l));
        out.merged.sem_high (it,m) = 2*std(d_h,'omitnan')/sqrt(numel(d_h));

        % ---------------- plugin ----------------
        p_l = [ mean(squeeze(pid_boa_low_pl (tidx,m,:)) - o_l_boa,'omitnan'), ...
            mean(squeeze(pid_hr_low_pl  (tidx,m,:)) - o_l_hr , 'omitnan'), ...
            mean(squeeze(pid_hs_low_pl  (tidx,m,:)) - o_l_hs , 'omitnan'), ...
            mean(squeeze(pid_uu_low_pl  (tidx,m,:)) - o_l_uu , 'omitnan') ];
        p_h = [ mean(squeeze(pid_boa_high_pl(tidx,m,:)) - o_h_boa,'omitnan'), ...
            mean(squeeze(pid_hr_high_pl (tidx,m,:)) - o_h_hr , 'omitnan'), ...
            mean(squeeze(pid_hs_high_pl (tidx,m,:)) - o_h_hs , 'omitnan'), ...
            mean(squeeze(pid_uu_high_pl (tidx,m,:)) - o_h_uu , 'omitnan') ];

        out.plugin.mean_low (it,m) = mean(p_l,'omitnan');
        out.plugin.mean_high(it,m) = mean(p_h,'omitnan');
        out.plugin.sem_low  (it,m) = 2*std(p_l,'omitnan')/sqrt(numel(p_l));
        out.plugin.sem_high (it,m) = 2*std(p_h,'omitnan')/sqrt(numel(p_h));
    end
end

    function [qe,sh] = get_qe_sh(Snode, mIdx, tIdx)
        qe = squeeze(Snode.PID_v_qe    (tIdx,mIdx,:));
        sh = squeeze(Snode.PID_v_shuff (tIdx,mIdx,:));
    end
end

function G = load_gauss(basePath, files, dim)
for i = 1:size(files,1)
    tmp = load(fullfile(basePath, files{i,1}));
    G.(files{i,2}) = tmp;
end
G.dim = dim;
end

function [out, xvals] = process_gauss(G, metrics, trialCats, opts)
% Pull pieces
M_vals = G.pid_boa.M_vals;
dimidx = find(M_vals==G.dim,1);

% Helper to extract sets (GT + sampled)
getSet = @(S) deal( ...
    squeeze(S.GT_results(:, 1, dimidx, 4, 1)), ... % high GT
    squeeze(S.GT_results(:, 1, dimidx, 4, 3)), ... % low  GT
    squeeze(S.sampled_results(:,:,dimidx,1,1,:)), ... % plugin high
    squeeze(S.sampled_results(:,:,dimidx,2,1,:)), ... % QE     high
    squeeze(S.sampled_results(:,:,dimidx,3,1,:)), ... % Shuff  high
    squeeze(S.sampled_results(:,:,dimidx,1,3,:)), ... % plugin low
    squeeze(S.sampled_results(:,:,dimidx,2,3,:)), ... % QE     low
    squeeze(S.sampled_results(:,:,dimidx,3,3,:))  ... % Shuff  low
    );

[gt_boa_h, gt_boa_l, pl_boa_h, qe_boa_h, sh_boa_h, pl_boa_l, qe_boa_l, sh_boa_l] = getSet(G.pid_boa);
[gt_hs_h,  gt_hs_l,  pl_hs_h,  qe_hs_h,  sh_hs_h,  pl_hs_l,  qe_hs_l,  sh_hs_l ] = getSet(G.pid_hs);
[gt_zs_h,  gt_zs_l,  pl_zs_h,  qe_zs_h,  sh_zs_h,  pl_zs_l,  qe_zs_l,  sh_zs_l ] = getSet(G.pid_zs);
[gt_bu_h,  gt_bu_l,  pl_bu_h,  qe_bu_h,  sh_bu_h,  pl_bu_l,  qe_bu_l,  sh_bu_l ] = getSet(G.pid_bu);

map13 = struct('Joint',1,'Syn',6,'Red',5);
map14 = struct('Joint',1,'Syn',7,'Red',6);

nT = numel(trialCats); nM = numel(metrics);
% Preallocate
z = zeros(nT,nM);
out.merged.mean_low  = z; out.merged.mean_high = z;
out.merged.sem_low   = z; out.merged.sem_high  = z;
out.plugin           = out.merged;

xvals = zeros(nT,1);

for it = 1:nT
    trial = trialCats(it);
    xvals(it) = 3*G.dim / trial;                     % 3·dim / nTrials
    tidx = find(G.pid_boa.ntrials_vals==trial,1);

    for m = 1:nM
        idx13 = map13.(metrics{m});
        idx14 = map14.(metrics{m});

        % Ground-truth offsets
        off_l = [sum(gt_boa_l(idx13)), sum(gt_zs_l(idx14)), sum(gt_hs_l(idx14)), sum(gt_bu_l(idx13))];
        off_h = [sum(gt_boa_h(idx13)), sum(gt_zs_h(idx14)), sum(gt_hs_h(idx14)), sum(gt_bu_h(idx13))];

        % ----- merged: weighted or simple average -----
        if false && opts.useWeightedMerged
            % weight determined from JOINT metric (idx13/idx14 Joint=1)
            max_info_low  = sum(gt_boa_l(map13.Joint));
            max_info_high = sum(gt_boa_h(map13.Joint));

            % joint at this trial across seeds
            qeLJ = squeeze(qe_boa_l(map13.Joint, tidx, :));
            shLJ = squeeze(sh_boa_l(map13.Joint, tidx, :));
            qeHJ = squeeze(qe_boa_h(map13.Joint, tidx, :));
            shHJ = squeeze(sh_boa_h(map13.Joint, tidx, :));

            wL = opts.weightFn( mean((qeLJ + shLJ)/2,'omitnan') / max_info_low );
            wH = opts.weightFn( mean((qeHJ + shHJ)/2,'omitnan') / max_info_high );

            d_l = collect_weighted_means(qe_boa_l, sh_boa_l, qe_zs_l, sh_zs_l, qe_hs_l, sh_hs_l, qe_bu_l, sh_bu_l, idx13, idx14, tidx, off_l, wL);
            d_h = collect_weighted_means(qe_boa_h, sh_boa_h, qe_zs_h, sh_zs_h, qe_hs_h, sh_hs_h, qe_bu_h, sh_bu_h, idx13, idx14, tidx, off_h, wH);
        else
            % original equal-mean merge
            d_l = collect_means((qe_boa_l+sh_boa_l)/2, (qe_zs_l+sh_zs_l)/2, (qe_hs_l+sh_hs_l)/2, (qe_bu_l+sh_bu_l)/2, idx13, idx14, tidx, off_l);
            d_h = collect_means((qe_boa_h+sh_boa_h)/2, (qe_zs_h+sh_zs_h)/2, (qe_hs_h+sh_hs_h)/2, (qe_bu_h+sh_bu_h)/2, idx13, idx14, tidx, off_h);
        end

        out.merged.mean_low (it,m) = mean(d_l,'omitnan');
        out.merged.mean_high(it,m) = mean(d_h,'omitnan');
        out.merged.sem_low  (it,m) = 2*std(d_l,'omitnan')/sqrt(numel(d_l));
        out.merged.sem_high (it,m) = 2*std(d_h,'omitnan')/sqrt(numel(d_h));

        % ----- plugin -----
        p_l = collect_means(pl_boa_l, pl_zs_l, pl_hs_l, pl_bu_l, idx13, idx14, tidx, off_l);
        p_h = collect_means(pl_boa_h, pl_zs_h, pl_hs_h, pl_bu_h, idx13, idx14, tidx, off_h);

        out.plugin.mean_low (it,m) = mean(p_l,'omitnan');
        out.plugin.mean_high(it,m) = mean(p_h,'omitnan');
        out.plugin.sem_low  (it,m) = 2*std(p_l,'omitnan')/sqrt(numel(p_l));
        out.plugin.sem_high (it,m) = 2*std(p_h,'omitnan')/sqrt(numel(p_h));
    end
end

% ---- local helpers ----
    function v = collect_means(A1,A2,A3,A4, idx13, idx14, tidx, off)
        v = [ mean(squeeze(A1(idx13,tidx,:)) - off(1), 'omitnan'), ...
            mean(squeeze(A2(idx14,tidx,:)) - off(2), 'omitnan'), ...
            mean(squeeze(A3(idx14,tidx,:)) - off(3), 'omitnan'), ...
            mean(squeeze(A4(idx13,tidx,:)) - off(4), 'omitnan') ];
    end
    function v = collect_weighted_means(qe1,sh1, qe2,sh2, qe3,sh3, qe4,sh4, idx13, idx14, tidx, off, w)
        v = [ mean(w*squeeze(qe1(idx13,tidx,:)) + (1-w)*squeeze(sh1(idx13,tidx,:)) - off(1), 'omitnan'), ...
              mean(w*squeeze(qe2(idx14,tidx,:)) + (1-w)*squeeze(sh2(idx14,tidx,:)) - off(2), 'omitnan'), ...
              mean(w*squeeze(qe3(idx14,tidx,:)) + (1-w)*squeeze(sh3(idx14,tidx,:)) - off(3), 'omitnan'), ...
              mean(w*squeeze(qe4(idx13,tidx,:)) + (1-w)*squeeze(sh4(idx13,tidx,:)) - off(4), 'omitnan') ];
    end
end

function plot_panel(D, xvals, panelRow, titleFun, rowLabel, metrics, opts, which)
% Draw a 1×3 row of subplots (for a given dataset/method)
% D has fields: mean_low, mean_high, sem_low, sem_high
if nargin < 8, which = 'broja'; end
xvals = xvals(:);
for m = 1:numel(metrics)
    nexttile((panelRow-1)*3 + m); hold on;

    % sort by x for monotonic lines
    [xs, ix] = sort(xvals);
    yL = abs(D.mean_low (ix,m));  eL = D.sem_low (ix,m);
    yH = abs(D.mean_high(ix,m));  eH = D.sem_high(ix,m);

    % data points
    hL = errorbar(xs, yL, eL, '-o', 'MarkerSize',opts.markerSize, 'CapSize',opts.capSize, ...
        'Color',opts.colors.low,  'MarkerFaceColor',opts.colors.low*0.85);
    hH = errorbar(xs, yH, eH, '-s', 'MarkerSize',opts.markerSize, 'CapSize',opts.capSize, ...
        'Color',opts.colors.high, 'MarkerFaceColor',opts.colors.high*0.85);

    % fit on average of Low/High
    yAll = (yL + yH)/2;
    eAll = (eL + eH)/2;
    [yfit, fitStr] = fit_curve(xs, yAll, eAll, opts.fitType, opts.weightedFit);
    xf = linspace(min(xs), max(xs), 300);
    switch lower(opts.fitType)
        case 'quadratic'
            yF = yfit(1)*xf.^2 + yfit(2)*xf;
        case 'linear0'
            yF = yfit(1)*xf;
        otherwise % 'linear'
            yF = yfit(1)*xf + yfit(2);
    end
    hF = plot(xf, yF, '-', 'Color',opts.colors.fit);

    % labels & cosmetics
    if strcmpi(which,'broja')
        xlabel('(bins)^3 / nTrials','FontSize',opts.labelSize);
    else
        xlabel('3·dim / nTrials','FontSize',opts.labelSize);
    end
    ylabel('Bias','FontSize',opts.labelSize);
    title(titleFun(metrics{m}), 'FontSize',opts.titleSize,'FontWeight','bold');

    grid on; box on;
    set(gca,'LineWidth',0.8,'TickDir','out');

    % Row tag on first subplot of the row
    % if m==1
    %     yl = ylim; xl = xlim;
    %     text(xl(1), yl(2) + 0.06*range(yl), rowLabel, ...
    %         'HorizontalAlignment','left','VerticalAlignment','bottom', ...
    %         'FontWeight','bold','FontSize',opts.titleSize);
    % end

    % Legend only on middle subplot of the row
    if m==2
        legend([hL hH hF], {'Low','High',['Fit: ' fitStr]}, ...
            'Location','best','Box','off','FontSize',opts.fontSize);
    else
        % --- write fit formula on every panel ---
        txt = ['Fit: ' fitStr];
        text(0.03, 0.95, txt, ...
            'Units','normalized', ...
            'HorizontalAlignment','left', ...
            'VerticalAlignment','top', ...
            'FontSize',opts.fontSize, ...
            'Color',opts.colors.fit, ...
            'Interpreter','none', ...
            'BackgroundColor','w', ...
            'Margin',2, ...
            'EdgeColor','none');
    end
end
end

function [coeffs, fitStr] = fit_curve(x, y, e, fitType, weighted)
x = x(:); y = y(:); e = e(:);
switch lower(fitType)
    case 'quadratic'
        X = [x.^2 x];      % ax^2 + bx
    case 'linear0'
        X = x;             % bx
    otherwise
        X = [x ones(size(x))]; % bx + c
end
if weighted
    w = 1./max(e, eps).^2;
    coeffs = (X'*(w.*X)) \ (X'*(w.*y));
else
    coeffs = X \ y;
end

switch lower(fitType)
    case 'quadratic', fitStr = sprintf('%.3fx^2 + %.3fx', coeffs(1), coeffs(2));
    case 'linear0',  fitStr = sprintf('%.3fx', coeffs(1));
    otherwise,       fitStr = sprintf('%.3fx + %.3f', coeffs(1), coeffs(2));
end
end

function s = slope_table(x, D, metrics)
% Linear no-intercept fit to avg(Low, High) for Syn/Red
x = x(:);
idxSyn = find(strcmp(metrics,'Syn'));
idxRed = find(strcmp(metrics,'Red'));

avgSyn = (abs(D.mean_low(:,idxSyn)) + abs(D.mean_high(:,idxSyn)))/2;
avgRed = (abs(D.mean_low(:,idxRed)) + abs(D.mean_high(:,idxRed)))/2;

b_syn = (x'*(x)) \ (x'*avgSyn);
b_red = (x'*(x)) \ (x'*avgRed);

s = struct('Syn',b_syn,'Red',b_red);
end

%% ===================== Optional utility (vectorized) =====================
function [weighted_mat, w_used] = apply_adaptive_weighting(qe_mat, shuff_mat, max_info, get_weight)
    % Applies an adaptive weight between QE and Shuffled estimates using
    % the JOINT column (assumed to be column 1) to compute the information
    % fraction. Each row is a trial setting; columns are metrics, rows/cols
    % may be squeezed depending on usage. Provided for completeness. Not
    % used directly in the main flow where we operate at per-trial scalars.
    if isempty(qe_mat)
        weighted_mat = qe_mat; w_used = []; return;
    end
    joint_qe    = qe_mat(:,1);
    joint_shuff = shuff_mat(:,1);
    joint_avg   = (joint_qe + joint_shuff)/2;

    info_frac = joint_avg ./ max_info;
    w_used = arrayfun(get_weight, info_frac);

    weighted_mat = w_used .* qe_mat + (1 - w_used) .* shuff_mat;
end
