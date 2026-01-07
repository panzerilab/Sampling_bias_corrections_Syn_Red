function PlotFigure2discr()
% One row, three panels (DISCRETE, plugin):
% 1) bias vs trials (bins=4, alpha=3)   — Broja/high
% 2) bias vs bins   (trials=512)        — Bin_sweep/low  [matches your bin code]
% 3) bias vs alpha  (trials=512, bins=4) — Broja/per-alpha files
%
% Bias = mean(measured) - GT, where GT is taken using the same
% conventions as your bin plot:
% - GroundTruth_value is a cell: we vertcat and use component positions:
%   Joint = sum(1:4), Syn = 4, Red = 1.

% ---------------- Visuals ----------------
% Turn off LaTeX everywhere (tex = default math-ish, none = literal)
set(groot,'defaultTextInterpreter','tex');
set(groot,'defaultAxesTickLabelInterpreter','tex');
set(groot,'defaultLegendInterpreter','tex');

% Or if you want literal (no TeX parsing), use 'none' for all three lines above.

% Force a common font that Inkscape will have
set(groot,'defaultAxesFontName','Arial');
set(groot,'defaultTextFontName','Arial');
set(groot,'defaultLegendFontName','Arial');

set(0,'DefaultTextFontName','Arial'); 
set(0,'DefaultAxesFontName','Arial'); 
set(0,'DefaultLegendFontName','Arial');
FontSize = 9; lw = 1;
cols = struct('Joint',[0.9290,0.6940,0.1250], ... % yellow
              'Syn',  [0.4660,0.6740,0.1880], ... % green
              'Red',  [0.0000,0.4470,0.7410]);    % blue

aFill = 0.18;
alphatrials = 16;
fig = figure('Units','centimeters','Position',[1 1 18 6]);
tl  = tiledlayout(fig,1,3,'TileSpacing','compact','Padding','compact');



% ---------------- Panel 1: Trials (bins=4, alpha=3) ----------------
ax1 = nexttile; hold(ax1,'on');
% Files (Broja/high), plugin arrays: [trials x comps x iters]
F_data = 'Results/Broja/Simuldata_bit_of_all_low_I_BROJA.mat';
F_gt   = 'Results/Broja/GroundTruth_bit_of_all_low.mat';
assert_exist(F_data); assert_exist(F_gt);

% load plugin
L = load(F_data,'PID_v_plugin'); 
PID = L.PID_v_plugin;                 % T x C x I
T = size(PID,1);
comp = @(c) squeeze(PID(:,c,:));      % T x I

% GT matching bins=4, alpha=3 using same approach as your bin plot
G = load(F_gt,'GroundTruth_value');
GTv = G.GroundTruth_value;            % cell or numeric
GT_vec = pick_GT_likeBinCode(GTv, 4, 3);  % 1x5
GT_joint = sum(GT_vec(1:4)); GT_syn = GT_vec(4); GT_red = GT_vec(1);

% mean±3*SEM (SEM over iters), Bias = mean(meas) - GT
[m_joint1, s_joint1] = m_sem3(comp(1));   % Joint measured vs trials
[m_syn1,   s_syn1  ] = m_sem3(comp(2));
[m_red1,   s_red1  ] = m_sem3(comp(3));
x_trials = 4* map_trials_axis(T);            % [16 32 64 128 256 512 ...]
plot_fill(ax1, x_trials, m_joint1 - GT_joint, s_joint1, cols.Joint, lw, aFill);
plot_fill(ax1, x_trials, m_syn1   - GT_syn,   s_syn1,   cols.Syn,   lw, aFill);
plot_fill(ax1, x_trials, m_red1   - GT_red,   s_red1,   cols.Red,   lw, aFill);
set(ax1,'XScale','log','FontSize',FontSize,'LineWidth',lw);
xticks([100, 1000])
xlabel(ax1,'Trials (N)'); ylabel(ax1,'Bias [bits]'); yline(ax1,0,'k-','LineWidth',lw,'HandleVisibility','off');
title(ax1,'Trials (R=4)');
legend({'Joint','Syn','Red'},'Box','off','Location','best');
ylim([0,.5])

% ---------------- Panel 2: Bins (trials=512) ----------------
ax2 = nexttile; hold(ax2,'on');
% EXACTLY mirror your bin-plot handling: Bin_sweep/low; trialIdx mapping
ResultFolder = 'Results/Bin_sweep';
F2_data = fullfile(ResultFolder,'Simuldata_bit_of_all_low_I_BROJA.mat');
F2_gt   = fullfile(ResultFolder,'GroundTruth_bit_of_all_low.mat');
assert_exist(F2_data); assert_exist(F2_gt);

% Load: trials x comps x iters x cats  -> we want a chosen trial slice
trialIdx = 2;          % your original setting (maps to N=512)
X = load(F2_data, 'PID_v_plugin');           % T x C x I x K
PID_all = X.PID_v_plugin;

% x-axis categories (bins = 2:9) and that trial slice
x_bins = 2:9;
PID_slice = squeeze(PID_all(trialIdx,:,:,:)); % C x I x K
PID_slice = permute(PID_slice,[3 1 2]);       % K x C x I

% GT from GroundTruth_value cell via vertcat (like your code)
Lgt = load(F2_gt,'GroundTruth_value');
PID_gt_all = Lgt.GroundTruth_value;           % cell, length K (bins)
GT_all = vertcat(PID_gt_all{:});              % K x 5
GT_joint = sum(GT_all(:,1:4),2); GT_syn = GT_all(:,4); GT_red = GT_all(:,1);

% Means/3*SEM across iters → bias per K
[m_joint2, s_joint2] = m_sem3_rows(PID_slice(:,1,:));  % K x 1
[m_syn2,   s_syn2  ] = m_sem3_rows(PID_slice(:,2,:));
[m_red2,   s_red2  ] = m_sem3_rows(PID_slice(:,3,:));
plot_fill(ax2, x_bins, m_joint2 - GT_joint, s_joint2, cols.Joint, lw, aFill);
plot_fill(ax2, x_bins, m_syn2   - GT_syn,   s_syn2,   cols.Syn,   lw, aFill);
plot_fill(ax2, x_bins, m_red2   - GT_red,   s_red2,   cols.Red,   lw, aFill);
% yline(ax2,0,'k-','LineWidth',lw,'HandleVisibility','off');
set(ax2,'FontSize',FontSize,'LineWidth',lw);
xlabel(ax2,'Number of bins (R)');% yticks(ax2,[]); 
title(ax2,'Bins (N=512)');
xlim([min(x_bins), max(x_bins)])
ylim([0, 0.71])
% ---------------- Panel 3: Alpha (trials=512, bins=4) ----------------
ax3 = nexttile; hold(ax3,'on');
% Per-alpha files in Broja: Simuldata_bit_of_all_<alpha>_I_BROJA.mat
% Each is trials x comps x iters. We pick the trial closest to 512.
alpha_list = [1:15 20:5:60];
m_joint3 = nan(numel(alpha_list),1);
m_syn3   = nan(numel(alpha_list),1);
m_red3   = nan(numel(alpha_list),1);
s_joint3 = m_joint3; s_syn3 = m_syn3; s_red3 = m_red3;
GT_joint3 = m_joint3; GT_syn3 = m_syn3; GT_red3 = m_red3;

for i = 1:numel(alpha_list)
    a = alpha_list(i);
    Fa_data = sprintf('Results/Broja/Simuldata_bit_of_all_%d_I_BROJA.mat', a);
    Fa_gt   = sprintf('Results/Broja/GroundTruth_bit_of_all_%d.mat', a);
    assert_exist(Fa_data); assert_exist(Fa_gt);

    S = load(Fa_data,'PID_v_plugin');              % T x C x I
    P = S.PID_v_plugin;
    t_sel = nearest_trial_idx(size(P,1), alphatrials);     % choose trial closest to 512
    Xc = squeeze(P(t_sel,:,:));                    % C x I

    % GT like bin code: vertcat across cell, but per-alpha GT can be numeric
    G = load(Fa_gt,'GroundTruth_value');
    GTa = pick_GT_likeBinCode(G.GroundTruth_value, 4, a);  % 1x5 (bins=4, alpha=a)
    GT_joint3(i) = sum(GTa(1:4));%GTa(6); 
    GT_syn3(i)   = GTa(4); 
    GT_red3(i)   = GTa(1);
    % if a == 60
    %     disp('hola')
    % end

    % mean/3*SEM across iters (columns)
    [m_joint3(i), s_joint3(i)] = m_sem3_col(sum(Xc(2:5,:),1));
    [m_syn3(i),   s_syn3(i)]   = m_sem3_col(Xc(2,:));
    [m_red3(i),   s_red3(i)]   = m_sem3_col(Xc(3,:));
end

plot_fill(ax3, GT_joint3, m_joint3 - GT_joint3, s_joint3, cols.Joint, lw, aFill);
plot_fill(ax3, GT_joint3, m_syn3   - GT_syn3,   s_syn3,   cols.Syn,   lw, aFill);
plot_fill(ax3, GT_joint3, m_red3   - GT_red3,   s_red3,   cols.Red,   lw, aFill);



% yline(ax3,0,'k-','LineWidth',lw,'HandleVisibility','off');
set(ax3,'FontSize',FontSize,'LineWidth',lw);
xlabel(ax3,'Information level [bits]'); %yticks(ax3,[]);
title(ax3,['\alpha ' sprintf('(N=%d, R=4)',4*alphatrials)] );
xlim([0, 1.1])
ylim([0, 0.71])
% Link y across the row for comparability
% linkaxes([ax1 ax2 ax3],'y');
set(gcf,'Renderer','painters');

% --- optional saves ---
%saveas(gcf,'Figures_mat/discrete_panels_binsCodeStyle.png');
% print(gcf,'-dsvg','discrete_panels_fig2', 'ContentType','vector');
outDir = 'Figures_mat';
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

exportgraphics(gcf, fullfile(outDir, 'Figure_2_discrete.svg'), 'ContentType','vector');



end % function


% ===================== Helpers (same spirit as your bin code) =====================

function assert_exist(f)
    if ~exist(f,'file'), error('File not found: %s', f); end
end

function tvals = map_trials_axis(T)
    % Map trial index → values like your conventions
    base = [16 32 64 128 256 512 1024 2048 4096];
    tvals = base(1:min(T,numel(base)));
end

function idx = nearest_trial_idx(T, target)
    vals = map_trials_axis(T);
    [~,idx] = min(abs(vals - target));
end

function v = pick_GT_likeBinCode(GTv, n_bins, alpha_val)
% Reproduce your bin GT handling for Broja files:
% - If cell: try structs with fields (bins/alpha/value); fallback to first numeric-like.
% - Return a row vector (1x5) so indexing [1:4] works.
    if iscell(GTv)
        % Try to find a struct with matching bins/alpha
        for k=1:numel(GTv)
            x = GTv{k};
            if isstruct(x) && isfield(x,'value')
                okb = ~isfield(x,'bins')  || x.bins  == n_bins;
                oka = ~isfield(x,'alpha') || x.alpha == alpha_val;
                if okb && oka
                    v = row5(x.value); return;
                end
            elseif isnumeric(x) && numel(x)>=4
                v = row5(x); return;
            end
        end
        % fallback: last numeric entry
        for k=numel(GTv):-1:1
            x = GTv{k};
            if isnumeric(x) && numel(x)>=4
                v = row5(x); return;
            end
        end
        error('Could not parse GroundTruth_value cell.');
    else
        v = row5(GTv);
    end
end

function r = row5(x)
    x = x(:).';
    if numel(x) < 6, x(6) = NaN; end   % pad if only 4 given
    r = x(1:6);
end

function [m, s3] = m_sem3(XTxi)
% XTxi: T x I   (rows: x-axis points; columns: iter)
    m  = mean(XTxi, 2, 'omitnan');
    ni = max(1,size(XTxi,2));
    s3 = 3 * std(XTxi, 0, 2, 'omitnan') / sqrt(ni);
end

function [m, s3] = m_sem3_rows(Kx1xI)
% K x 1 x I -> treat as K x I
    X = squeeze(Kx1xI);
    if isvector(X), X = X(:); end
    m  = mean(X, 2, 'omitnan');
    ni = max(1,size(X,2));
    s3 = 3 * std(X, 0, 2, 'omitnan') / sqrt(ni);
end

function [m, s3] = m_sem3_col(vi)
% 1 x I vector of iters
    v = vi(:);
    m  = mean(v,'omitnan');
    ni = max(1,numel(v));
    s3 = 3 * std(v,0,'omitnan') / sqrt(ni);
end

function plot_fill(ax, x, m, s3, col, lw, a)
% Robust shaded band (sorted x, deduped), then mean line
    x = x(:).'; m = m(:).'; s = s3(:).';
    mask = isfinite(x) & isfinite(m) & isfinite(s);
    x = x(mask); m = m(mask); s = s(mask);
    if numel(x)<2, plot(ax,x,m,'LineWidth',lw,'Color',col); return; end
    [x,ord] = sort(x); m = m(ord); s = s(ord);
    [x,iu]  = unique(x,'stable'); m = m(iu); s = s(iu);
    xx = [x, fliplr(x)];
    yy = [m+s, fliplr(m-s)];
    fill(ax, xx, yy, col, 'FaceAlpha', a, 'EdgeColor','none','HandleVisibility','off');
    plot(ax, x, m, 'LineWidth', lw, 'Color', col);
end
