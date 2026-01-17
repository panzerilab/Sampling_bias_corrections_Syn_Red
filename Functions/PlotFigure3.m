function PlotFigure3()
% ---------------- Visuals ----------------
% Turn off LaTeX everywhere (tex = default math-ish, none = literal)
set(groot,'defaultTextInterpreter','tex');
set(groot,'defaultAxesTickLabelInterpreter','tex');
set(groot,'defaultLegendInterpreter','tex');

% Or if you want literal (no TeX parsing), use 'none' for all three lines above.
%%
% Force a common font that Inkscape will have
set(groot,'defaultAxesFontName','Arial');
set(groot,'defaultTextFontName','Arial');
set(groot,'defaultLegendFontName','Arial');

set(0,'DefaultTextFontName','Arial'); 
set(0,'DefaultAxesFontName','Arial'); 
set(0,'DefaultLegendFontName','Arial');
FontSize = 9; lw = 1;
cols = struct('Joint', [0.9290,0.6940,0.1250], ... % yellow
              'Ind',   [0.4, 0.26, 0.13], ... 
              'Union', [0.8500, 0.3250, 0.0980], ...
              'Single',[0.4940, 0.1840, 0.5560]);
aFill = 0.18;

fig = figure('Units','centimeters','Position',[1 1 16.5 10]);
tl  = tiledlayout(fig,2,3,'TileSpacing','compact','Padding','compact');

%% I_MMI
ax1 = nexttile(1); hold(ax1,'on');
% Files (Broja/high), plugin arrays: [trials x comps x iters]
F_data = 'Results/IMMI/Simuldata_bit_of_all_high_I_MMI.mat';
F_gt   = 'Results/IMMI/GroundTruth_bit_of_all_high.mat';
assert_exist(F_data); assert_exist(F_gt);

% load plugin
L = load(F_data,'PID_v_plugin'); 
PID_v = L.PID_v_plugin;                 % T x C x I
T = size(PID_v,1);
n_iter = size(PID_v,3);
comp = @(c) squeeze(PID_v(:,c,:));      % T x I

% GT matching bins=4, alpha=3 using same approach as your bin plot
G = load(F_gt,'GroundTruth_value');
GTv = G.GroundTruth_value;            % cell or numeric
GT_vec = pick_GT_likeBinCode(GTv, 4, 3);  % 1x5
GT_joint = sum(GT_vec(1:4)); 
GT_union = sum(GT_vec(1:3)); 
GT_ind = GT_vec(5);
GT = GT_vec;
Joint_GT = sum(GT(1:4));
Union_GT = sum(GT(2:3))+GT(1);

Info_singNeuro_GT = max(GT(2:3))+GT(1);
Info_singNeuro = max(PID_v(:,4:5,:), [], 2)+(PID_v(:,3,:));
Info_singNeuro_corr = Info_singNeuro - Info_singNeuro_GT;

JointCorr = PID_v(:,1,:)-Joint_GT;
UnionCorr = PID_v(:,6,:)-Union_GT;

mean_singNeuro = mean(Info_singNeuro_corr , 3, 'omitnan');
mean_UNION = mean(UnionCorr , 3, 'omitnan');
mean_Joint = mean(JointCorr , 3, 'omitnan');

SEM_singNeuro = 3 * std(Info_singNeuro_corr,1,3, 'omitnan') / sqrt(n_iter);
SEM_Joint = 3 * std(JointCorr,1,3, 'omitnan') / sqrt(n_iter);
SEM_UNION = 3 * std(UnionCorr,1,3, 'omitnan') / sqrt(n_iter);

x_trials = 4* map_trials_axis(T);            % [16 32 64 128 256 512 ...]
plot_fill(ax1, x_trials, mean_Joint, SEM_Joint, cols.Joint, lw, aFill);
plot_fill(ax1, x_trials, mean_UNION, SEM_UNION,   cols.Union, lw, aFill);
plot_fill(ax1, x_trials, mean_singNeuro, SEM_singNeuro,   cols.Ind, lw, aFill);
set(ax1,'XScale','log','FontSize',FontSize,'LineWidth',lw); 
xlabel(ax1,'Trials (N)'); ylabel(ax1,'Bias [bits]');
% yline(ax1,0,'k-','LineWidth',lw,'HandleVisibility','off');
title(ax1,'I_{MMI}');
legend({'Joint','Union','Single'},'Box','off','Location','best');
ylim([0,.5])

ax2 = nexttile(4); hold(ax2,'on');
% EXACTLY mirror your bin-plot handling: Bin_sweep/low; trialIdx mapping
ResultFolder = 'Results/Bin_sweep';
F2_data = fullfile(ResultFolder,'Simuldata_bit_of_all_low_I_MMI.mat');
F2_gt   = fullfile(ResultFolder,'GroundTruth_bit_of_all_low_IMMI.mat');
assert_exist(F2_data); assert_exist(F2_gt);

% Load: trials x comps x iters x cats  -> we want a chosen trial slice
% trialIdx = 2;          % your original setting (maps to N=512)
% X = load(F2_data, 'PID_v_plugin');           % T x C x I x K
% PID_all = X.PID_v_plugin;
% 
% % x-axis categories (bins = 2:9) and that trial slice
% x_bins = 2:9;
% PID_slice = squeeze(PID_all(trialIdx,:,:,:)); % C x I x K
% PID_slice = permute(PID_slice,[3 1 2]);       % K x C x I
% 
% % GT from GroundTruth_value cell via vertcat (like your code)
% Lgt = load(F2_gt,'GroundTruth_value');
% PID_gt_all = Lgt.GroundTruth_value;           % cell, length K (bins)
% GT_all = vertcat(PID_gt_all{:});              % K x 5
% GT_joint = sum(GT_all(:,1:4),2); GT_syn = GT_all(:,4); GT_red = GT_all(:,1);
trialIdx = 2;

PID_v = load(F2_data, 'PID_v_plugin');   
PID_v = PID_v.PID_v_plugin;
PID_v = squeeze(PID_v(trialIdx,:, :, :));

GroundTruth = load(F2_gt,'GroundTruth_value');
GroundTruth = GroundTruth.GroundTruth_value;

GT_Joint = zeros(1,length(GroundTruth));
GT_Union = zeros(1,length(GroundTruth));
GT_singleNeuron = zeros(1,length(GroundTruth));

for idxGt = 1:length(GroundTruth)
    GT_tmp = GroundTruth{idxGt};
    GT_Joint(idxGt) = sum(GT_tmp(1:4));
    GT_Union(idxGt) = max(GT_tmp(2:3))+GT_tmp(1);
    GT_singleNeuron(idxGt) = max(GT_tmp(2:3))+GT_tmp(1);
end

Inf_singleNeuron = max(PID_v(4:5,:,:), [], 1)+(PID_v(3,:,:));

mean_singleNeuron = squeeze(mean(Inf_singleNeuron,2, 'omitnan'))  - GT_singleNeuron';
mean_Union = squeeze(mean(PID_v(6,:,:),2, 'omitnan'))     -GT_Union';
mean_Joint = squeeze(mean(PID_v(1,:,:),2, 'omitnan'))     -GT_Joint';

SEM_singleNeuron = squeeze(3 * std(Inf_singleNeuron,1,2,'omitnan') / sqrt(n_iter));
SEM_Joint = squeeze(3 * std(PID_v(1,:,:),1,2,'omitnan') / sqrt(n_iter));
SEM_Union = squeeze(3 * std(PID_v(6,:,:),1,2,'omitnan') / sqrt(n_iter));

x_bins = 2:9;
% Means/3*SEM across iters → bias per K
% [m_joint2, s_joint2] = m_sem3_rows(PID_slice(:,1,:));  % K x 1
% [m_syn2,   s_syn2  ] = m_sem3_rows(PID_slice(:,2,:));
% [m_red2,   s_red2  ] = m_sem3_rows(PID_slice(:,3,:));
plot_fill(ax2, x_bins, mean_Joint, SEM_Joint, cols.Joint, lw, aFill);
plot_fill(ax2, x_bins, mean_Union,   SEM_Union,   cols.Union,   lw, aFill);
plot_fill(ax2, x_bins, mean_singleNeuron,   SEM_singleNeuron,   cols.Ind,   lw, aFill);
% yline(ax2,0,'k-','LineWidth',lw,'HandleVisibility','off');
set(ax2,'FontSize',FontSize,'LineWidth',lw);
xlabel(ax2,'Number of bins (R)');% yticks(ax2,[]); 
ylabel(ax2,'Bias [bits]')
xlim([min(x_bins), max(x_bins)])


%% I_BROJA
ax1 = nexttile(2); hold(ax1,'on');
% Files (Broja/high), plugin arrays: [trials x comps x iters]
F_data = 'Results/Broja/Simuldata_bit_of_all_high_I_BROJA.mat';
F_gt   = 'Results/Broja/GroundTruth_bit_of_all_high.mat';
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
GT_joint = sum(GT_vec(1:4)); 
GT_union = sum(GT_vec(1:3)); 
GT_ind = GT_vec(5);

% mean±3*SEM (SEM over iters), Bias = mean(meas) - GT
[m_joint, s_joint] = m_sem3(comp(1));   % Joint measured vs trials
[m_union, s_union] = m_sem3(comp(8));
[m_ind,     s_ind] = m_sem3(comp(7));
x_trials = 4* map_trials_axis(T);            % [16 32 64 128 256 512 ...]
plot_fill(ax1, x_trials, m_joint - GT_joint, s_joint, cols.Joint, lw, aFill);
plot_fill(ax1, x_trials, m_union - GT_union, s_union,   cols.Union, lw, aFill);
plot_fill(ax1, x_trials, m_ind   - GT_ind,     s_ind,   cols.Ind, lw, aFill);
set(ax1,'XScale','log','FontSize',FontSize,'LineWidth',lw);
xlabel(ax1,'Trials (N)'); 
% ylabel(ax1,'Bias [bits]');
% yline(ax1,0,'k-','LineWidth',lw,'HandleVisibility','off');
title(ax1,'BROJA');
legend({'Joint','Union','I_{Ind}'},'Box','off','Location','best');
ylim([0,.5])

ax2 = nexttile(5); hold(ax2,'on');
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

PID_v = load(F2_data, 'PID_v_plugin'); 
PID_v = PID_v.PID_v_plugin;
PID_v = squeeze(PID_v(trialIdx,:, :, :));

GroundTruth = load(F2_gt, 'GroundTruth_value');
GroundTruth = GroundTruth.GroundTruth_value;

GT_Joint = zeros(1,length(GroundTruth));
GT_Union = zeros(1,length(GroundTruth));
GT_Union_IMMI = zeros(1,length(GroundTruth));
GT_IInd = zeros(1,length(GroundTruth));

for idxGt = 1:length(GroundTruth)
    GT_tmp = GroundTruth{idxGt};
    GT_Joint(idxGt) = sum(GT_tmp(1:4));
    GT_Union(idxGt) = sum(GT_tmp(2:3))+GT_tmp(1);
    GT_Union_IMMI(idxGt) = max(GT_tmp(2:3))+GT_tmp(1);
    GT_IInd(idxGt) = GT_tmp(5);
end

mean_Union = squeeze(mean(PID_v(6,:,:),2, 'omitnan'))     -GT_Union';
mean_Joint = squeeze(mean(PID_v(1,:,:),2, 'omitnan'))     -GT_Joint';
mean_IInd  = squeeze(mean(PID_v(7,:,:),2, 'omitnan'))     -GT_IInd';

SEM_Joint = squeeze(2 * std(PID_v(1,:,:),1,2,'omitnan') / sqrt(n_iter));
SEM_Union = squeeze(2 * std(PID_v(6,:,:),1,2,'omitnan') / sqrt(n_iter));
SEM_IInd  = squeeze(2 * std(PID_v(7,:,:),1,2,'omitnan') / sqrt(n_iter));
plot_fill(ax2, x_bins, mean_Joint, SEM_Joint, cols.Joint, lw, aFill);
plot_fill(ax2, x_bins, mean_Union,   SEM_Union,   cols.Union,   lw, aFill);
plot_fill(ax2, x_bins, mean_IInd,   SEM_IInd,   cols.Ind,   lw, aFill);
% yline(ax2,0,'k-','LineWidth',lw,'HandleVisibility','off');
set(ax2,'FontSize',FontSize,'LineWidth',lw);
xlabel(ax2,'Number of bins (R)');% yticks(ax2,[]); 
% ylabel(ax2,'Bias [bits]')

xlim([min(x_bins), max(x_bins)])








%% I_gPID
filename = 'Results/Gauss/Finalresults_across_M_and_ntrials.mat';
data = load(filename);
GroundTruth =  data.GT_results;
PID_v =  data.sampled_results;
M_vals =  data.M_vals;
n_iter = size(PID_v,6);
dim = 20;
thumbRule = 4*dim;
dimidx = find(M_vals==dim);
info_amount = 3;
PID_v = squeeze(PID_v(:, :, dimidx, 1, info_amount, :));
GroundTruth = squeeze(GroundTruth(:, :, dimidx,1, info_amount));

GT_Union = GroundTruth(2,:,:)+GroundTruth(5,:,:);
GT_Joint = GroundTruth(1,:,:);
GT_PInd = GroundTruth(10,:,:);

Union = PID_v(2,:,:)+PID_v(5,:,:);
mean_UNION =  mean(Union ,3, 'omitnan')-GT_Union;
mean_Joint =  mean(PID_v(1,:,:),3, 'omitnan')-GT_Joint;
mean_PInd =   mean(PID_v(10,:,:),3, 'omitnan')-GT_PInd;

SEM_Joint = 3 * std(PID_v(1,:,:),1,3,'omitnan') / sqrt(n_iter);
SEM_UNION = 3 * std(Union,1,3,'omitnan') / sqrt(n_iter);
SEM_PInd =  3 * std(PID_v(10,:,:),1,3,'omitnan') / sqrt(n_iter);

ax1 = nexttile(3); hold(ax1,'on');
% % Files (Broja/high), plugin arrays: [trials x comps x iters]
% F_data = 'Results/Gauss/Finalresults_across_M_and_ntrials.mat';
% assert_exist(F_data); 
% 
% % load plugin
% L = load(F_data,'sampled_results'); 
% PID = L.sampled_results;                 % T x C x I
% T = size(PID,1);
% comp = @(c) squeeze(PID(:,c,:));      % T x I
% 
% % GT matching bins=4, alpha=3 using same approach as your bin plot
% G = load(F_data,'GT_results');
% GTv = G.GT_results;            % cell or numeric
% GT_vec = pick_GT_likeBinCode(GTv, 4, 3);  % 1x5
% GT_joint = sum(GT_vec(1:4)); 
% GT_union = sum(GT_vec(1:3)); 
% GT_ind = GT_vec(5);

% mean±3*SEM (SEM over iters), Bias = mean(meas) - GT
% [m_joint, s_joint] = m_sem3(comp(1));   % Joint measured vs trials
% [m_union, s_union] = m_sem3(comp(8));
% [m_ind,     s_ind] = m_sem3(comp(7));
x_trials = data.ntrials_vals;            % [16 32 64 128 256 512 ...]
plot_fill(ax1, x_trials, mean_Joint, SEM_Joint, cols.Joint, lw, aFill);
plot_fill(ax1, x_trials, mean_UNION, SEM_UNION, cols.Union, lw, aFill);
plot_fill(ax1, x_trials, mean_PInd,  SEM_PInd,  cols.Ind, lw, aFill);
set(ax1,'XScale','log','FontSize',FontSize,'LineWidth',lw);
xlabel(ax1,'Trials (N)'); 
% ylabel(ax1,'Bias [bits]'); 
% yline(ax1,0,'k-','LineWidth',lw,'HandleVisibility','off');
xticks([100, 1000])
title(ax1,'gPID');
legend({'Joint','Union','I_{Ind}'},'Box','off','Location','best');
% ylim([0,.4])

ax2 = nexttile(6); hold(ax2,'on');

data = load(filename);
GroundTruth =  data.GT_results;
PID_v =  data.sampled_results;
trial_vals =  data.ntrials_vals;

trialidx = find(trial_vals==128);

PID_v = squeeze(PID_v(:,trialidx,:,1,info_amount,:));
GroundTruth = squeeze(GroundTruth(:,trialidx,:,1,info_amount));

thumbRule = floor(trialIdx/4);

% PID_v = squeeze(PID_v(:, 1:length(trial_vals), 1, info_amount, :));
% GroundTruth = squeeze(GroundTruth(:, 1:length(trial_vals), 1, info_amount));

GT_Union = GroundTruth(2,:,:)+GroundTruth(5,:,:);
GT_Joint = GroundTruth(1,:,:);
GT_PInd = GroundTruth(10,:,:);

Union = PID_v(2,:,:)+PID_v(5,:,:);
mean_UNION =  mean(Union ,3, 'omitnan')-GT_Union;
mean_Joint =  mean(PID_v(1,:,:),3, 'omitnan')-GT_Joint;
mean_PInd =   mean(PID_v(10,:,:),3, 'omitnan')-GT_PInd;

SEM_Joint = 3 * std(PID_v(1,:,:),1,3,'omitnan') / sqrt(n_iter);
SEM_UNION = 3 * std(Union,1,3,'omitnan') / sqrt(n_iter);
SEM_PInd =  3 * std(PID_v(10,:,:),1,3,'omitnan') / sqrt(n_iter);


plot_fill(ax2, M_vals, mean_Joint, SEM_Joint, cols.Joint, lw, aFill);
plot_fill(ax2, M_vals, mean_UNION,   SEM_UNION,   cols.Union,   lw, aFill);
plot_fill(ax2, M_vals, mean_PInd,   SEM_PInd,   cols.Ind,   lw, aFill);
% yline(ax2,0,'k-','LineWidth',lw,'HandleVisibility','off');
set(ax2,'FontSize',FontSize,'LineWidth',lw);
xlabel(ax2,'Number of dim (d)');% yticks(ax2,[]); 
% ylabel(ax2,'Bias [bits]');
xlim([min(M_vals), max(M_vals)])

set(gcf, 'Renderer', 'painters');

outDir = 'Figures_mat';
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

exportgraphics(gcf, fullfile(outDir, 'Figure_3.svg'), 'ContentType','vector');
end

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
    if numel(x) < 5, x(5) = NaN; end   % pad if only 4 given
    r = x(1:5);
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
    m  = median(v,'omitnan');
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
