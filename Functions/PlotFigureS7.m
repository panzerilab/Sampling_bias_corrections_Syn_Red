function PlotFigureS7(datatype, varargin)
% plotAcrossInfoLevels_alphaX
% Plots bias vs information level (alpha on x-axis) with fixed settings:
%   DISCRETE: bins = 4, trials = 128
%   GAUSSIAN: M = 20, trials = 512
%
% Panels: Joint | Syn | Red
% Lines: one per bias method
%   - DISCRETE  : plugin, qe, shuff, qeShuff, weighted (if weight matrix available), venkatesh
%   - GAUSSIAN  : plugin, resample, shuffle, shuff-resamp, venkatesh
%
% Usage:
%   plotAcrossInfoLevels_alphaX('discr', 'Alphas', 1:15, 'ResultFolder','Broja');
%   plotAcrossInfoLevels_alphaX('gauss', 'Alphas', [1 1.1 1.2 1.5 2 10000], 'ResultFolder','Gauss');
%
% Name-Value options:
%   'Alphas'        : vector of alphas to sweep (default: 1:15 for DISCRETE;
%                                         [1 1.1 1.2 1.5 2 10000] for GAUSSIAN)
%   'ResultFolder'  : base results folder (default: 'Broja' for DISCRETE; 'Gauss' for GAUSSIAN)
%   'SavePrefix'    : output file prefix (default auto)
%   'LineAlpha'     : CI face alpha (default 0.18)
%   'Colors'        : Nx3 colormap for series (optional; default = MATLAB lines)
%
% Output files:
%   <SavePrefix>.png (300 dpi), <SavePrefix>.svg

% ---------- Parse & defaults ----------
if nargin < 1 || isempty(datatype), datatype = 'discr'; end
datatype = char(lower(string(datatype)));
% --- fixed color order by bias method ---
% cmap = [
%     0.9290, 0.6940, 0.1250;  % 1: plugin
%     0.0000, 0.4470, 0.7410;  % 2: qe (discr) / resample (gauss)
%     0.8500, 0.3250, 0.0980;  % 3: shuff-sub / shuffle
%     0.4940, 0.1840, 0.5560;  % 4: merged (qeShuff / shuff-resamp / weighted)
%     0.4660, 0.6740, 0.1880;  % 5: venkatesh
%     0.5,    0.5,    0.5;     % 6: (spare)
%     0.4,    0.26,   0.13;    % 7: (spare)
%     0.4940, 0.1840, 0.5560]; % 8: use this colors
cmap = [
    0.9290, 0.6940, 0.1250;  % 1: Joint
    0,      0.4470, 0.7410;  % 2: Red
    0.8500, 0.3250, 0.0980;  % 3: Union
    0.4940, 0.1840, 0.5560;  % 4: SR
    0.4660, 0.6740, 0.1880;  % 5: Syn
    0.5,    0.5,    0.5;     % 6: U1+U2
    0.4,    0.26,   0.13;    % 7: PInd
    0.4940, 0.1840, 0.5560]; % 8: PInd (alt)

p = inputParser;
addParameter(p, 'Alphas', [], @(x)isnumeric(x)&&~isempty(x));
addParameter(p, 'ResultFolder', '', @(s)ischar(s)||isstring(s));
addParameter(p, 'SavePrefix', '', @(s)ischar(s)||isstring(s));
addParameter(p, 'LineAlpha', 0.18, @(x)isnumeric(x)&&isscalar(x)&&x>=0&&x<=1);
addParameter(p, 'Colors', [], @(x)(isnumeric(x)&&size(x,2)==3)||isempty(x));
parse(p, varargin{:});

alphas       = p.Results.Alphas;
ResultFolder = char(string(p.Results.ResultFolder));
save_prefix  = char(string(p.Results.SavePrefix));
aFill        = p.Results.LineAlpha;
cmap_in      =cmap;% p.Results.Colors;

if strcmp(datatype,'discr')
    if isempty(alphas), alphas = 1:15; end
    if isempty(ResultFolder), ResultFolder = 'Broja'; end
    if isempty(save_prefix),  save_prefix  = 'acrossInfo_discrete_alphaX'; end
else
    if isempty(alphas), alphas = [1 1.1 1.2 1.5 2 10000]; end
    if isempty(ResultFolder), ResultFolder = 'Gauss'; end
    if isempty(save_prefix),  save_prefix  = 'acrossInfo_gaussian_alphaX'; end
end

% ---------- Visual defaults ----------
set(0,'DefaultTextFontName','Arial'); 
set(0,'DefaultAxesFontName','Arial'); 
set(0,'DefaultLegendFontName','Arial');
FontSize = 8; lw = 1;

comp_names = {'Joint','Syn','Red', 'Unq'};
% series (bias methods) color order
cmap = lines(8);
if ~isempty(cmap_in), cmap = cmap_in; end

% ---------- Assemble series (each bias method = one line vs alpha) ----------
series = struct('name',{},'color',{},'mean',{},'sem',{});
if strcmp(datatype,'discr')
    % ----- DISCRETE (bins=4, trials=128) -----
    bins_fixed   = 4;
    trials_fixed = 128;

    % Candidate bias list (keep order = legend order)
    bias_list = {'qe','shuff','weighted'}; %'qeshuff'}; %'plugin',

    ser_ix = 0;
    for bi = 1:numel(bias_list)
        bname = bias_list{bi};
        try
            [mJ, sJ, mS, sS, mR, sR, gtJ] = build_discrete_alpha_line(ResultFolder, alphas, bins_fixed, trials_fixed, bname);
            ser_ix = ser_ix + 1;
            series(ser_ix).name  = canonical_bias_name(bname);
            series(ser_ix).color = color_for_bias(datatype, bname, cmap);
            series(ser_ix).mean  = [mJ(:), mS(:), mR(:)];  % nAlpha x 3
            series(ser_ix).sem   = [sJ(:), sS(:), sR(:)];  % nAlpha x 3
        catch
            % Skip silently if unavailable (e.g., no weighted matrix present)
        end
    end
    x_vals  = gtJ(:)'; %alphas(:)'; 
    x_label = 'Information level [bits]';

else
    % ----- GAUSSIAN (M=20, trials=512) -----
    M_fixed      = 20;
    trials_fixed = 128;

    [Samp, GT, ntrials_vals, M_vals] = load_gauss_files(ResultFolder);
    [~, Midx] = min(abs(M_vals - M_fixed));
    [~, tidx] = min(abs(ntrials_vals - trials_fixed));

    % Bias list
    bias_list = {'resample','shuffle','shuff-resamp','venkatesh'}; %'plugin',
    ser_ix = 0;

    for bi = 1:numel(bias_list)
        bname = bias_list{bi};
        try
            [mJ, sJ, mS, sS, mR, sR, mU, sU,gtJ] = build_gauss_alpha_line(Samp, GT, tidx, Midx, alphas, bname);
            ser_ix = ser_ix + 1;
            series(ser_ix).name  = canonical_bias_name(bname);
            series(ser_ix).color = color_for_bias(datatype, bname, cmap);
            series(ser_ix).mean  = [mJ(:), mS(:), mR(:), mU(:)];  % nAlpha x 3
            series(ser_ix).sem   = [sJ(:), sS(:), sR(:), sU(:)];  % nAlpha x 3
        catch
            disp('hola')% Skip missing bias gracefully
        end
    end
    x_vals  = gtJ(:)'; %1./alphas(:)'; 
    x_label = 'Information level [bits]';
end

if isempty(series)
    error('No series could be assembled for %s backend with the requested settings.', datatype);
end

% ---------- PLOTTING ----------
fig = figure('Units','centimeters','Position',[1 1 18 6]);
tl  = tiledlayout(fig,1,3,'TileSpacing','compact','Padding','compact');
axs = gobjects(1,3);  % collect handles

for c = 1:3
    ax = nexttile; axs(c) = ax; hold(ax,'on');
    yline(ax, 0, 'k-', 'LineWidth', lw, 'HandleVisibility','off');

    for k = 1:numel(series)
        m = series(k).mean(:,c);
        s = series(k).sem(:,c);
        xv = [x_vals(:); flipud(x_vals(:))];
        yv = [m + s; flipud(m - s)];
        fill(ax, xv, yv, series(k).color, 'FaceAlpha', aFill, 'EdgeColor','none','HandleVisibility','off');
        plot(ax, x_vals, m, 'LineWidth', lw, 'Color', series(k).color, 'DisplayName', series(k).name);
    end

    set(ax,'LineWidth',lw,'FontSize',FontSize);
    xlabel(ax, x_label);
    if c==1, ylabel(ax,'Bias [bits]');end %
    title(ax, comp_names{c}, 'FontSize', FontSize, 'FontWeight','bold');
    xlim(ax, [min(x_vals) max(x_vals)]);
    if strcmp(datatype,'discr')
        ylim(ax, [-.1, .1]);
        xlim(ax, [min(x_vals) 1.1]);
    end
    
end
linkaxes(axs, 'y');
% legend once (left-most panel)
axes(tl.Children(end));
lgd = legend('Location','northwest','Box','off'); 
lgd.ItemTokenSize = [12 6];
lgd.FontSize = 8;

set(gcf,'Renderer','painters');
print(gcf,'-dsvg',[save_prefix '.svg']);
end

% ====================== DISCRETE builder ======================
function [mJ, sJ, mS, sS, mR, sR, gtJ] = build_discrete_alpha_line(ResultFolder, alphas, bins_fixed, trials_fixed, bias_name)
% One line vs alpha for fixed bins & trials; bias_name ∈ {'plugin','qe','shuff','qeShuff','weighted','venkatesh'}
nA = numel(alphas);
mJ = nan(nA,1); sJ = mJ; mS = mJ; sS = mJ; mR = mJ; sR = mJ;
gtJ = [];

for i = 1:nA
    a = alphas(i);
    Fd = fullfile('Results', ResultFolder, sprintf('Simuldata_bit_of_all_%g_I_BROJA.mat', a));
    Fg = fullfile('Results', ResultFolder, sprintf('GroundTruth_bit_of_all_%g.mat', a));
    assert_exist(Fd); assert_exist(Fg);

    % Load PID arrays (trials x comps x iters); pick nearest trial to trials_fixed
    Lp = load(Fd, 'PID_v_plugin');
    Pp = Lp.PID_v_plugin;
    tsel = nearest_trial_idx(size(Pp,1), trials_fixed);

    % Try load other bias variants if requested
    switch lower(bias_name)
        case 'plugin'
            X = squeeze(Pp(tsel,:,:));          % comps x iters

        case 'qe'
            L = load(Fd, 'PID_v_qe'); X = squeeze(L.PID_v_qe(tsel,:,:));

        case 'shuff'
            L = load(Fd, 'PID_v_shuff'); X = squeeze(L.PID_v_shuff(tsel,:,:));

        case 'qeshuff'
            % if exist(Fd,'file')
            %     if ~exist_var_in(Fd,'PID_v_qeShuff')
            %         error('PID_v_qeShuff not found in %s', Fd);
            %     end
            % end
            % L = load(Fd, 'PID_v_qeShuff'); 
            L1 = load(Fd, 'PID_v_qe'); X1 = squeeze(L1.PID_v_qe(tsel,:,:));
            L2 = load(Fd, 'PID_v_shuff'); X2 = squeeze(L2.PID_v_shuff(tsel,:,:));
            X = (X1+X2)/2;

        case 'weighted'
            % build weighted = w*QE + (1-w)*SHUFF using adaptive weight matrix (if available)
            Lqe = load(Fd,'PID_v_qe');     Lsh = load(Fd,'PID_v_shuff');
            Xqe = squeeze(Lqe.PID_v_qe(tsel,:,:));      % comps x iters
            Xsh = squeeze(Lsh.PID_v_shuff(tsel,:,:));   % comps x iters
            % weight selection based on joint info fraction at bins=4
            try
                Sw = load('adaptive_weight_matrix_allcases.mat','weight_matrix','info_levels_fine');
            catch
                Sw = load('adaptive_weight_matrix_allcases_2.mat','weight_matrix','info_levels_fine');
            end
            % Estimate info fraction from averaged QE/SHUFF joint
            J_est = mean((Xqe(1,:)+Xsh(1,:))/2,'omitnan');
            max_info = log2(bins_fixed);
            frac = min(max(J_est/max_info,0),1);
            [~, idx] = min(abs(Sw.info_levels_fine(:)-frac));
            % weight_matrix dims may differ across your files; choose a central slice safely
            w = Sw.weight_matrix(min(bins_fixed-1, size(Sw.weight_matrix,1)), ...
                                 min(tsel, size(Sw.weight_matrix,2)), ...
                                 idx);
            if ~isfinite(w), w = 0.5; end
            X = w.*Xqe + (1-w).*Xsh;

        case 'venkatesh'
            % build syn/red from plugin using union clamp with qeShuff GT at same alpha
            Lm = load(Fd,'PID_v_qeShuff');
            P  = squeeze(Pp(tsel,:,:));         % comps x iters
            GT = squeeze(Lm.PID_v_qeShuff(end,:,:)); % comps x iters (use large-trial estimate as GT-ish)
            % Joint GT & plugin debias factor
            joint_gt = mean(GT(1,:), 'omitnan');
            joint_pl = mean(P(1,:),  'omitnan');
            debias   = joint_gt / max(joint_pl,eps);
            % Clamp union between bounds from GT (component 3..5 are Red,U1,U2 in your layout)
            % Here, we reconstruct Syn/Red similar to your earlier venkatesh helper.
            % For fixed bins=4, use average across iter:
            i1_gt = mean(GT(3,:), 'omitnan') + mean(GT(4,:), 'omitnan'); % Red + U1
            i2_gt = mean(GT(3,:), 'omitnan') + mean(GT(5,:), 'omitnan'); % Red + U2
            union_est = sum(P(3:5,:),1);            % per-iter
            union_est = union_est * debias;
            union_est = min(max(union_est, i1_gt), min(i1_gt+i2_gt, joint_gt));
            syn_line  = joint_gt - union_est;       % per-iter
            red_line  = i1_gt + i2_gt - union_est;  % per-iter
            X = zeros(size(P));                     % comps x iters
            X(1,:) = joint_gt;                      % keep joint at GT
            X(2,:) = syn_line;
            X(3,:) = red_line;
            X(4,:) = union_est - i2_gt;
            X(5,:) = union_est - i1_gt;

        otherwise
            error('Unknown discrete bias: %s', bias_name);
    end

    % Ground truth at bins=4 for this alpha
    G = load(Fg,'GroundTruth_value');
    GTa = pick_GT_likeBinCode(G.GroundTruth_value, bins_fixed, a); % 1x5
    GT_joint = sum(GTa(1:4)); GT_syn = GTa(4); GT_red = GTa(1);

    % Means & 3*SEM over iters → bias
    [mJ(i), sJ(i)] = m_sem3_col(X(1,:));   mJ(i) = mJ(i) - GT_joint;
    [mS(i), sS(i)] = m_sem3_col(X(2,:));   mS(i) = mS(i) - GT_syn;
    [mR(i), sR(i)] = m_sem3_col(X(3,:));   mR(i) = mR(i) - GT_red;

    gtJ = [gtJ GT_joint];
end
end

% ====================== GAUSSIAN builder ======================
function [mJ, sJ, mS, sS, mR, sR, mU, sU, gtJ] = build_gauss_alpha_line(Samp, GT, tidx, Midx, alphas, bias_name)
% Samp,GT dims: comps x nTrials x nM x nBias x nAlpha x nIters (GT's last dim optional)
% Returns line vs alpha (length = numel(alphas)) at fixed trial & M

if strcmp(bias_name, 'shuff-resamp')
    b1 = pick_bias_index_gauss('resample');
    b2 = pick_bias_index_gauss('shuffle');
else
    b = pick_bias_index_gauss(bias_name);
end


nAlphaAvail = size(Samp,5);
% Decide which alpha indices to sample:
a_is_indexish = all(abs(alphas - round(alphas))<1e-9 & alphas>=1 & alphas<=nAlphaAvail);
if a_is_indexish
    aidx = round(alphas(:))';
else
    % If true alpha values aren't provided in file, sample first N alphas
    aidx = 1:min(numel(alphas), nAlphaAvail);
end

nA = numel(aidx);
mJ = nan(nA,1); sJ = mJ; mS = mJ; sS = mJ; mR = mJ; sR = mR; 
mU = mJ; sU = mR;
gtJ = [];
for ii = 1:nA
    ai = aidx(ii);
    if strcmp(bias_name, 'shuff-resamp')
    P = squeeze((Samp(:, tidx, Midx, b1, ai, :)+Samp(:, tidx, Midx, b2, ai, :))/2);  % comps x iters
    else
        P = squeeze(Samp(:, tidx, Midx, b, ai, :));  % comps x iters
    end
    if isvector(P), P = P(:); end
    Gi= squeeze(GT(:,   tidx, Midx, 1, ai, :));  % comps x [iters]

    if ndims(Gi)==2, Gm = mean(Gi,2,'omitnan'); else, Gm = Gi; end
    Gm = Gm(:);

    % comps: Joint=1, Syn=6, Red=5
    [mJ(ii), sJ(ii)] = m_sem3_col(squeeze(P(1,:)));  mJ(ii) = mJ(ii) - Gm(1);
    [mS(ii), sS(ii)] = m_sem3_col(squeeze(P(6,:)));  mS(ii) = mS(ii) - Gm(6);
    [mR(ii), sR(ii)] = m_sem3_col(squeeze(P(5,:)));  mR(ii) = mR(ii) - Gm(5);
    [mU(ii), sU(ii)] = m_sem3_col(squeeze(P(2,:)));  mU(ii) = mU(ii) - Gm(2);
    gtJ = [gtJ Gm(1)];
end
end

% ====================== LOADERS & HELPERS ======================
function [Samp, GT, ntrials_vals, M_vals] = load_gauss_files(ResultFolder)
file_main = fullfile('Results', ResultFolder, 'Finalresults_across_M_and_ntrials.mat');
file_d80  = fullfile('Results', ResultFolder, 'Finalresults_across_M_and_ntrials_d80.mat');
if isfile(file_main)
    D = load(file_main);
elseif isfile(file_d80)
    D = load(file_d80);
else
    error('Could not find Gaussian results in Results/%s.', ResultFolder);
end
Samp = D.sampled_results;
GT   = D.GT_results;
ntrials_vals = D.ntrials_vals;
M_vals       = D.M_vals;
end

function assert_exist(f)
if ~exist(f,'file'), error('File not found: %s', f); end
end

function idx = nearest_trial_idx(T, target)
vals = map_trials_axis(T);
[~,idx] = min(abs(vals - target));
end

function tvals = map_trials_axis(T)
base = [16 32 64 128 256 512 1024 2048 4096];
tvals = base(1:min(T,numel(base)));
end

function name = canonical_bias_name(b)
b = string(b); %lower(string(b));
switch b
    case 'qeShuff', name = 'qeShuff';
    case 'qeshuff', name = 'qeShuff';
    case 'shuff-resamp', name = 'shuff-resamp';
    otherwise, name = char(b);
end
end

function [m, s3] = m_sem3_col(v)
v  = v(:);
m  = mean(v,'omitnan');
nI = max(1,numel(v));
s3 = 3 * std(v,0,'omitnan') / sqrt(nI);
end

function v = pick_GT_likeBinCode(GTv, n_bins, alpha_val)
% Robust GT extractor (works with cells or numeric)
if iscell(GTv)
    for k=1:numel(GTv)
        x = GTv{k};
        if isstruct(x) && isfield(x,'value')
            okb = ~isfield(x,'bins')  || x.bins  == n_bins;
            oka = ~isfield(x,'alpha') || x.alpha == alpha_val;
            if okb && oka, v = row5(x.value); return; end
        elseif isnumeric(x) && numel(x)>=4
            v = row5(x); return;
        end
    end
    for k=numel(GTv):-1:1
        x = GTv{k};
        if isnumeric(x) && numel(x)>=4, v = row5(x); return; end
    end
    error('Could not parse GroundTruth_value cell.');
else
    v = row5(GTv);
end
end

function r = row5(x)
x = x(:).';
if numel(x)<5, x(5) = NaN; end
r = x(1:5);
end

function b = pick_bias_index_gauss(name)
switch lower(string(name))
    case "plugin",      b = 1;
    case "resample",    b = 2;
    case "shuffle",     b = 3;
    case {"shuff-resamp","shuffresamp","shuff_resamp"}, b = 2; % (averaging not needed along alpha)
    case "venkatesh",   b = 4;
    otherwise, error('Unknown Gaussian bias: %s', name);
end
end
function col = color_for_bias(datatype, bias_name, cmap)
% Returns the RGB row from cmap for a given bias method.
b = lower(string(bias_name));
if strcmp(datatype,'discr')
    % discrete names
    switch b
        case "plugin",         idx = 1;
        case "qe",             idx = 2;
        case {"shuff","shuffsub"}, idx = 3;
        case {"qeshuff","qeShuff","weighted","merged"}, idx = 4;
        case "venkatesh",      idx = 5;
        otherwise,             idx = 6;
    end
else
    % gaussian names
    switch b
        case "plugin",                   idx = 1;
        case "resample",                 idx = 2;
        case {"shuffle","shuffsub"},     idx = 3;
        case {"shuff-resamp","shuffresamp","shuff_resamp","merged"}, idx = 4;
        case "venkatesh",                idx = 5;
        otherwise,                       idx = 6;
    end
end
col = cmap(idx,:);
end
