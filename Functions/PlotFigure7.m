%% ===========================================================
%  Redundancy/Synergy Bias — Discrete (top) + Gaussian (bottom)
%  Figure 1: BIAS (RAW ONLY)
%  Figure 2: BIAS (RAW + Corrections: Resample, QE, Shuffle-Sub)
%  Discrete ground truth: noisy mod-sum (not maximum information)
%  ===========================================================

clear; clc;
rng(1); % Reproducibility
DISCRETE_FILE = 'PID_discrete_bias.mat';
GAUSS_FILE    = 'PID_gaussian_bias.mat';

%% ---------------- GLOBAL PARAMS ----------------
NSOURCES          = 3;
TOTAL_VARS        = NSOURCES + 1;  % sources + target
NUM_REPS          = 100;            % adjust up for final runs
N_SHUFFLES_CORR   = 20;            % resample/QE/shuffle-sub internal shuffles

%% ------------ DISCRETE DATA (NOISY MOD-SUM GROUND TRUTH) ---
% Y = (X1 + X2 + X3 + N) mod R, with P(N=0)=1-eps, P(N=k)=eps/(R-1), k=1..R-1
sample_sizes    = 2.^(10:14);
cardinalities   = 2:6;            % includes R=1 edge case
DISCRETE_EPS    = 0.25;           % <-- noise level (reduces information below max)
nsources        = NSOURCES;
total_vars      = TOTAL_VARS;
num_reps        = NUM_REPS;

nr = numel(sample_sizes); nc = numel(cardinalities);

% Means of *bias* (estimate - ground truth)
red_bias_raw = zeros(nr, nc);
syn_bias_raw = zeros(nr, nc);
red_bias_res = zeros(nr, nc);
syn_bias_res = zeros(nr, nc);
red_bias_qe  = zeros(nr, nc);
syn_bias_qe  = zeros(nr, nc);
red_bias_ss  = zeros(nr, nc);
syn_bias_ss  = zeros(nr, nc);

% Per-repetition *bias* (for SEM)
red_bias_all_raw = zeros(nr, nc, num_reps);
syn_bias_all_raw = zeros(nr, nc, num_reps);
red_bias_all_res = zeros(nr, nc, num_reps);
syn_bias_all_res = zeros(nr, nc, num_reps);
red_bias_all_qe  = zeros(nr, nc, num_reps);
syn_bias_all_qe  = zeros(nr, nc, num_reps);
red_bias_all_ss  = zeros(nr, nc, num_reps);
syn_bias_all_ss  = zeros(nr, nc, num_reps);

if isfile(DISCRETE_FILE)
    S = load(DISCRETE_FILE);

    % restore everything the figures need
    red_bias_raw     = S.red_bias_raw;
    syn_bias_raw     = S.syn_bias_raw;
    red_bias_res     = S.red_bias_res;
    syn_bias_res     = S.syn_bias_res;
    red_bias_qe      = S.red_bias_qe;
    syn_bias_qe      = S.syn_bias_qe;
    red_bias_ss      = S.red_bias_ss;
    syn_bias_ss      = S.syn_bias_ss;

    red_bias_all_raw = S.red_bias_all_raw;
    syn_bias_all_raw = S.syn_bias_all_raw;
    red_bias_all_res = S.red_bias_all_res;
    syn_bias_all_res = S.syn_bias_all_res;
    red_bias_all_qe  = S.red_bias_all_qe;
    syn_bias_all_qe  = S.syn_bias_all_qe;
    red_bias_all_ss  = S.red_bias_all_ss;
    syn_bias_all_ss  = S.syn_bias_all_ss;

    fprintf('[LOAD] DISCRETE results loaded from %s\n', DISCRETE_FILE);

else
    for i = 1:nr
        n = sample_sizes(i);
        for j = 1:nc
            R = cardinalities(j);
    
            % Exact discrete ground truth from closed-form PMF
            [red_star, syn_star] = truth_discr_noisy_modsum(R, DISCRETE_EPS);
    
            rb_raw = zeros(1, num_reps); sb_raw = zeros(1, num_reps);
            rb_res = zeros(1, num_reps); sb_res = zeros(1, num_reps);
            rb_qe  = zeros(1, num_reps); sb_qe  = zeros(1, num_reps);
            rb_ss  = zeros(1, num_reps); sb_ss  = zeros(1, num_reps);
    
            for rep = 1:num_reps
                % Sources iid uniform in {1..R}
                X = randi(R, nsources, n);
    
                % Target via noisy mod-sum
                X0 = X - 1;                               % zero-based
                sum0 = mod(sum(X0,1), R);                 % 1 x n
                if R == 1
                    N = zeros(1,n);                       % degenerate
                else
                    flip = rand(1,n) >= (1 - DISCRETE_EPS);
                    N = zeros(1,n);
                    idx = find(flip);
                    N(idx) = randi(R-1, 1, numel(idx));   % 1..R-1
                end
                Y0 = mod(sum0 + N, R);
                Y  = Y0 + 1;
    
                data = [X; Y].'; % [n x (nsources+1)], last col is Y
    
                % RAW
                [r0, s0] = compute_syn_red_discr(data, nsources, R, total_vars);
    
                % Corrections
                [rR, sR] = resample_syn_red_discr(data, nsources, R, total_vars, N_SHUFFLES_CORR);
                [rQ, sQ] = qe_syn_red_discr(       data, nsources, R, total_vars, N_SHUFFLES_CORR);
                [rS, sS] = shuffsub_syn_red_discr( data, nsources, R, total_vars, N_SHUFFLES_CORR);
    
                % Bias = estimate - ground truth
                rb_raw(rep) = r0 - red_star; sb_raw(rep) = s0 - syn_star;
                rb_res(rep) = rR - red_star; sb_res(rep) = sR - syn_star;
                rb_qe(rep)  = rQ - red_star; sb_qe(rep)  = sQ - syn_star;
                rb_ss(rep)  = rS - red_star; sb_ss(rep)  = sS - syn_star;
            end
    
            % Store per-rep (for SEM) and means
            red_bias_all_raw(i,j,:) = rb_raw; syn_bias_all_raw(i,j,:) = sb_raw;
            red_bias_all_res(i,j,:) = rb_res; syn_bias_all_res(i,j,:) = sb_res;
            red_bias_all_qe(i,j,:)  = rb_qe;  syn_bias_all_qe(i,j,:)  = sb_qe;
            red_bias_all_ss(i,j,:)  = rb_ss;  syn_bias_all_ss(i,j,:)  = sb_ss;
    
            red_bias_raw(i,j) = mean(rb_raw); syn_bias_raw(i,j) = mean(sb_raw);
            red_bias_res(i,j) = mean(rb_res); syn_bias_res(i,j) = mean(sb_res);
            red_bias_qe(i,j)  = mean(rb_qe);  syn_bias_qe(i,j)  = mean(sb_qe);
            red_bias_ss(i,j)  = mean(rb_ss);  syn_bias_ss(i,j)  = mean(sb_ss);
    
            fprintf('[DISCR] n=%d, R=%d -> Bias Red(raw)=%.4f  Syn(raw)=%.4f | Res=%.4f  QE=%.4f  SS=%.4f\n', ...
                n, R, red_bias_raw(i,j), syn_bias_raw(i,j), ...
                syn_bias_res(i,j), syn_bias_qe(i,j), syn_bias_ss(i,j));
        end
    end
    save(DISCRETE_FILE,'-v7.3', ...
            'red_bias_raw','syn_bias_raw','red_bias_res','syn_bias_res', ...
            'red_bias_qe','syn_bias_qe','red_bias_ss','syn_bias_ss', ...
            'red_bias_all_raw','syn_bias_all_raw','red_bias_all_res','syn_bias_all_res', ...
            'red_bias_all_qe','syn_bias_all_qe','red_bias_all_ss','syn_bias_all_ss');

    fprintf('[SAVE] DISCRETE results saved to %s\n', DISCRETE_FILE);
end
%% ------------ GAUSSIAN DATA (KNOWN COV / GROUND TRUTH) -----
% Latent S model per dimension: Xi = S + Ni,  Y = S + Ny
sample_sizes_g = 2.^(6:10);
dims           = [1 3:3:15];
nsources       = NSOURCES;
total_vars     = TOTAL_VARS;
num_reps       = NUM_REPS;

% latent variances
vS  = 1.0; vNx = 1.0; vNy = 1.0;

ng = numel(sample_sizes_g); nd = numel(dims);

% Means of *bias*
red_bias_g_raw = zeros(ng, nd);
syn_bias_g_raw = zeros(ng, nd);
red_bias_g_res = zeros(ng, nd);
syn_bias_g_res = zeros(ng, nd);
red_bias_g_qe  = zeros(ng, nd);
syn_bias_g_qe  = zeros(ng, nd);
red_bias_g_ss  = zeros(ng, nd);
syn_bias_g_ss  = zeros(ng, nd);

% Per-rep *bias* (SEM)
red_bias_all_g_raw = zeros(ng, nd, num_reps);
syn_bias_all_g_raw = zeros(ng, nd, num_reps);
red_bias_all_g_res = zeros(ng, nd, num_reps);
syn_bias_all_g_res = zeros(ng, nd, num_reps);
red_bias_all_g_qe  = zeros(ng, nd, num_reps);
syn_bias_all_g_qe  = zeros(ng, nd, num_reps);
red_bias_all_g_ss  = zeros(ng, nd, num_reps);
syn_bias_all_g_ss  = zeros(ng, nd, num_reps);

if isfile(GAUSS_FILE)
    S = load(GAUSS_FILE);

    % restore everything the figures need
    red_bias_g_raw     = S.red_bias_g_raw;
    syn_bias_g_raw     = S.syn_bias_g_raw;
    red_bias_g_res     = S.red_bias_g_res;
    syn_bias_g_res     = S.syn_bias_g_res;
    red_bias_g_qe      = S.red_bias_g_qe;
    syn_bias_g_qe      = S.syn_bias_g_qe;
    red_bias_g_ss      = S.red_bias_g_ss;
    syn_bias_g_ss      = S.syn_bias_g_ss;

    red_bias_all_g_raw = S.red_bias_all_g_raw;
    syn_bias_all_g_raw = S.syn_bias_all_g_raw;
    red_bias_all_g_res = S.red_bias_all_g_res;
    syn_bias_all_g_res = S.syn_bias_all_g_res;
    red_bias_all_g_qe  = S.red_bias_all_g_qe;
    syn_bias_all_g_qe  = S.syn_bias_all_g_qe;
    red_bias_all_g_ss  = S.red_bias_all_g_ss;
    syn_bias_all_g_ss  = S.syn_bias_all_g_ss;

    fprintf('[LOAD] GAUSSIAN results loaded from %s\n', GAUSS_FILE);

else
    for i = 1:ng
        n = sample_sizes_g(i);
        for j = 1:nd
            d = dims(j);
    
            % Exact block covariance for one dimension
            B = [vS+vNx, vS,      vS,      vS; ...
                 vS,     vS+vNx,  vS,      vS; ...
                 vS,     vS,      vS+vNx,  vS; ...
                 vS,     vS,      vS,      vS+vNy];
            Sigma_true = kron(eye(d), B);  % p x p
            p = size(Sigma_true,1);
    
            % Ground-truth redundancy & synergy from Sigma_true
            src_idx = cell(1,nsources);
            for s = 1:nsources, src_idx{s} = (s-1)*d + (1:d); end
            tgt_idx = nsources*d + (1:d);
    
            mi_single = zeros(1,nsources);
            for s = 1:nsources
                mi_single(s) = gaussianMI(Sigma_true, src_idx{s}, tgt_idx);
            end
            red_star_g = min(mi_single);
    
            pair_mi = [];
            for a = 1:nsources
                for b = a+1:nsources
                    pair_mi(end+1) = gaussianMI(Sigma_true, [src_idx{a}, src_idx{b}], tgt_idx); %#ok<AGROW>
                end
            end
            all_mi   = gaussianMI(Sigma_true, [src_idx{:}], tgt_idx);
            syn_star_g = all_mi - max(pair_mi);
    
            % Run reps and compute bias
            rb_raw = zeros(1, num_reps); sb_raw = zeros(1, num_reps);
            rb_res = zeros(1, num_reps); sb_res = zeros(1, num_reps);
            rb_qe  = zeros(1, num_reps); sb_qe  = zeros(1, num_reps);
            rb_ss  = zeros(1, num_reps); sb_ss  = zeros(1, num_reps);
    
            for rep = 1:num_reps
                data  = mvnrnd(zeros(1,p), Sigma_true, n);
    
                % RAW
                [r0, s0] = compute_syn_red_gauss(data, nsources, d);
    
                % Corrections
                [rR, sR] = resample_syn_red_gauss(data, nsources, d, N_SHUFFLES_CORR);
                [rQ, sQ] = qe_syn_red_gauss(       data, nsources, d, N_SHUFFLES_CORR);
                [rS, sS] = shuffsub_syn_red_gauss( data, nsources, d, N_SHUFFLES_CORR);
    
                rb_raw(rep) = r0 - red_star_g; sb_raw(rep) = s0 - syn_star_g;
                rb_res(rep) = rR - red_star_g; sb_res(rep) = sR - syn_star_g;
                rb_qe(rep)  = rQ - red_star_g; sb_qe(rep)  = sQ - syn_star_g;
                rb_ss(rep)  = rS - red_star_g; sb_ss(rep)  = sS - syn_star_g;
            end
    
            % Save
            red_bias_all_g_raw(i,j,:) = rb_raw; syn_bias_all_g_raw(i,j,:) = sb_raw;
            red_bias_all_g_res(i,j,:) = rb_res; syn_bias_all_g_res(i,j,:) = sb_res;
            red_bias_all_g_qe(i,j,:)  = rb_qe;  syn_bias_all_g_qe(i,j,:)  = sb_qe;
            red_bias_all_g_ss(i,j,:)  = rb_ss;  syn_bias_all_g_ss(i,j,:)  = sb_ss;
    
            red_bias_g_raw(i,j) = mean(rb_raw); syn_bias_g_raw(i,j) = mean(sb_raw);
            red_bias_g_res(i,j) = mean(rb_res); syn_bias_g_res(i,j) = mean(sb_res);
            red_bias_g_qe(i,j)  = mean(rb_qe);  syn_bias_g_qe(i,j)  = mean(sb_qe);
            red_bias_g_ss(i,j)  = mean(rb_ss);  syn_bias_g_ss(i,j)  = mean(sb_ss);
    
            fprintf('[GAUSS] n=%d, d=%d -> Bias Red(raw)=%.4f  Syn(raw)=%.4f | Res=%.4f  QE=%.4f  SS=%.4f\n', ...
                n, d, red_bias_g_raw(i,j), syn_bias_g_raw(i,j), ...
                syn_bias_g_res(i,j), syn_bias_g_qe(i,j), syn_bias_g_ss(i,j));
        end
    end
    save(GAUSS_FILE,'-v7.3', ...
        'red_bias_g_raw','syn_bias_g_raw','red_bias_g_res','syn_bias_g_res', ...
        'red_bias_g_qe','syn_bias_g_qe','red_bias_g_ss','syn_bias_g_ss', ...
        'red_bias_all_g_raw','syn_bias_all_g_raw','red_bias_all_g_res','syn_bias_all_g_res', ...
        'red_bias_all_g_qe','syn_bias_all_g_qe','red_bias_all_g_ss','syn_bias_all_g_ss');

    fprintf('[SAVE] GAUSSIAN results saved to %s\n', GAUSS_FILE);
end
%% --- Figure sizing + typography (A4 textwidth & legible font) ---
textwidth_cm = 17.6;     % typical LaTeX textwidth on A4
figW = textwidth_cm;     
figH = 8;               % adjust if you want a taller grid
baseFont = 'Arial';      % or 'Arial' if you prefer sans-serif

fsAxis  = 8.5;
fsLab   = 9;
fsTitle = 9;
fsLegend= 8;

lw_raw = 1.7; ms = 5; capSize = 6;

% Apply defaults per-figure later with: set(fig,'DefaultAxesFontName',baseFont, ...)


%% -------------------- FIGURE 1: BIAS (RAW ONLY) --------------
fig_raw = figure('Name','PID Bias — Raw Only');
set(fig_raw,'Units','centimeters','Position',[2 2 figW 2*figH],'Color','w', ...
    'DefaultAxesFontName',baseFont,'DefaultTextFontName',baseFont, ...
    'DefaultAxesFontSize',fsAxis);

% Palettes (unchanged)
baseDisc   = [0.00 0.35 0.70];
baseGauss  = [0.85 0.33 0.10];
colDisc_R  = makeSequential(baseDisc,  length(cardinalities));
colDisc_n  = makeSequential(baseDisc,  length(sample_sizes));
colGauss_d = makeSequential(baseGauss, length(dims));
colGauss_n = makeSequential(baseGauss, length(sample_sizes_g));

% Semantic bases
clr_syn = [0.4660 0.6740 0.1880];   % green-ish (synergy)
clr_red = [0.0000 0.4470 0.7410];   % blue-ish  (redundancy)
% --- Tonal palettes (light -> dark) per panel count ---
cols_red_R  = makeSequential(clr_red, length(cardinalities));
cols_syn_R  = makeSequential(clr_syn, length(cardinalities));

cols_red_n  = makeSequential(clr_red, length(sample_sizes));
cols_syn_n  = makeSequential(clr_syn, length(sample_sizes));

cols_red_d  = makeSequential(clr_red, length(dims));            % Gaussian: by d
cols_syn_d  = makeSequential(clr_syn, length(dims));

cols_red_ng = makeSequential(clr_red, length(sample_sizes_g));  % Gaussian: by n
cols_syn_ng = makeSequential(clr_syn, length(sample_sizes_g));


t1 = tiledlayout(fig_raw,4, 4, 'TileSpacing','compact', 'Padding','compact');

% ===================== TOP ROW: DISCRETE ======================
% 1) Discrete: Red bias vs n
ax1 = nexttile(1); hold on;
for j = 1:length(cardinalities)
    c = cols_red_R(j,:);
    sem = squeeze(std(red_bias_all_raw(:, j, :), 0, 3)) / sqrt(num_reps);
    % plotMeanBand(sample_sizes, red_bias_raw(:, j), 3*sem, '-', ...
    %     'LineWidth', lw_raw, 'MarkerSize', ms, 'CapSize', capSize, ...
    %     'Color', c, 'MarkerFaceColor', c, 'MarkerEdgeColor', c);
    plotMeanBand(sample_sizes, red_bias_raw(:, j), sem, c, '-', lw_raw, ms, 0.20);

end
yline(0,'k-'); set(gca,'XScale','log');
xlabel('Sample Size','FontSize',fsLab); ylabel('Bias [bits]','FontSize',fsLab);
title('Discrete Redundancy','FontSize',fsTitle);
legend(arrayfun(@(R) sprintf('R = %d', R), cardinalities, 'UniformOutput', false), ...
       'Location','best','Box','off','FontSize',fsLegend); grid off; hold off;

% 2) Discrete: Syn bias vs n
ax2 = nexttile(2); hold on;
for j = 1:length(cardinalities)
    c =  cols_syn_R(j,:);
    sem = squeeze(std(syn_bias_all_raw(:, j, :), 0, 3)) / sqrt(num_reps);
    % plotMeanBand(sample_sizes, syn_bias_raw(:, j), 3*sem, '-', ...
    %     'LineWidth', lw_raw, 'MarkerSize', ms, 'CapSize', capSize, ...
    %     'Color', c, 'MarkerFaceColor', c, 'MarkerEdgeColor', c);
    plotMeanBand(sample_sizes, syn_bias_raw(:, j), sem, c, '-', lw_raw, ms, 0.20);

end
yline(0,'k-'); set(gca,'XScale','log');
xlabel('Sample Size','FontSize',fsLab);
ylabel([],'FontSize',fsLab);
ylim([0, 0.5]);%yticklabels([]);
title('Discrete Synergy','FontSize',fsTitle);
legend(arrayfun(@(R) sprintf('R = %d', R), cardinalities, 'UniformOutput', false), ...
       'Location','best','Box','off','FontSize',fsLegend); grid off; hold off;

% 3) Discrete: Red bias vs R
ax3 = nexttile(3); hold on;
for i = 1:length(sample_sizes)
    c = cols_red_n(i,:);
    sem = squeeze(std(red_bias_all_raw(i, :, :), 0, 3))' / sqrt(num_reps);
    % plotMeanBand(cardinalities, red_bias_raw(i, :), 3*sem, '-', ...
    %     'LineWidth', lw_raw, 'MarkerSize', ms, 'CapSize', capSize, ...
    %     'Color', c, 'MarkerFaceColor', c, 'MarkerEdgeColor', c);
    plotMeanBand(cardinalities, red_bias_raw(i, :)', sem, c, '-', lw_raw, ms, 0.20);
end
yline(0,'k-');
xlim([min(cardinalities), max(cardinalities)]);
xlabel('R','FontSize',fsLab); %yticklabels([]);
title('Discrete Redundancy','FontSize',fsTitle);
legend(arrayfun(@(s) n2label(s), sample_sizes, 'UniformOutput', false), ...
       'Location','best','Box','off','FontSize',fsLegend); grid off; hold off;

% 4) Discrete: Syn bias vs R
ax4 = nexttile(4); hold on;
for i = 1:length(sample_sizes)
    c =  cols_syn_n(i,:);
    sem = squeeze(std(syn_bias_all_raw(i, :, :), 0, 3))' / sqrt(num_reps);
    % plotMeanBand(cardinalities, syn_bias_raw(i, :), 3*sem, '-', ...
    %     'LineWidth', lw_raw, 'MarkerSize', ms, 'CapSize', capSize, ...
    %     'Color', c, 'MarkerFaceColor', c, 'MarkerEdgeColor', c);
    plotMeanBand(cardinalities, syn_bias_raw(i, :)', sem, c, '-', lw_raw, ms, 0.20);

end
yline(0,'k-');
xlim([min(cardinalities), max(cardinalities)]);
xlabel('R','FontSize',fsLab);
ylabel([],'FontSize',fsLab);
ylim([0, 0.5]);%yticklabels([]);
linkaxes([ax1 ax2 ax3 ax4],'y')
yticklabels(ax2,[])
yticklabels(ax3,[])
yticklabels(ax4,[])
title('Discrete Synergy','FontSize',fsTitle);
legend(arrayfun(@(s) n2label(s), sample_sizes, 'UniformOutput', false), ...
       'Location','best','Box','off','FontSize',fsLegend); grid off; hold off;












% =================== Third ROW: GAUSSIAN =====================
% 5) Gaussian: Red bias vs n
ax9 = nexttile(5+4); hold on;
for j = 1:length(dims)
    c = cols_red_d(j,:);
    sem = squeeze(std(red_bias_all_g_raw(:, j, :), 0, 3)) / sqrt(num_reps);
    % plotMeanBand(sample_sizes_g, red_bias_g_raw(:, j), 3*sem, '-', ...
    %     'LineWidth', lw_raw, 'MarkerSize', ms, 'CapSize', capSize, ...
    %     'Color', c, 'MarkerFaceColor', c, 'MarkerEdgeColor', c);
    plotMeanBand(sample_sizes_g, red_bias_g_raw(:, j), sem, c, '-', lw_raw, ms, 0.20);

end
yline(0,'k-'); set(gca,'XScale','log');
xlabel('Sample Size','FontSize',fsLab);
ylabel('Bias [bits]','FontSize',fsLab);
xticks([100, 1000])
title('Gaussian Redundancy','FontSize',fsTitle);
legend(arrayfun(@(d) sprintf('d = %d', d), dims, 'UniformOutput', false), ...
       'Location','best','Box','off','FontSize',fsLegend); grid off; hold off;

% 6) Gaussian: Syn bias vs n
ax10 = nexttile(6+4); hold on;
for j = 1:length(dims)
    c =  cols_syn_d(j,:);
    sem = squeeze(std(syn_bias_all_g_raw(:, j, :), 0, 3)) / sqrt(num_reps);
    % plotMeanBand(sample_sizes_g, syn_bias_g_raw(:, j), 3*sem, '-', ...
    %     'LineWidth', lw_raw, 'MarkerSize', ms, 'CapSize', capSize, ...
    %     'Color', c, 'MarkerFaceColor', c, 'MarkerEdgeColor', c);
    plotMeanBand(sample_sizes_g, syn_bias_g_raw(:, j), sem, c, '-', lw_raw, ms, 0.20);

end
yline(0,'k-'); set(gca,'XScale','log');
xlabel('Sample Size','FontSize',fsLab);
xticks([100, 1000])
ylabel([],'FontSize',fsLab);%yticklabels([]);
title('Gaussian Synergy','FontSize',fsTitle);
legend(arrayfun(@(d) sprintf('d = %d', d), dims, 'UniformOutput', false), ...
       'Location','best','Box','off','FontSize',fsLegend); grid off; hold off;

% 7) Gaussian: Red bias vs d
ax11 = nexttile(7+4); hold on;
for i = 1:length(sample_sizes_g)
    c = cols_red_ng(i,:);
    sem = squeeze(std(red_bias_all_g_raw(i, :, :), 0, 3))' / sqrt(num_reps);
    % plotMeanBand(dims, red_bias_g_raw(i, :), 3*sem, '-', ...
    %     'LineWidth', lw_raw, 'MarkerSize', ms, 'CapSize', capSize, ...
    %     'Color', c, 'MarkerFaceColor', c, 'MarkerEdgeColor', c);
    plotMeanBand(dims, red_bias_g_raw(i, :)', sem, c, '-', lw_raw, ms, 0.20);
end
yline(0,'k-');
xlim([min(dims), max(dims)]);
xlabel('d','FontSize',fsLab);
ylabel([],'FontSize',fsLab);%yticklabels([]);
title('Gaussian Redundancy','FontSize',fsTitle);
legend(arrayfun(@(s) n2label(s), sample_sizes_g, 'UniformOutput', false), ...
       'Location','best','Box','off','FontSize',fsLegend); grid off; hold off;

% 8) Gaussian: Syn bias vs d
ax12 = nexttile(8+4); hold on;
for i = 1:length(sample_sizes_g)
    c =  cols_syn_ng(i,:);
    sem = squeeze(std(syn_bias_all_g_raw(i, :, :), 0, 3))' / sqrt(num_reps);
    % plotMeanBand(dims, syn_bias_g_raw(i, :), 3*sem, '-', ...
    %     'LineWidth', lw_raw, 'MarkerSize', ms, 'CapSize', capSize, ...
    %     'Color', c, 'MarkerFaceColor', c, 'MarkerEdgeColor', c);
    plotMeanBand(dims, syn_bias_g_raw(i, :)', sem, c, '-', lw_raw, ms, 0.20);

end
yline(0,'k-');
xlim([min(dims), max(dims)]);
xlabel('d','FontSize',fsLab);
ylabel([],'FontSize',fsLab);%yticklabels([]);
linkaxes([ax9, ax10, ax11, ax12], 'y')
yticklabels(ax10,[])
yticklabels(ax11,[])
yticklabels(ax12,[])

title('Gaussian Synergy','FontSize',fsTitle);
legend(arrayfun(@(s) n2label(s), sample_sizes_g, 'UniformOutput', false), ...
       'Location','best','Box','off','FontSize',fsLegend); grid off; hold off;

% Export RAW-BIAS at textwidth
% exportgraphics(fig_raw,'PID_Bias_RawOnly.pdf','ContentType','vector');
% saveas(fig_raw,'PID_Bias_RawOnly', 'svg', 'Resolution')


% exportgraphics(fig_raw,'PID_Bias_RawOnly.png','Resolution',600);
% set(fig_raw,'Renderer','painters');
% print(fig_raw,'-dsvg','PID_Bias_RawOnly');

%% --------------- FIGURE 2: BIAS (ALL CORRECTIONS) -----------
% fig_corr = figure('Name','PID Bias — All Corrections');
% set(fig_corr,'Units','centimeters','Position',[2 2 figW figH],'Color','w', ...
%     'DefaultAxesFontName',baseFont,'DefaultTextFontName',baseFont, ...
%     'DefaultAxesFontSize',fsAxis);

lw_alt = 1.2; % for corrected curves

% helper: Discrete (QE ':' , SS '-.' , merged '-' )
% plotAlldiscr = @(x, y_raw, sem_raw, y_res, sem_res, y_qe, sem_qe, y_ss, sem_ss, color) [ ...
%     plotMeanBand(x, y_qe,  3*sem_qe,  ':s',  'LineWidth', lw_alt, 'MarkerSize', ms, 'CapSize', capSize, ...
%              'Color', color, 'MarkerFaceColor', 'none',  'MarkerEdgeColor', color,'HandleVisibility','off'); ...
%     plotMeanBand(x, y_ss,  3*sem_ss,  '-.d', 'LineWidth', lw_alt, 'MarkerSize', ms, 'CapSize', capSize, ...
%              'Color', color, 'MarkerFaceColor', 'none',  'MarkerEdgeColor', color,'HandleVisibility','off'); ...
%     plotMeanBand(x, (y_qe+y_ss)/2,  3*(sem_ss+sem_qe)/2,  '-', 'LineWidth', lw_alt, 'MarkerSize', ms, 'CapSize', capSize, ...
%              'Color', color, 'MarkerFaceColor', 'none',  'MarkerEdgeColor', color); ...         
% ];
% 
% % helper: Gaussian (Res '--' , SS '-.' , merged '-' )
% plotAllgauss = @(x, y_raw, sem_raw, y_res, sem_res, y_qe, sem_qe, y_ss, sem_ss, color) [ ...
%     plotMeanBand(x, y_res, 3*sem_res, '--o', 'LineWidth', lw_alt, 'MarkerSize', ms, 'CapSize', capSize, ...
%              'Color', color, 'MarkerFaceColor', 'none',  'MarkerEdgeColor', color,'HandleVisibility','off'); ...
%     plotMeanBand(x, y_ss,  3*sem_ss,  '-.d', 'LineWidth', lw_alt, 'MarkerSize', ms, 'CapSize', capSize, ...
%              'Color', color, 'MarkerFaceColor', 'none',  'MarkerEdgeColor', color,'HandleVisibility','off'); ...
%     plotMeanBand(x, (y_res+y_ss)/2,  3*(sem_ss+sem_res)/2,  '-', 'LineWidth', lw_alt, 'MarkerSize', ms, 'CapSize', capSize, ...
%              'Color', color, 'MarkerFaceColor', 'none',  'MarkerEdgeColor', color); ...         
% ];

% Discrete (QE ':' , SS '-.' , merged '-')
plotAlldiscr = @(x, y_raw, sem_raw, y_res, sem_res, y_qe, sem_qe, y_ss, sem_ss, color) [ ...
    plotMeanBand(x, y_qe, sem_qe, color,  ':',  lw_alt, ms, 0.18); ...
    plotMeanBand(x, y_ss, sem_ss, color, '-.',  lw_alt, ms, 0.18); ...
    plotMeanBand(x, (y_qe+y_ss)/2, (sem_ss+sem_qe)/2, color, '-', lw_alt, ms, 0.22) ...
];

% Gaussian (Res '--' , SS '-.' , merged '-')
plotAllgauss = @(x, y_raw, sem_raw, y_res, sem_res, y_qe, sem_qe, y_ss, sem_ss, color) [ ...
    plotMeanBand(x, y_res, sem_res, color, '--', lw_alt, ms, 0.18); ...
    plotMeanBand(x, y_ss, sem_ss, color, '-.',  lw_alt, ms, 0.18); ...
    plotMeanBand(x, (y_res+y_ss)/2, (sem_ss+sem_res)/2, color, '-', lw_alt, ms, 0.22) ...
];

% ---- picks (unchanged) ----
idxR4     = find(cardinalities==4, 1);      % R = 4
idxN_2p11 = find(sample_sizes==2^11, 1);    % n = 2^11
idxd12    = find(dims==12, 1);              % d = 12
idxNg_2p8 = find(sample_sizes_g==2^8, 1);   % n = 2^8

% t2 = tiledlayout(fig_corr, 2, 4, 'TileSpacing','compact', 'Padding','compact');

% ===================== TOP ROW: DISCRETE ======================
% Red bias vs n (R=4)
ax5 = nexttile(1+4); hold on;
j = idxR4; c = clr_red;%colDisc_R(j,:);
sem_raw = squeeze(std(red_bias_all_raw(:, j, :), 0, 3)) / sqrt(num_reps);
sem_res = squeeze(std(red_bias_all_res(:, j, :), 0, 3)) / sqrt(num_reps);
sem_qe  = squeeze(std(red_bias_all_qe(:,  j, :), 0, 3)) / sqrt(num_reps);
sem_ss  = squeeze(std(red_bias_all_ss(:,  j, :), 0, 3)) / sqrt(num_reps);
plotAlldiscr(sample_sizes, red_bias_raw(:, j), sem_raw, ...
                      red_bias_res(:, j), sem_res, ...
                      red_bias_qe(:,  j), sem_qe, ...
                      red_bias_ss(:,  j), sem_ss, c);
yline(0,'k-'); set(gca,'XScale','log');
xlabel('Sample Size','FontSize',fsLab); ylabel('Bias [bits]','FontSize',fsLab);
title('Discrete Redundancy','FontSize',fsTitle);
legend({'R = 4'}, 'Location','best','Box','off','FontSize',fsLegend); grid off; hold off;

% Syn bias vs n (R=4)
ax6 = nexttile(2+4); hold on; c = clr_syn;
plotAlldiscr(sample_sizes, syn_bias_raw(:, j), sem_raw, ...
                      syn_bias_res(:, j), sem_res, ...
                      syn_bias_qe(:,  j), sem_qe, ...
                      syn_bias_ss(:,  j), sem_ss, c);
yline(0,'k-'); set(gca,'XScale','log');
xlabel('Sample Size','FontSize',fsLab); ylabel([],'FontSize',fsLab);
% %yticklabels([]);
title('Discrete Synergy','FontSize',fsTitle);
legend({'R = 4'}, 'Location','best','Box','off','FontSize',fsLegend); grid off; hold off;

% Red bias vs R (n=2^11)
ax7 = nexttile(3+4); hold on;
i = idxN_2p11; c = clr_red;%colDisc_n(i,:);
sem_raw = squeeze(std(red_bias_all_raw(i, :, :), 0, 3))' / sqrt(num_reps);
sem_res = squeeze(std(red_bias_all_res(i, :, :), 0, 3))' / sqrt(num_reps);
sem_qe  = squeeze(std(red_bias_all_qe(i,  :, :), 0, 3))' / sqrt(num_reps);
sem_ss  = squeeze(std(red_bias_all_ss(i,  :, :), 0, 3))' / sqrt(num_reps);
plotAlldiscr(cardinalities, red_bias_raw(i, :), sem_raw, ...
                     red_bias_res(i, :), sem_res, ...
                       red_bias_qe(i,  :), sem_qe, ...
                       red_bias_ss(i,  :), sem_ss, c);
yline(0,'k-');
xlim([min(cardinalities), max(cardinalities)]);
xlabel('R','FontSize',fsLab); 
%yticklabels([]);
title('Discrete Redundancy','FontSize',fsTitle);
legend({n2label(2^11)}, 'Location','best','Box','off','FontSize',fsLegend); grid off; hold off;

% Syn bias vs R (n=2^11)
ax8 = nexttile(4+4); hold on;c = clr_syn;%
plotAlldiscr(cardinalities, syn_bias_raw(i, :), sem_raw, ...
                       syn_bias_res(i, :), sem_res, ...
                       syn_bias_qe(i,  :), sem_qe, ...
                       syn_bias_ss(i,  :), sem_ss, c);
yline(0,'k-');
xlim([min(cardinalities), max(cardinalities)]);
xlabel('R','FontSize',fsLab); ylabel([],'FontSize',fsLab);
%yticklabels([]);
title('Discrete Synergy','FontSize',fsTitle);
linkaxes([ax5 ax6 ax7 ax8], 'y')
yticklabels(ax6,[])
yticklabels(ax7,[])
yticklabels(ax8,[])
legend({n2label(2^11)}, 'Location','best','Box','off','FontSize',fsLegend); grid off; hold off;

% =================== BOTTOM ROW: GAUSSIAN =====================
% Red bias vs n (d=12)
ax13 = nexttile(5+8); hold on;
j = idxd12; c = clr_red ;%colGauss_d(j,:);
sem_raw = squeeze(std(red_bias_all_g_raw(:, j, :), 0, 3)) / sqrt(num_reps);
sem_res = squeeze(std(red_bias_all_g_res(:, j, :), 0, 3)) / sqrt(num_reps);
sem_qe  = squeeze(std(red_bias_all_g_qe(:,  j, :), 0, 3)) / sqrt(num_reps);
sem_ss  = squeeze(std(red_bias_all_g_ss(:,  j, :), 0, 3)) / sqrt(num_reps);
plotAllgauss(sample_sizes_g, red_bias_g_raw(:, j), sem_raw, ...
                         red_bias_g_res(:, j), sem_res, ...
                         red_bias_g_qe(:,  j), sem_qe, ...
                         red_bias_g_ss(:,  j), sem_ss, c);
yline(0,'k-'); set(gca,'XScale','log');
xlabel('Sample Size','FontSize',fsLab); ylabel('Bias [bits]','FontSize',fsLab);
xticks([100, 1000])
title('Gaussian Redundancy','FontSize',fsTitle);
legend({'d = 12'}, 'Location','best','Box','off','FontSize',fsLegend); grid off; hold off;

% Syn bias vs n (d=12)  **BUGFIX: use syn_bias_g_* arrays (not *_all_)**
ax14 = nexttile(6+8); hold on;c = clr_syn;
plotAllgauss(sample_sizes_g, syn_bias_g_raw(:, j), sem_raw, ...
                         syn_bias_g_res(:, j), sem_res, ...
                         syn_bias_g_qe(:,  j), sem_qe, ...
                         syn_bias_g_ss(:,  j), sem_ss, c);
yline(0,'k-'); set(gca,'XScale','log');
xticks([100, 1000])
xlabel('Sample Size','FontSize',fsLab); ylabel([],'FontSize',fsLab);%yticklabels([]);
title('Gaussian Synergy','FontSize',fsTitle);
legend({'d = 12'}, 'Location','best','Box','off','FontSize',fsLegend); grid off; hold off;

% Red bias vs d (n=2^8)
ax15 = nexttile(7+8); hold on;
i = idxNg_2p8; c =clr_red;% colGauss_n(i,:);
sem_raw = squeeze(std(red_bias_all_g_raw(i, :, :), 0, 3))' / sqrt(num_reps);
sem_res = squeeze(std(red_bias_all_g_res(i, :, :), 0, 3))' / sqrt(num_reps);
sem_qe  = squeeze(std(red_bias_all_g_qe(i,  :, :), 0, 3))' / sqrt(num_reps);
sem_ss  = squeeze(std(red_bias_all_g_ss(i,  :, :), 0, 3))' / sqrt(num_reps);
plotAllgauss(dims, red_bias_g_raw(i, :), sem_raw, ...
                        red_bias_g_res(i, :), sem_res, ...
                        red_bias_g_qe(i,  :), sem_qe, ...
                        red_bias_g_ss(i,  :), sem_ss, c);
yline(0,'k-');
xlim([min(dims), max(dims)]);
xlabel('d','FontSize',fsLab); ylabel([],'FontSize',fsLab);%yticklabels([]);
title('Gaussian Redundancy','FontSize',fsTitle);
legend({n2label(2^8)}, 'Location','best','Box','off','FontSize',fsLegend); grid off; hold off;

% Syn bias vs d (n=2^8)
ax16 = nexttile(8+8); hold on; c = clr_syn;
plotAllgauss(dims, syn_bias_g_raw(i, :), sem_raw, ...
                        syn_bias_g_res(i, :), sem_res, ...
                        syn_bias_g_qe(i,  :), sem_qe, ...
                        syn_bias_g_ss(i,  :), sem_ss, c);
yline(0,'k-');
xlim([min(dims), max(dims)]);
xlabel('d','FontSize',fsLab); ylabel([],'FontSize',fsLab);%yticklabels([]);
title('Gaussian Synergy','FontSize',fsTitle);
linkaxes([ax13 ax14 ax15 ax16], 'y')
yticklabels(ax14,[])
yticklabels(ax15,[])
yticklabels(ax16,[])
legend({n2label(2^8)}, 'Location','best','Box','off','FontSize',fsLegend); grid off; hold off;

% % Line-style legend (unchanged text but label order matches helpers)
% axes('Position',[0 0 1 1],'Visible','off'); hold on;
% plot(nan,nan,'--','Color','k','LineWidth',lw_alt);   % Resample
% plot(nan,nan,':','Color','k','LineWidth',lw_alt);    % QE
% plot(nan,nan,'-.','Color','k','LineWidth',lw_alt);   % Shuffle-Sub
% plot(nan,nan,'-','Color','k','LineWidth',lw_raw);    % Merged
% legend({'Resample','QE','Shuffle-Sub','Merged'}, 'Box','off', ...
%        'Orientation','horizontal', 'FontSize', fsLegend, ...
%        'Position',[0.34 0.96 0.32 0.03]);

% ---- Simple, non-overlapping legend at the bottom ----
figure(fig_raw);  % make sure we're targeting the right figure
axes('Position',[0 0 1 1],'Visible','off'); % placeholder, not overlay
hold on;

% Create dummy lines for legend
h1 = plot(nan,nan,'--k','LineWidth',lw_alt);
h2 = plot(nan,nan,':k','LineWidth',lw_alt);
h3 = plot(nan,nan,'-.k','LineWidth',lw_alt);
h4 = plot(nan,nan,'-k','LineWidth',lw_raw);

legend([h1 h2 h3 h4], ...
       {'Resample','QE','Shuffle-Sub','Merged'}, ...
       'Orientation','horizontal', ...
       'Location','southoutside', ...
       'Box','off', ...
       'FontSize',fsLegend);

hold off;
outDir = fullfile(fileparts(pwd), 'Figures_mat');  

if ~exist(outDir, 'dir')
    mkdir(outDir);
end

exportgraphics(gcf, fullfile(outDir, 'Figure_7.svg'), 'ContentType','vector');



%% -------------------- FUNCTIONS ----------------------------
function MI = gaussianMI(cov_mat, source_idx, target_idx)
    idx = [source_idx, target_idx];
    cov_xy = cov_mat(idx, idx);
    cov_x  = cov_mat(source_idx, source_idx);
    cov_y  = cov_mat(target_idx, target_idx);
    MI = 0.5 * log2(det(cov_x) * det(cov_y) / det(cov_xy));
end

function MI = mutualInformationLast(p)
    p = p / sum(p(:));
    dims = ndims(p);
    z_dim = dims;
    x_dims = 1:(dims-1);
    pz = sum(p, x_dims);
    px = sum(p, z_dim);
    px_exp = repmat(px, [ones(1, dims-1), size(p, dims)]);
    pz_exp = repmat(reshape(pz, [ones(1, dims-1), numel(pz)]), size(px));
    valid = p > eps & px_exp > eps & pz_exp > eps;
    MI = sum(p(valid) .* log2(p(valid) ./ (px_exp(valid) .* pz_exp(valid))));
end

function [red, syn] = compute_syn_red_discr(input_data_t, nsources, card, total_vars)
    counts = accumarray(input_data_t, 1, repmat(card, 1, total_vars));
    % REDUNDANCY: min I(Xi; Y)
    mi_values = zeros(1, nsources);
    for s = 1:nsources
        dims = sort([s, total_vars]);  % Xi and Y
        joint_counts = sum(counts, setdiff(1:total_vars, dims));
        mi_values(s) = mutualInformationLast(joint_counts);
    end
    red = min(mi_values);
    % SYNERGY: I(all sources; Y) - max I(Xi, Xj; Y)
    info_all = mutualInformationLast(counts);
    pairwise_mis = [];
    for a = 1:nsources
        for b = a+1:nsources
            dims = sort([a, b, total_vars]);
            joint_counts_pair = sum(counts, setdiff(1:total_vars, dims));
            pairwise_mis(end+1) = mutualInformationLast(joint_counts_pair);
        end
    end
    syn = info_all - max(pairwise_mis);
end

function [red, syn] = compute_syn_red_gauss(data, nsources, dim_per_var)
    cov_est = cov(data);
    total_dims = (nsources + 1) * dim_per_var;
    % Indices
    source_idxs = cell(1, nsources);
    for s = 1:nsources
        source_idxs{s} = ((s-1)*dim_per_var + 1):(s*dim_per_var);
    end
    target_idx = (nsources*dim_per_var + 1):total_dims;
    % REDUNDANCY: min MI(Xi; Y)
    mi_vals = zeros(1, nsources);
    for s = 1:nsources
        mi_vals(s) = gaussianMI(cov_est, source_idxs{s}, target_idx);
    end
    red = min(mi_vals);
    % SYNERGY: MI(all sources; Y) - max MI(Xi,Xj; Y)
    all_sources_idx = [source_idxs{:}];
    info_all = gaussianMI(cov_est, all_sources_idx, target_idx);
    pairwise_mis = [];
    for a = 1:nsources
        for b = a+1:nsources
            pair_idx = [source_idxs{a}, source_idxs{b}];
            pairwise_mis(end+1) = gaussianMI(cov_est, pair_idx, target_idx);
        end
    end
    syn = info_all - max(pairwise_mis);
end

function cols = makeSequential(baseRGB, n)
    w = linspace(0.15, 0.85, n);
    cols = zeros(n,3);
    for k = 1:n
        cols(k,:) = (1 - w(k))*baseRGB + w(k)*[1 1 1];
    end
    cols = cols(end:-1:1, :);
end

function s = n2label(n)
    k = log2(n);
    superscripts = '⁰¹²³⁴⁵⁶⁷⁸⁹';
    k_str = num2str(k);
    s = 'n = 2';
    for ch = k_str
        s = [s superscripts(str2double(ch)+1)];
    end
end

%% --------- Corrections (Gaussian) ---------
function [red_corr, syn_corr] = resample_syn_red_gauss(data, nsources, dim_per_var, nshuffles)
    if nargin < 4 || isempty(nshuffles), nshuffles = 200; end
    [nTrials, p] = size(data);
    mu = mean(data, 1); Sigma = cov(data);
    R = chol_with_jitter(Sigma); % lower

    [red_obs, syn_obs] = compute_syn_red_gauss(data, nsources, dim_per_var);

    red_sim = zeros(1, nshuffles);
    syn_sim = zeros(1, nshuffles);
    for k = 1:nshuffles
        Z = randn(nTrials, p) * R' + mu;
        [r_k, s_k] = compute_syn_red_gauss(Z, nsources, dim_per_var);
        red_sim(k) = r_k; syn_sim(k) = s_k;
    end
    red_corr = 2*red_obs - mean(red_sim);
    syn_corr = 2*syn_obs - mean(syn_sim);
end

function [red_corr, syn_corr] = qe_syn_red_gauss(data, nsources, dim_per_var, nshuffles)
    if nargin < 4 || isempty(nshuffles), nshuffles = 200; end
    n = size(data,1);
    [red_full, syn_full] = compute_syn_red_gauss(data, nsources, dim_per_var);

    red2 = zeros(2, nshuffles); syn2 = zeros(2, nshuffles);
    red4 = zeros(4, nshuffles); syn4 = zeros(4, nshuffles);

    for k = 1:nshuffles
        idx = randperm(n);
        n2 = floor(n/2); h1 = idx(1:n2); h2 = idx(n2+1:2*n2);
        [red2(1,k), syn2(1,k)] = compute_syn_red_gauss(data(h1,:), nsources, dim_per_var);
        [red2(2,k), syn2(2,k)] = compute_syn_red_gauss(data(h2,:), nsources, dim_per_var);

        n4 = floor(n/4);
        q1 = idx(1:n4); q2 = idx(n4+1:2*n4); q3 = idx(2*n4+1:3*n4); q4 = idx(3*n4+1:4*n4);
        [red4(1,k), syn4(1,k)] = compute_syn_red_gauss(data(q1,:), nsources, dim_per_var);
        [red4(2,k), syn4(2,k)] = compute_syn_red_gauss(data(q2,:), nsources, dim_per_var);
        [red4(3,k), syn4(3,k)] = compute_syn_red_gauss(data(q3,:), nsources, dim_per_var);
        [red4(4,k), syn4(4,k)] = compute_syn_red_gauss(data(q4,:), nsources, dim_per_var);
    end

    m_red2 = mean(red2, 'all'); m_syn2 = mean(syn2, 'all');
    m_red4 = mean(red4, 'all'); m_syn4 = mean(syn4, 'all');

    x = [1/n, 2/n, 4/n];
    p_red = polyfit(x, [red_full, m_red2, m_red4], 2);
    p_syn = polyfit(x, [syn_full, m_syn2, m_syn4], 2);

    red_corr = p_red(3); syn_corr = p_syn(3);
end

function [red_corr, syn_corr] = shuffsub_syn_red_gauss(data, nsources, dim_per_var, nshuffles)
    if nargin < 4 || isempty(nshuffles), nshuffles = 200; end
    [nTrials, p] = size(data);
    [red_obs, syn_obs] = compute_syn_red_gauss(data, nsources, dim_per_var);
    red_null = zeros(1, nshuffles); syn_null = zeros(1, nshuffles);

    tgt_cols = (nsources*dim_per_var + 1) : p;
    for k = 1:nshuffles
        rp = randperm(nTrials);
        Z = data; Z(:, tgt_cols) = data(rp, tgt_cols);
        [r_k, s_k] = compute_syn_red_gauss(Z, nsources, dim_per_var);
        red_null(k) = r_k; syn_null(k) = s_k;
    end
    red_corr = red_obs - mean(red_null);
    syn_corr = syn_obs - mean(syn_null);
end

%% --------- Corrections (Discrete) ----------
function [red_corr, syn_corr] = resample_syn_red_discr(data_t, nsources, card, total_vars, nshuffles)
    if nargin < 5 || isempty(nshuffles), nshuffles = 200; end
    n = size(data_t,1);

    [red_obs, syn_obs] = compute_syn_red_discr(data_t, nsources, card, total_vars);

    dimsizes = repmat(card, 1, total_vars);
    counts = accumarray(data_t, 1, dimsizes);
    pmf = counts / sum(counts(:));
    w = pmf(:); M = numel(w);

    red_sim = zeros(1, nshuffles);
    syn_sim = zeros(1, nshuffles);

    for k = 1:nshuffles
        idx = randsample(M, n, true, w);
        [subs{1:total_vars}] = ind2sub(dimsizes, idx);
        sim_t = zeros(n, total_vars);
        for c = 1:total_vars, sim_t(:, c) = subs{c}(:); end
        [r, s] = compute_syn_red_discr(sim_t, nsources, card, total_vars);
        red_sim(k) = r; syn_sim(k) = s;
    end

    red_corr = 2*red_obs - mean(red_sim);
    syn_corr = 2*syn_obs - mean(syn_sim);
end

function [red_corr, syn_corr] = qe_syn_red_discr(data_t, nsources, card, total_vars, nshuffles)
    if nargin < 5 || isempty(nshuffles), nshuffles = 200; end
    n = size(data_t,1);
    [red_full, syn_full] = compute_syn_red_discr(data_t, nsources, card, total_vars);

    red2 = zeros(2, nshuffles); syn2 = zeros(2, nshuffles);
    red4 = zeros(4, nshuffles); syn4 = zeros(4, nshuffles);

    for k = 1:nshuffles
        idx = randperm(n);

        n2 = floor(n/2); h1 = idx(1:n2); h2 = idx(n2+1:2*n2);
        [red2(1,k), syn2(1,k)] = compute_syn_red_discr(data_t(h1,:), nsources, card, total_vars);
        [red2(2,k), syn2(2,k)] = compute_syn_red_discr(data_t(h2,:), nsources, card, total_vars);

        n4 = floor(n/4);
        q1 = idx(1:n4); q2 = idx(n4+1:2*n4); q3 = idx(2*n4+1:3*n4); q4 = idx(3*n4+1:4*n4);
        [red4(1,k), syn4(1,k)] = compute_syn_red_discr(data_t(q1,:), nsources, card, total_vars);
        [red4(2,k), syn4(2,k)] = compute_syn_red_discr(data_t(q2,:), nsources, card, total_vars);
        [red4(3,k), syn4(3,k)] = compute_syn_red_discr(data_t(q3,:), nsources, card, total_vars);
        [red4(4,k), syn4(4,k)] = compute_syn_red_discr(data_t(q4,:), nsources, card, total_vars);
    end

    m_red2 = mean(red2, 'all'); m_syn2 = mean(syn2, 'all');
    m_red4 = mean(red4, 'all'); m_syn4 = mean(syn4, 'all');

    x = [1/n, 2/n, 4/n];
    p_red = polyfit(x, [red_full, m_red2, m_red4], 2);
    p_syn = polyfit(x, [syn_full, m_syn2, m_syn4], 2);

    red_corr = p_red(3); syn_corr = p_syn(3);
end

function [red_corr, syn_corr] = shuffsub_syn_red_discr(data_t, nsources, card, total_vars, nshuffles)
    if nargin < 5 || isempty(nshuffles), nshuffles = 200; end
    n = size(data_t,1);
    [red_obs, syn_obs] = compute_syn_red_discr(data_t, nsources, card, total_vars);

    red_null = zeros(1, nshuffles); syn_null = zeros(1, nshuffles);
    tgt_col = total_vars; % last column as target

    for k = 1:nshuffles
        rp = randperm(n);
        Z = data_t; Z(:, tgt_col) = data_t(rp, tgt_col);
        [r, s] = compute_syn_red_discr(Z, nsources, card, total_vars);
        red_null(k) = r; syn_null(k) = s;
    end

    red_corr = red_obs - mean(red_null);
    syn_corr = syn_obs - mean(syn_null);
end

% ---------- helpers ----------
function R = chol_with_jitter(S)
    [R, p] = chol(S, 'lower');
    if p == 0, return; end
    jitter = 1e-10;
    I = eye(size(S));
    for k = 1:8
        [R, p] = chol(S + jitter*I, 'lower');
        if p == 0, return; end
        jitter = jitter * 10;
    end
    error('Covariance not SPD even after jittering.');
end

%% ---------- discrete truth (noisy mod-sum) ----------
function [red_star, syn_star] = truth_discr_noisy_modsum(R, eps)
    % Exact PMF for (X1,X2,X3,Y) with noisy mod-sum target.
    if R <= 1
        red_star = 0; syn_star = 0; return;
    end
    p = zeros(R,R,R,R);                 % dims: X1 X2 X3 Y
    pN = ones(1,R) * (eps/(R-1));       % noise over {0..R-1}
    pN(1) = 1 - eps;                    % N=0 weight (index 1)
    invR3 = 1 / (R^3);
    for x1 = 1:R
        for x2 = 1:R
            for x3 = 1:R
                sum0 = mod((x1-1)+(x2-1)+(x3-1), R);
                for y = 1:R
                    n = mod((y-1) - sum0, R) + 1; % 1..R index for noise
                    p(x1,x2,x3,y) = invR3 * pN(n);
                end
            end
        end
    end
    % Red* = min I(Xi;Y); Syn* = I(X1,X2,X3;Y) - max I(Xi,Xj;Y)
    nsources = 3; total_vars = 4; % last dim is Y
    mi_values = zeros(1, nsources);
    for s = 1:nsources
        joint_xy = sum(p, setdiff(1:total_vars, [s total_vars]));
        mi_values(s) = mutualInformationLast(joint_xy);
    end
    red_star = min(mi_values);
    info_all = mutualInformationLast(p);
    pairwise_mis = [];
    for a = 1:nsources
        for b = a+1:nsources
            joint_xxy = sum(p, setdiff(1:total_vars, [a b total_vars]));
            pairwise_mis(end+1) = mutualInformationLast(joint_xxy); %#ok<AGROW>
        end
    end
    syn_star = info_all - max(pairwise_mis);
end

function [hFill,hLine] = plotMeanBand(x, y, sem, color, lineSpec, lw, ms, alphaFace)
    % Plots a mean line with a shaded band of ±3*SEM.
    % x,y,sem are vectors of same length. lineSpec e.g. '-', ':s', '-.d'
    if isrow(x),  x = x(:);  end
    if isrow(y),  y = y(:);  end
    if isrow(sem), sem = sem(:); end

    up = y + 3*sem;
    lo = y - 3*sem;

    X = [x; flipud(x)];
    Y = [up; flipud(lo)];

    holdState = ishold(gca);
    hold on;
    hFill = fill(X, Y, color, 'FaceAlpha', alphaFace, 'EdgeColor', 'none', 'HandleVisibility','off');
    hLine = plot(x, y, lineSpec, 'Color', color, 'LineWidth', lw, ...
                 'MarkerSize', ms, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', color);
    if ~holdState, hold off; end
end
