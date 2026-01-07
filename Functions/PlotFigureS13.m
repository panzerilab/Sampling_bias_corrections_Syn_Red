results_folder = fullfile('..','Results');
figures_folder = fullfile('..','Figures');
% Parameter setup
dims = 1:42;
trial_sizes = 2*64 * 2.^(0:5);  % [64, 128, 256, 512, 1024, 2048]
n_simulations = 100;

% Preallocate for mean and std of MI values
MI_joint_raw = zeros(length(dims), length(trial_sizes), n_simulations);
mean_MI_joint = zeros(length(dims), length(trial_sizes));
std_MI_joint = zeros(length(dims), length(trial_sizes));

MI_single_raw = zeros(2, length(dims), length(trial_sizes), n_simulations);
mean_MI_single = zeros(2, length(dims), length(trial_sizes));
std_MI_single = zeros(2, length(dims), length(trial_sizes));

entropy_raw  = zeros(length(dims), length(trial_sizes), n_simulations);
mean_entropy = zeros(length(dims), length(trial_sizes));
std_entropy  = zeros(length(dims), length(trial_sizes));

log2 = log(2);

% Loop over dimensionalities and trial sizes
for i_d = 1:length(dims)
    d = dims(i_d);

    for i_t = 1:length(trial_sizes)
        n_trials = trial_sizes(i_t);

        if ~((d == 10 || d == 20) || (n_trials == 256 || n_trials == 128))
            continue;
        end

        % Temporary storage for MI values for this combo
        MI_joint = zeros(1, n_simulations);
        MI_single = zeros(2, n_simulations);
        entropy_tmp = zeros(1, n_simulations);

        for sim = 1:n_simulations
            % Generate stimulus
            S = randn(n_trials, d);  % Gaussian, identity covariance

            % Linear transformations (optional; used only if R has signal)
            A1 = randn(d);
            A2 = randn(d);

            % Noise
            eps1 = randn(n_trials, d);
            eps2 = randn(n_trials, d);

            % Responses (add signal here if desired)
            R1 =eps1; %  S * A1' + 
            R2 =eps2; %  S * A2' + 

            % Covariances
            Sigma_S = cov(S);
            Sigma_R = cov([R1, R2]);
            Sigma_SR = cov([S, R1, R2]);

            % Compute MI(S; R1, R2)
            MI_joint(sim)     = (1 / (2 * log2)) * log(det(Sigma_S) * det(Sigma_R) / det(Sigma_SR));
            entropy_tmp(sim)     = 0.5 / log(2) * log(det(Sigma_S));

            % MI(S; R1)
            Sigma_R1 = cov(R1);
            Sigma_SR1 = cov([S, R1]);
            MI_single(1, sim) = (1 / (2 * log2)) * log(det(Sigma_S) * det(Sigma_R1) / det(Sigma_SR1));

            % MI(S; R2)
            Sigma_R2 = cov(R2);
            Sigma_SR2 = cov([S, R2]);
            MI_single(2, sim) = (1 / (2 * log2)) * log(det(Sigma_S) * det(Sigma_R2) / det(Sigma_SR2));
        end

        % Store statistics
        entropy_raw(i_d, i_t, :) = entropy_tmp;
        mean_entropy(i_d, i_t) = mean(entropy_tmp);
        std_entropy(i_d, i_t) = std(entropy_tmp);

        MI_joint_raw(i_d, i_t, :) = MI_joint;
        mean_MI_joint(i_d, i_t) = mean(MI_joint);
        std_MI_joint(i_d, i_t) = std(MI_joint);

        MI_single_raw(:,i_d, i_t, :) = MI_single;
        mean_MI_single(:, i_d, i_t) = mean(MI_single, 2);
        std_MI_single(:, i_d, i_t) = std(MI_single, 0, 2);
    end

    %fprintf('Finished d = %d/%d\n', i_d, length(dims));
end

%% Numerical information


% Preallocate result matrices: dims x ntrials
Joint_goodman_all = zeros(length(dims), length(trial_sizes));
Single_goodman_all = zeros(length(dims), length(trial_sizes));
entropy_goodman_all = zeros(length(dims), length(trial_sizes));

Joint_cai_all = zeros(length(dims), length(trial_sizes));
Single_cai_all = zeros(length(dims), length(trial_sizes));
entropy_cai_all = zeros(length(dims), length(trial_sizes));

% Loop over all dimensions and trial sizes
for i_d = 1:length(dims)
    M = dims(i_d);
    ks  = 1:M;
    ks2 = 1:(2*M);
    ks3 = 1:(3*M);

    for i_t = 1:length(trial_sizes)
        ntrials = trial_sizes(i_t);
        if ~((M == 10 || M == 20  ||M ==30) || (ntrials == 256 || ntrials == 128))
            continue;
        end
        %% --- Goodman Bias ---
        biasH_M_goodman   = M * log(2 / (ntrials - 1)) + sum(psi((ntrials - ks) / 2));
        biasH_XY_goodman  = 2 * M * log(2 / (ntrials - 1)) + sum(psi((ntrials - ks2) / 2));
        biasH_MXY_goodman = 3 * M * log(2 / (ntrials - 1)) + sum(psi((ntrials - ks3) / 2));

        Joint_goodman = (1 / (2 * log2)) * ...
            (biasH_M_goodman + biasH_XY_goodman - biasH_MXY_goodman);

        Single_goodman = (1 / (2 * log2)) * ...
            (2 * biasH_M_goodman - biasH_XY_goodman);

        %% --- Cai Bias ---
        biasH_M_cai   = sum(log(1 - ks / ntrials));
        biasH_XY_cai  = sum(log(1 - ks2 / ntrials));
        biasH_MXY_cai = sum(log(1 - ks3 / ntrials));

        Joint_cai = (1 / (2 * log2)) * ...
            (biasH_M_cai + biasH_XY_cai - biasH_MXY_cai);

        Single_cai = (1 / (2 * log2)) * ...
            (2 * biasH_M_cai - biasH_XY_cai);

        %% Store results
        entropy_goodman_all(i_d, i_t) = (1 / (2 * log2)) *biasH_M_goodman;
        Joint_goodman_all(i_d, i_t) = Joint_goodman;
        Single_goodman_all(i_d, i_t) = Single_goodman;

        entropy_cai_all(i_d, i_t) = (1 / (2 * log2)) *biasH_M_cai;
        Joint_cai_all(i_d, i_t) = Joint_cai;
        Single_cai_all(i_d, i_t) = Single_cai;
    end
end
%%
fig = figure('Units','centimeters','Position',[1 1 16 16]);
tl  = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
set(groot, 'defaultAxesFontSize', 9);
set(groot, 'defaultAxesTitleFontSizeMultiplier', 1.2);
set(groot, 'defaultAxesLabelFontSizeMultiplier', 1.2);

colors = lines(3);  % Use consistent colors across panels

% ========== Panel 1 (1,1): Entropy vs. Trial Size ==========
nexttile(1);
d_vals = [10, 20];
d_entropy = [10, 20];
for i = 1:length(d_entropy)
    d_idx = find(dims == d_entropy(i)); 
    label  = sprintf('Numerical (d=%d)',d_entropy(i));
    plot(trial_sizes, mean_entropy(d_idx,:), 'LineWidth', 2, 'LineStyle', '-', 'Color', colors(i,:), 'DisplayName',label); hold on;
end
for i = 1:length(d_entropy)
    d_idx = find(dims == d_entropy(i));
    label  = sprintf('Goodman (d=%d)',d_entropy(i));
    plot(trial_sizes, entropy_goodman_all(d_idx,:), '--s', 'Color', colors(i,:), 'LineWidth', 1.5, 'DisplayName',label); hold on;
end
for i = 1:length(d_entropy)
    d_idx = find(dims == d_entropy(i));
    label  = sprintf('Cai (d=%d)',d_entropy(i));
    plot(trial_sizes, entropy_cai_all(d_idx,:), '-.^', 'Color', colors(i,:), 'LineWidth', 1.5, 'DisplayName',label);
end
for i = 1:length(d_entropy)
    d_idx = find(dims == d_entropy(i));
    mu = mean_entropy(d_idx,:);
    err = 2 * std_entropy(d_idx,:) / sqrt(n_simulations);
    fill([trial_sizes, fliplr(trial_sizes)], ...
         [mu - err, fliplr(mu + err)], ...
         colors(i,:), 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility','off'); hold on;
end
% for i = 1:length(d_vals)
%     d_idx = find(dims == d_vals(i));
%     mu_goodman = entropy_goodman_all(d_idx,:);
%     mu_cai = entropy_cai_all(d_idx,:);
%     % err = 2 * std_MI_joint(d_idx,:) / sqrt(n_simulations);
%     err = 2 * std_entropy(d_idx,:) / sqrt(n_simulations);
%     fill([trial_sizes, fliplr(trial_sizes)], ...
%          [mu_goodman - err, fliplr(mu_goodman + err)], ...
%          colors(i,:), 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility','off'); hold on;
%     fill([trial_sizes, fliplr(trial_sizes)], ...
%          [mu_cai - err, fliplr(mu_cai + err)], ...
%          colors(i,:), 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility','off'); hold on;
% end

%xlabel('Trial size'); 
ylabel('Entropy bias (bits)');
% ylim([-3, 0.01])
title('Entropy Bias vs. Trial Size');
legend()
% legend('Goodman (d=10)', 'Goodman (d=20)', 'Goodman (d=30)', ...% 'Numerical (d=10)', 'Numerical (d=20)', 'Numerical (d=30)',...
%     'Cai (d=10)', 'Cai (d=20)','Cai (d=30)', ...
%     Location='southeast'); 
grid off; box off; set(gca,'XScale','log');

% ========== Panel 2 (1,2): Entropy vs. Dimension ==========
nexttile(2);
n_vals = [128, 256];
for i = 1:length(n_vals)
    n_idx = find(trial_sizes == n_vals(i));
    label  = sprintf('Numerical (N=%d)',n_vals(i));
    plot(dims, mean_entropy(:,n_idx), 'LineWidth', 2, 'LineStyle', '-', 'Color', colors(i,:), 'DisplayName',label); hold on;
end
for i = 1:length(n_vals)
    n_idx = find(trial_sizes == n_vals(i));
    label  = sprintf('Goodman (N=%d)',n_vals(i));
    plot(dims, entropy_goodman_all(:,n_idx), '--s', 'Color', colors(i,:), 'LineWidth', 1.5, 'DisplayName',label); hold on;
end
for i = 1:length(n_vals)
    n_idx = find(trial_sizes == n_vals(i));
    label  = sprintf('Cai (N=%d)',n_vals(i));
    plot(dims, entropy_cai_all(:,n_idx), '-.^', 'Color', colors(i,:), 'LineWidth', 1.5, 'DisplayName',label);
end
for i = 1:length(n_vals)
    n_idx = find(trial_sizes == n_vals(i));
    mu = mean_entropy(:,n_idx);
    err = 2 * std_entropy(:,n_idx) / sqrt(n_simulations);
    fill([dims, fliplr(dims)], ...
         [mu' - err', fliplr(mu' + err')], ...
         colors(i,:), 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility','off'); hold on;
end
% for i = 1:length(n_vals)
%     n_idx = find(trial_sizes == n_vals(i));
%     % mu_goodman = mean_MI_joint(d_idx,:)-Joint_goodman_all(d_idx,:);
%     % mu_cai = mean_MI_joint(d_idx,:)-Joint_cai_all(d_idx,:);
%     mu_goodman = entropy_goodman_all(:,n_idx);
%     mu_cai = entropy_cai_all(:,n_idx);
%     % err = 2 * std_MI_joint(d_idx,:) / sqrt(n_simulations);
%     err = 2 * std_entropy(:,n_idx) / sqrt(n_simulations);
%     fill([dims, fliplr(dims)], ...
%          [mu_goodman' - err', fliplr(mu_goodman' + err')], ...
%          colors(i,:), 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility','off'); hold on;
%     fill([dims, fliplr(dims)], ...
%          [mu_cai' - err', fliplr(mu_cai' + err')], ...
%          colors(i,:), 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility','off'); hold on;
% end

%xlabel('Dimension'); % ylabel('Entropy Bias (bits)');
xlim([2, 42]);
% ylim([-6, 0.01])
title('Entropy Bias vs. Dimension'); 
legend(Location='southwest');  
grid off; box off;
nexttile(3);
for i = 1:length(d_entropy)
    d_idx = find(dims == d_entropy(i));
    % plot(trial_sizes, mean_MI_joint(d_idx,:), '-', 'Color', colors(i,:), 'LineWidth', 2); hold on;
    % plot(trial_sizes, Joint_goodman_all(d_idx,:), '--', 'Color', colors(i,:), 'LineWidth', 1.5);
    % plot(trial_sizes, Joint_cai_all(d_idx,:), ':', 'Color', colors(i,:), 'LineWidth', 1.5);
    label  = sprintf('Goodman (d=%d)',d_entropy(i));
    plot(trial_sizes, mean_MI_joint(d_idx,:)-Joint_goodman_all(d_idx,:), '--s', 'Color', colors(i,:), 'LineWidth', 1.5, 'DisplayName',label); hold on;
    % plot(trial_sizes, squeeze(mean_MI_single(1,d_idx,:))'-Single_goodman_all(d_idx,:), '--s', 'Color', colors(i,:), 'LineWidth', 1.5); hold on;

end
for i = 1:length(d_entropy)
    d_idx = find(dims == d_entropy(i));
    % plot(trial_sizes, mean_MI_joint(d_idx,:), '-', 'Color', colors(i,:), 'LineWidth', 2); hold on;
    % plot(trial_sizes, Joint_goodman_all(d_idx,:), '--', 'Color', colors(i,:), 'LineWidth', 1.5);
    % plot(trial_sizes, Joint_cai_all(d_idx,:), ':', 'Color', colors(i,:), 'LineWidth', 1.5);
    label  = sprintf('Cai (d=%d)',d_entropy(i));
    plot(trial_sizes,  mean_MI_joint(d_idx,:)-Joint_cai_all(d_idx,:),    '-.^', 'Color', colors(i,:), 'LineWidth', 1.5, 'DisplayName',label); hold on;
    % plot(trial_sizes,  squeeze(mean_MI_single(1,d_idx,:))'-Single_cai_all(d_idx,:),    '-.^', 'Color', colors(i,:), 'LineWidth', 1.5);
end
for i = 1:length(d_entropy)
    d_idx = find(dims == d_entropy(i));
    mu_goodman = mean_MI_joint(d_idx,:)-Joint_goodman_all(d_idx,:);
    mu_cai = mean_MI_joint(d_idx,:)-Joint_cai_all(d_idx,:);
    % mu_goodman = squeeze(mean_MI_single(1,d_idx,:))'-Single_goodman_all(d_idx,:);
    % mu_cai = squeeze(mean_MI_single(1,d_idx,:))'-Single_cai_all(d_idx,:);
    err = 2 * std_MI_joint(d_idx,:) / sqrt(n_simulations);
    % err = 2 * squeeze(std_MI_single(1,d_idx,:))' / sqrt(n_simulations);
    fill([trial_sizes, fliplr(trial_sizes)], ...
         [mu_goodman - err, fliplr(mu_goodman + err)], ...
         colors(i,:), 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility','off'); hold on;
    fill([trial_sizes, fliplr(trial_sizes)], ...
         [mu_cai - err, fliplr(mu_cai + err)], ...
         colors(i,:), 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility','off'); hold on;
end
xlabel('Trial size'); ylabel('Residual MI Bias (bits)');
title('Residual MI Bias vs. Trial Size');
legend(Location='southeast'); 
grid off; box off; set(gca,'XScale','log');

% ========== Panel 6 (3,2): Joint MI vs. Dimension ==========
nexttile(4);
for i = 1:length(n_vals)
    n_idx = find(trial_sizes == n_vals(i));
    % plot(dims, mean_MI_joint(:,n_idx), '-', 'Color', colors(i,:), 'LineWidth', 2); hold on;
    label  = sprintf('Goodman (N=%d)',n_vals(i));
    plot(dims, mean_MI_joint(:,n_idx)-Joint_goodman_all(:,n_idx), '--s', 'Color', colors(i,:), 'LineWidth', 1.5, 'DisplayName',label); hold on;
    % plot(dims, squeeze(mean_MI_single(1,:,n_idx))'-Single_goodman_all(:,n_idx), '--s', 'Color', colors(i,:), 'LineWidth', 1.5); hold on;
end
for i = 1:length(n_vals)
    n_idx = find(trial_sizes == n_vals(i));
    label  = sprintf('Cai (N=%d)',n_vals(i));
    % plot(dims, mean_MI_joint(:,n_idx), '-', 'Color', colors(i,:), 'LineWidth', 2); hold on;
    plot(dims,  mean_MI_joint(:,n_idx)-Joint_cai_all(:,n_idx), '-.^', 'Color', colors(i,:), 'LineWidth', 1.5, 'DisplayName',label);
    % plot(dims,  squeeze(mean_MI_single(1,:,n_idx))'-Single_cai_all(:,n_idx), '-.^', 'Color', colors(i,:), 'LineWidth', 1.5);
end
for i = 1:length(n_vals)
    n_idx = find(trial_sizes == n_vals(i));
    mu_goodman = mean_MI_joint(:,n_idx)-Joint_goodman_all(:,n_idx);
    mu_cai = mean_MI_joint(:,n_idx)-Joint_cai_all(:,n_idx);
    % mu_goodman = squeeze(mean_MI_single(1,:,n_idx))'-Single_goodman_all(:,n_idx);
    % mu_cai = squeeze(mean_MI_single(1,:,n_idx))'-Single_cai_all(:,n_idx);
    err = 2 * std_MI_joint(:,n_idx) / sqrt(n_simulations);
    % err = 2 * squeeze(std_MI_single(1,:,n_idx))' / sqrt(n_simulations);
    fill([dims, fliplr(dims)], ...
         [mu_goodman' - err', fliplr(mu_goodman' + err')], ...
         colors(i,:), 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility','off'); hold on;
    fill([dims, fliplr(dims)], ...
         [mu_cai' - err', fliplr(mu_cai' + err')], ...
         colors(i,:), 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility','off'); hold on;
end
xlabel('Dimension'); %ylabel('Residual MI Bias (bits)')
xlim([2, 42])
title('Residual MI Bias vs. Dimension');
legend(Location='northwest');
grid off; box off;
outDir = fullfile(fileparts(pwd), 'Figures_mat');  

if ~exist(outDir, 'dir')
    mkdir(outDir);
end

exportgraphics(gcf, fullfile(outDir, 'Figure_S13.svg'), 'ContentType','vector');





