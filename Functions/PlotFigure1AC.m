function PlotFigure1AC()
    %% Discrete example
    clc, clear;
    rng(20);
    results_folder = fullfile('Results'); % '..',
    figures_folder = fullfile('Figures'); % '..',
    
    n_simulations = 10000;
    MIJoint_values  = zeros(2,n_simulations);
    MISingle_values = zeros(2,2,n_simulations);
    PID_values    = zeros(2,4,n_simulations);
    P_dist_all = zeros(2, 2, 4, 4);
    SYnsh_values = zeros(2,n_simulations);
    max_marg = 0;
    trialsperStim_all = [50, 500];%  [50, 500];
    marg = {1,trialsperStim_all};
    for trialIdx = 1:length(trialsperStim_all)
        trialsperStim = trialsperStim_all(trialIdx);
        nstim = 2;
        nR =4;
        MI_opts.bias = 'plugin';
        MI_opts.bin_method = {'none'};
        MI_opts.suppressWarnings = true;
    
        PID_opts.bias='plugin';
        PID_opts.bin_method = {'none'};
        PID_opts.suppressWarnings = true;
    
        p_distr = zeros(2,nstim,nR,nR);
        fprintf('Calculating values for %d trials per stim\n', trialsperStim)
        neuron1 = randi(nR, 1, nstim*trialsperStim);
        neuron2 = randi(nR, 1, nstim*trialsperStim);
        stimuli = [];
        for idx_stim = 1:nstim
            stimuli = [stimuli, idx_stim*ones(1, trialsperStim)];
        end
    
        tmp_pdf = prob_estimator({stimuli, neuron1, neuron2}, {'P(A,B,C)'});
        P_dist_all(trialIdx, :, :, :, :) =  cell2mat(tmp_pdf);
        stim_cond1 = squeeze(tmp_pdf{1}(1,:,:))*2/1;
        stim_cond2 = squeeze(tmp_pdf{1}(2,:,:))*2/1;
    
        marg11 = sum(stim_cond1,1);
        marg12 = sum(stim_cond1,2);
        marg21 = sum(stim_cond2,1);
        marg22 = sum(stim_cond2,2);
    
        margtrial= {marg11 marg21;
            marg12 marg22};
    
        marg{trialIdx} = margtrial;
    
        max_margtmp = max([marg11; marg12'; marg21; marg22'], [], 'all');
        max_marg = max([max_marg, max_margtmp]);
    
        parfor i = 1:n_simulations
            neuron1 = randi(nR, 1, nstim*trialsperStim);
            neuron2 = randi(nR, 1, nstim*trialsperStim);
            inf_tm = MI({[neuron1;neuron2], stimuli}, {'I(A;B)'}, MI_opts);
            MIJoint_values(trialIdx,i) = inf_tm{1};
            inf_tm = MI({neuron1, stimuli}, {'I(A;B)'}, MI_opts);
            misingle1=inf_tm{1};
            inf_tm = MI({neuron2, stimuli}, {'I(A;B)'}, MI_opts);
            misingle2=inf_tm{1};
            MISingle_values(trialIdx,:,i) = [misingle1 misingle2];
            PID_tmp = PID({neuron1, neuron2, stimuli}, {'PID_atoms'}, PID_opts);
            PID_values(trialIdx,:,i) = [PID_tmp{2}, PID_tmp{3}, PID_tmp{4}, PID_tmp{1}];
        end
    end
    save(fullfile(results_folder, 'results_Figure_1A.mat'))
    
    %%
    load(fullfile(results_folder, 'results_Figure_1A.mat'))
    min_size=2;
    ratio = 4;
    col_width = min_size+ratio;
    set(0, 'DefaultTextFontName', 'Arial');
    set(0, 'DefaultAxesFontName', 'Arial');
    set(0, 'DefaultUicontrolFontName', 'Arial');
    set(0, 'DefaultUitableFontName', 'Arial');
    set(0, 'DefaultUipanelFontName', 'Arial');
    set(0, 'DefaultLegendFontName', 'Arial');
    
    figure_handle = figure('Units', 'inches', 'Position', [1, 1, 20, 12.5]);
    fontSize = 20;
    t = tiledlayout(figure_handle, 2*col_width,3*col_width,"TileSpacing","tight");
    row_length = 3*col_width;
    row_heigth = col_width;
    
    
    for n_trials_idx = 1:2
        position = col_width*2+1+row_heigth*row_length*(n_trials_idx-1);
        axes_hist(n_trials_idx) = nexttile(position, [row_heigth col_width]);
    
        % histogram(MI_values(n_trials_idx,:));
    
        [fjoint,xfjoint] = histcounts(MIJoint_values(n_trials_idx,:)); %kde(MIJoint_values(n_trials_idx,:));
        mi_sing = MISingle_values(n_trials_idx,:,:);
        [fsingle,xfsingle] = histcounts(mi_sing(:)); %kde(mi_sing(:));
        [fsyn,xfsyn] = histcounts(PID_values(n_trials_idx, 4, :)); %kde(PID_values(n_trials_idx, 4, :));
        [fred,xfred] = histcounts(PID_values(n_trials_idx, 1, :)); %kde(PID_values(n_trials_idx, 4, :));
        pid_unq = PID_values(n_trials_idx, 2:3, :);
        [funq,xfunq] = histcounts(pid_unq(:)); %kde(pid_unq);
    
        hold on;
        color_list = [92.9, 70.2, 14.5; % joint
            47.8, 67.8, 20.4; % syn
            50.2, 50.2, 50.2;
            0, 0.4470, 0.7410;       % Color for 'Red'
            ];% unq
        lw=1.7;
        plot((xfjoint(2:end)+xfjoint(1:end-1))/2, fjoint/sum(fjoint, 'all'),'DisplayName','joint','Color',color_list(1,:)./100,"LineWidth",lw)
        plot((xfsingle(2:end)+xfsingle(1:end-1))/2, fsingle/sum(fsingle, 'all'),'DisplayName','single',"LineWidth",lw)
        plot((xfsyn(2:end)+xfsyn(1:end-1))/2, fsyn/sum(fsyn, 'all'),'DisplayName','syn','Color',color_list(2,:)./100,"LineWidth",lw)
        plot((xfred(2:end)+xfred(1:end-1))/2, fred/sum(fred, 'all'),'DisplayName','red','Color',color_list(4,:),"LineWidth",lw)
        plot((xfunq(2:end)+xfunq(1:end-1))/2, funq/sum(funq, 'all'),'DisplayName','unq','Color',color_list(3,:)./100,"LineWidth",lw)
    
        xlim([-0.01, 0.4])
        legend({'Joint','Single', 'Syn', 'Red','Unq'}, 'Box', 'off',  'FontSize', fontSize)
        if n_trials_idx==2
            xlabel('Information [bits]', 'FontSize', fontSize+3)
        end
        ax = gca;
        ax.FontSize = fontSize;
        ax.XLabel.FontSize = fontSize;
        ax.YLabel.FontSize = fontSize;
    end
    
    
    linkaxes(axes_hist,'x');
    
    fileName = 'Figure_1Ahistograms.svg';
    filePath = fullfile(figures_folder, fileName);
    saveas(t, filePath, 'svg')
    %%
    
    yellowColors = [ ...
        255, 255, 229;
        255, 247, 188;
        254, 227, 145;
        254, 196, 79;
        254, 153, 41;
        217, 95, 14;
        153, 52, 4
        ] / 255;
    
    purpleColors = [ ...
        252, 251, 253;
        239, 237, 245;
        218, 218, 235;
        188, 189, 220;
        158, 154, 200;
        128, 125, 186;
        106, 81, 163;
        84, 39, 143;
        63, 0, 125
    ] / 255;
    
    numColors = 256;
    yellowColormap = interp1(linspace(0, 1, size(yellowColors, 1)), yellowColors, linspace(0, 1, numColors));
    purpleColormap = interp1(linspace(0, 1, size(purpleColors, 1)), purpleColors, linspace(0, 1, numColors));
    
    min_size=2;
    ratio = 4;
    col_width = min_size+ratio;
    set(0, 'DefaultTextFontName', 'Arial');
    set(0, 'DefaultAxesFontName', 'Arial');
    set(0, 'DefaultUicontrolFontName', 'Arial');
    set(0, 'DefaultUitableFontName', 'Arial');
    set(0, 'DefaultUipanelFontName', 'Arial');
    set(0, 'DefaultLegendFontName', 'Arial');
    
    close all
    figure_handle = figure('Units', 'inches', 'Position', [1, 1, 20, 12.5]);
    
    t = tiledlayout(figure_handle, 2*col_width,3*col_width,"TileSpacing","tight");
    row_length = 3*col_width;
    row_heigth = col_width;
    
    max_val_prob = max(P_dist_all, [],'all');
    max_max_all = max([max_val_prob max_marg]);
    k=0;
    fontSize = 20;
    for n_trials_idx = 1:2
        for stim_idx = 1:2
            position = row_length*row_heigth*(n_trials_idx-1)+col_width*(stim_idx-1)+min_size+1;
            ax(n_trials_idx,stim_idx) = nexttile(position, [ratio ratio]); 
            h=heatmap(squeeze(P_dist_all(n_trials_idx,stim_idx,:,:)));
            colormap(gca, yellowColormap);
            s=struct(h);
            s.XAxis.Visible='off';
            s.YAxis.Visible='off';
            clim([0, 0.07]);
            h.FontSize = fontSize;
            h.CellLabelColor = 'none';
            if n_trials_idx==1
                title_text = sprintf('Stimulus {%d}', stim_idx);
                h.Title = title_text;
                h.FontSize = fontSize+2;
            end
            xlabel('R1')
            ylabel('R2')
        end
    end
    
    for n_trials_idx = 1:2
        for stim_idx = 1:2
            position = row_length*row_heigth*(n_trials_idx-1)+ row_length*ratio+(stim_idx-1)*col_width+min_size+1;
            ax(n_trials_idx,stim_idx) = nexttile(position, [min_size ratio]); 
            marg_tmp = marg{n_trials_idx};
            % h(k)=bar(1:nR,marg_tmp{1,stim_idx});
            h=heatmap(marg_tmp{1,stim_idx});
            % axis equal
            % axis off
            colormap(gca, purpleColormap)
            colorbar off
    
            s=struct(h);
            s.YAxis.Visible='off';
            clim([0, max_marg]);
            h.FontSize = fontSize;
            h.CellLabelColor = 'none';
            xlabel('R1')
        end
    end
    
    
    
    
    for n_trials_idx = 1:2
        for stim_idx = 1:2
            position = row_length*row_heigth*(n_trials_idx-1)+(stim_idx-1)*col_width+1;
            ax(n_trials_idx,stim_idx) = nexttile(position, [ratio min_size]); 
            marg_tmp = marg{n_trials_idx};
            h=heatmap(marg_tmp{2,stim_idx});
            % axis equal
            % axis off
            colormap(gca, purpleColormap)
            s=struct(h);
            s.XAxis.Visible='off';
            h.FontSize = fontSize;
            h.CellLabelColor = 'none';
            clim([0, max_marg]);
            ylabel('R2')
            % h(k)=bar(1:nR,marg_tmp{2,stim_idx});
            % camroll(90)
        end
    end
    fileName = 'Figure_1PanelA.svg';
    filePath = fullfile(figures_folder, fileName);
    saveas(gca, filePath, 'svg');
    
    
    
    %% Gaussian example
    n_trials = 50;
    n_simulations = 10000;
    d = 5;
    
    MI_joint = zeros(1, n_simulations);
    MI_single = zeros(2, n_simulations);
    
    for sim = 1:n_simulations
        % Generate 5D Gaussian stimulus
        S = randn(n_trials, d);  % Covariance is identity
    
        % Linear transformation matrices
        A1 = randn(d);  % R1 = A1*S + noise
        A2 = randn(d);  % R2 = A2*S + noise
    
        % Noise (uncorrelated Gaussian)
        eps1 = randn(n_trials, d);
        eps2 = randn(n_trials, d);
    
        R1 = eps1;
        R2 = eps2;
    
        % Joint data matrices
        data_joint = [S, R1, R2];
        data_r1 = [S, R1];
        data_r2 = [S, R2];
    
        % Covariance matrices
        cov_joint = cov(data_joint);
        cov_r1 = cov(data_r1);
        cov_r2 = cov(data_r2);
        cov_s = cov(S);
    
        % Compute MI(S; R1, R2)
        Sigma_S = cov_s;
        Sigma_R = cov([R1, R2]);
        Sigma_SR = cov_joint;
        MI_joint(sim) = 0.5 * log(det(Sigma_S) * det(Sigma_R) / det(Sigma_SR));
    
        % MI(S; R1)
        Sigma_R1 = cov(R1);
        Sigma_SR1 = cov_r1;
        MI_single(1, sim) = 0.5 * log(det(Sigma_S) * det(Sigma_R1) / det(Sigma_SR1));
    
        % MI(S; R2)
        Sigma_R2 = cov(R2);
        Sigma_SR2 = cov_r2;
        MI_single(2, sim) = 0.5 * log(det(Sigma_S) * det(Sigma_R2) / det(Sigma_SR2));
    end
    
    % Number of trials
    % n_trials = 128;
    short_trials = 20;
    
    % Assume last simulation's data
    S_all = S; R1_all = R1; R2_all = R2;
    
    % Estimate full covariance matrix from all trials
    data_all = [S_all, R1_all, R2_all];
    cov_all = cov(data_all);
    
    % Estimate covariance from only first 60 trials
    data_short = data_all(1:short_trials, :);
    cov_short = cov(data_short);
    
    %%
    % Variable block size
    block_size = d;
    
    % Total dimensions
    dim_labels = {'S1','S2','S3','S4','S5', ...
                  'R1_1','R1_2','R1_3','R1_4','R1_5', ...
                  'R2_1','R2_2','R2_3','R2_4','R2_5'};
    
    % Define yellow-orange colormap
    yellowColors = [ ...
        255, 247, 188;
        % 254, 227, 145;
        % 254, 196, 79;
        254, 196, 79;
        254, 153, 41;
        217,  95, 14;
        153,  52,  4
    ] / 255;
    
    numColors = 256;
    yellowColormap = interp1(linspace(0, 1, size(yellowColors, 1)), ...
                             yellowColors, ...
                             linspace(0, 1, numColors));
    
    % Determine common color scale
    cmin = min([cov_short(:); cov_all(:)]);
    cmax = max([cov_short(:); cov_all(:)]);
    
    % Plot using tiled layout
    figure;
    t = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    ax1 = nexttile;
    imagesc(cov_short);
    title(sprintf('Covariance Matrix (N = %d)', short_trials));
    axis square;
    caxis([cmin, cmax]);
    
    ax2 = nexttile;
    imagesc(cov_all);
    title(sprintf('Covariance Matrix (N = %d)', n_trials));
    axis square;
    caxis([cmin, cmax]);
    % axes = [ax1, ax2];
    % Apply formatting to both subplots
    for ax = [ax1, ax2]
        axes(ax); %#ok<LAXES>
        colormap(yellowColormap);
    
        % Dividing lines
        hold on;
        for k = block_size:block_size:(3*block_size)-1
            xline(k+0.5, 'k', 'LineWidth', 1.2);
            yline(k+0.5, 'k', 'LineWidth', 1.2);
        end
        hold off;
    
        % Only show middle block label
        tick_positions = (block_size+1)/2 + [0, block_size, 2*block_size];  % [3, 8, 13] for block size 5
        tick_labels = {'S', 'R1', 'R2'};
    
        xticks(tick_positions);
        yticks(tick_positions);
        xticklabels(tick_labels);
        yticklabels(tick_labels);
        xtickangle(0);
    end
    
    % Add one shared colorbar
    cb = colorbar;
    cb.Layout.Tile = 'east';
    
    cb.Label.String = 'Covariance';
    cb.Label.FontSize = 12;
    
    % Save
    saveas(gcf, fullfile(figures_folder, 'Figure1C.svg'));
end
    %%

function MIX = mutualInformationXYZ(pxyz)
% Calculate marginal probabilities
px = squeeze(sum(pxyz, [2, 3])); % Marginal probability density of Z
pyz = squeeze(sum(pxyz, 1)); % Joint probability density of X and Y

% Initialize mutual information
MIX = 0;

% Loop through each value of x, y, and z
for i = 1:size(pxyz, 1)
    for j = 1:size(pxyz, 2)
        for k = 1:size(pxyz, 3)
            if pxyz(i,j,k) ~= 0 && px(i) ~= 0 && pyz(j,k) ~= 0
                MIX = MIX + pxyz(i,j,k) * log2(pxyz(i,j,k) / (px(i) * pyz(j,k)));
            end
        end
    end
end
end