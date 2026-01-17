function runSimulation_venkatesh(nRbins, trial_categories, numRep, simul_case, alpha, info_amount, isalphasweep, isBinSweep, redundancyMeasure)
n_neurons = 2;

% Compute Rates
switch simul_case
    case 'uncorr_unique'
        B = 5; beta = 0; gamma_flat = 0;
    case 'unique_high_red'
        B = 5; beta = 0.4; gamma_flat = 2;
    case 'unique_high_syn'
        B = 5; beta = 0.1; gamma_flat = 20;
    case 'bit_of_all'
        B = 5; beta = 0.25; gamma_flat = 17.5;
end

opts_PID.bin_method = {'none'};
opts_PID.suppressWarnings = true;
opts_PID.computeNulldist = false;
opts_PID.pid_constrained = true;
opts_PID.shuff = 20;
opts_PID.xtrp = 20;
opts_PID.redundancy_measure = redundancyMeasure;

opts_PID_plugin = opts_PID;
opts_PID_plugin.bias = 'plugin';

opts_MI.bin_method = {'none'};
opts_MI.suppressWarnings = true;
opts_MI.computeNulldist = false;
opts_MI.bias = 'plugin';

opts_MI_qe = opts_MI;
opts_MI_qe.bias = 'qe';

if strcmp(redundancyMeasure, 'I_BROJA')
    reqOutputs_all = {'Joint', 'PID_atoms', 'Union', 'q_dist', 'p_ind'};
else
    reqOutputs_all = {'Joint', 'PID_atoms', 'Union'};
end

reqOutputs = {'Joint', 'PID_atoms', 'Union'};

if isBinSweep
    PID_v_venka = zeros(length(trial_categories), 8, numRep, length(nRbins));

    for binIdx = 1:length(nRbins)
        curr_nRbins = nRbins(binIdx);

        eps = 1;
        rate1 = B + alpha*beta*[1-eps,1-eps,eps,eps] + alpha*(1-beta)*[1-eps,eps,1-eps,eps];
        rate2 = B + alpha*beta*[1-eps,eps,1-eps,eps] + alpha*(1-beta)*[1-eps,1-eps,eps,eps];
        rateShared = gamma_flat*[1,1,1,1];

        num_trials = 2048;
        X1 = [poissrnd(rate1(1),1,num_trials), poissrnd(rate1(2),1,num_trials), ...
              poissrnd(rate1(3),1,num_trials), poissrnd(rate1(4),1,num_trials)];
        X2 = [poissrnd(rate2(1),1,num_trials), poissrnd(rate2(2),1,num_trials), ...
              poissrnd(rate2(3),1,num_trials), poissrnd(rate2(4),1,num_trials)];
        Xshared = [poissrnd(rateShared(1),1,num_trials), poissrnd(rateShared(2),1,num_trials), ...
                   poissrnd(rateShared(3),1,num_trials), poissrnd(rateShared(4),1,num_trials)];
        X1 = X1 + Xshared;
        X2 = X2 + Xshared;

        bin_opts.bin_method = {'eqpop'};
        bin_opts.n_bins = {curr_nRbins};
        [~, edges_spikes1] = binning({X1}, bin_opts);
        [~, edges_spikes2] = binning({X2}, bin_opts);
        edges_spikes1(1) = 0; edges_spikes2(1) = 0;
        edges_spikes1(end+1) = Inf; edges_spikes2(end+1) = Inf;

        parfor rep = 1:numRep
            PID_tmp = zeros(length(trial_categories), 8);
            for trials = 1:length(trial_categories)
                num_trials = trial_categories(trials);
                X1 = [poissrnd(rate1(1),1,num_trials), poissrnd(rate1(2),1,num_trials), ...
                      poissrnd(rate1(3),1,num_trials), poissrnd(rate1(4),1,num_trials)];
                X2 = [poissrnd(rate2(1),1,num_trials), poissrnd(rate2(2),1,num_trials), ...
                      poissrnd(rate2(3),1,num_trials), poissrnd(rate2(4),1,num_trials)];
                Xshared = [poissrnd(rateShared(1),1,num_trials), poissrnd(rateShared(2),1,num_trials), ...
                           poissrnd(rateShared(3),1,num_trials), poissrnd(rateShared(4),1,num_trials)];
                X1 = X1 + Xshared;
                X2 = X2 + Xshared;

                Y = [ones(1,num_trials), 2*ones(1,num_trials), 3*ones(1,num_trials), 4*ones(1,num_trials)];
                Yrnd = hShuffle({Y}, {'A'});

                X1_binned = discretize(X1, edges_spikes1);
                X2_binned = discretize(X2, edges_spikes2);
                input = {X1_binned, X2_binned, Y};

                outputs = MI({[X1_binned; X2_binned], Y}, {'I(A;B)'}, opts_MI);
                Ijoint = outputs{1};

                outputs = MI({[X1_binned; X2_binned], Y}, {'I(A;B)'}, opts_MI_qe);
                Ijoint_qe = max([0, outputs{1}]);
                debias_factor = Ijoint_qe / Ijoint;

                i1 = MI({X1_binned, Y}, {'I(A;B)'}, opts_MI_qe);
                i2 = MI({X2_binned, Y}, {'I(A;B)'}, opts_MI_qe);

                imx_qe = max([0, i1{1}]);
                imy_qe = max([0, i2{1}]);

                Ijoint_qe = max([Ijoint_qe, imx_qe, imy_qe]);

                PID_res_plugin = PID(input, reqOutputs_all, opts_PID_plugin);

                union_info = (PID_res_plugin{3} + PID_res_plugin{4} + PID_res_plugin{5}) * debias_factor;
                union_info = max([union_info, imx_qe, imy_qe]);
                union_info = min([union_info, imx_qe + imy_qe, Ijoint_qe]);

                PID_res_plugin{1} = Ijoint_qe;
                PID_res_plugin{2} = Ijoint_qe - union_info;
                PID_res_plugin{3} = imx_qe + imy_qe - union_info;
                PID_res_plugin{4} = union_info - imy_qe;
                PID_res_plugin{5} = union_info - imx_qe;
                PID_res_plugin{6} = union_info;

                PID_tmp(trials, :) = [PID_res_plugin{1}, PID_res_plugin{2}, PID_res_plugin{3}, ...
                                      PID_res_plugin{4}, PID_res_plugin{5}, PID_res_plugin{6}, 0, 0];
            end
            PID_v_venka(:,:,rep,binIdx) = PID_tmp;
        end
    end

    folder_name = 'BinSweep_venka';
    if ~isfolder(['Results/', folder_name])
        mkdir(['Results/', folder_name]);
    end
    save(['Results/', folder_name, '/', 'Simuldata_', char(simul_case), '_', char(info_amount), '_', ...
        char(opts_PID.redundancy_measure), '.mat'], 'PID_v_venka');
    return;
end

% Default (non-sweep) case
eps = 1;
rate1 = B + alpha*beta*[1-eps,1-eps,eps,eps] + alpha*(1-beta)*[1-eps,eps,1-eps,eps];
rate2 = B + alpha*beta*[1-eps,eps,1-eps,eps] + alpha*(1-beta)*[1-eps,1-eps,eps,eps];
rateShared = gamma_flat*[1,1,1,1];

num_trials = 2048;
X1 = [poissrnd(rate1(1),1,num_trials), poissrnd(rate1(2),1,num_trials), ...
      poissrnd(rate1(3),1,num_trials), poissrnd(rate1(4),1,num_trials)];
X2 = [poissrnd(rate2(1),1,num_trials), poissrnd(rate2(2),1,num_trials), ...
      poissrnd(rate2(3),1,num_trials), poissrnd(rate2(4),1,num_trials)];
Xshared = [poissrnd(rateShared(1),1,num_trials), poissrnd(rateShared(2),1,num_trials), ...
           poissrnd(rateShared(3),1,num_trials), poissrnd(rateShared(4),1,num_trials)];
X1 = X1 + Xshared;
X2 = X2 + Xshared;

bin_opts.bin_method = {'eqpop'};
bin_opts.n_bins = {nRbins};
[~, edges_spikes1] = binning({X1}, bin_opts);
[~, edges_spikes2] = binning({X2}, bin_opts);
edges_spikes1(1) = 0; edges_spikes2(1) = 0;
edges_spikes1(end+1) = Inf; edges_spikes2(end+1) = Inf;

PID_v_venka = zeros(length(trial_categories), 8, numRep);

parfor rep = 1:numRep
    PID_tmp = zeros(length(trial_categories), 8);
    for trials = 1:length(trial_categories)
        num_trials = trial_categories(trials);
        X1 = [poissrnd(rate1(1),1,num_trials), poissrnd(rate1(2),1,num_trials), ...
              poissrnd(rate1(3),1,num_trials), poissrnd(rate1(4),1,num_trials)];
        X2 = [poissrnd(rate2(1),1,num_trials), poissrnd(rate2(2),1,num_trials), ...
              poissrnd(rate2(3),1,num_trials), poissrnd(rate2(4),1,num_trials)];
        Xshared = [poissrnd(rateShared(1),1,num_trials), poissrnd(rateShared(2),1,num_trials), ...
                   poissrnd(rateShared(3),1,num_trials), poissrnd(rateShared(4),1,num_trials)];
        X1 = X1 + Xshared;
        X2 = X2 + Xshared;

        Y = [ones(1,num_trials), 2*ones(1,num_trials), 3*ones(1,num_trials), 4*ones(1,num_trials)];
        Yrnd = hShuffle({Y}, {'A'});

        X1_binned = discretize(X1, edges_spikes1);
        X2_binned = discretize(X2, edges_spikes2);
        input = {X1_binned, X2_binned, Y};

        outputs = MI({[X1_binned; X2_binned], Y}, {'I(A;B)'}, opts_MI);
        Ijoint = outputs{1};

        outputs = MI({[X1_binned; X2_binned], Y}, {'I(A;B)'}, opts_MI_qe);
        Ijoint_qe = max([0, outputs{1}]);
        debias_factor = Ijoint_qe / Ijoint;

        i1 = MI({X1_binned, Y}, {'I(A;B)'}, opts_MI_qe);
        i2 = MI({X2_binned, Y}, {'I(A;B)'}, opts_MI_qe);

        imx_qe = max([0, i1{1}]);
        imy_qe = max([0, i2{1}]);

        Ijoint_qe = max([Ijoint_qe, imx_qe, imy_qe]);

        PID_res_plugin = PID(input, reqOutputs_all, opts_PID_plugin);

        union_info = (PID_res_plugin{3} + PID_res_plugin{4} + PID_res_plugin{5}) * debias_factor;
        % union_info = max([union_info, imx_qe, imy_qe]);
        % union_info = min([union_info, imx_qe + imy_qe, Ijoint_qe]);

        PID_res_plugin{1} = Ijoint_qe;
        PID_res_plugin{2} = Ijoint_qe - union_info;
        PID_res_plugin{3} = imx_qe + imy_qe - union_info;
        PID_res_plugin{4} = union_info - imy_qe;
        PID_res_plugin{5} = union_info - imx_qe;
        PID_res_plugin{6} = union_info;

        PID_tmp(trials, :) = [PID_res_plugin{1}, PID_res_plugin{2}, PID_res_plugin{3}, ...
                              PID_res_plugin{4}, PID_res_plugin{5}, PID_res_plugin{6}, 0, 0];
    end
    PID_v_venka(:, :, rep) = PID_tmp;
end

folder_name = 'Broja_venka';
if ~isfolder(['Results/', folder_name])
    mkdir(['Results/', folder_name]);
end
save(['Results/', folder_name, '/', 'Simuldata_', char(simul_case), '_', char(info_amount), '_', ...
      char(opts_PID.redundancy_measure), '.mat'], 'PID_v_venka');
end
