function runSimulation(nRbins, trial_categories, numRep, simul_case, alpha,info_amount, isalphasweep, isBinSweep, redundancyMeasure)
n_neurons = 2;
% Compute Rates
switch simul_case
    case 'uncorr_unique'
        B = 5; beta = 0;gamma_flat = 0;
    case 'unique_high_red'
        B = 5; beta = 0.4;gamma_flat = 2;
    case 'unique_high_syn'
        B = 5; beta = 0.1;gamma_flat = 20;
    case 'bit_of_all'
        B = 5; beta = 0.25;gamma_flat = 17.5;
end


opts_PID.bin_method = {'none'};
opts_PID.suppressWarnings = true;
opts_PID.computeNulldist = false;
opts_PID.pid_constrained = false;
opts_PID.shuff = 20;
opts_PID.xtrp = 20;
opts_PID.redundancy_measure = redundancyMeasure;

opts_PID_plugin        = opts_PID;
opts_PID_plugin.bias   = 'plugin';
opts_PID_shuff         = opts_PID;
opts_PID_shuff.bias    = 'shuffSub';
opts_PID_qe            = opts_PID;
opts_PID_qe.bias       = 'qe';
opts_PID_infoCorr      = opts_PID;
opts_PID_infoCorr.bias = 'resample';

if strcmp(redundancyMeasure, 'I_BROJA')
    reqOutputs_all = {'Joint', 'PID_atoms', 'Union', 'q_dist', 'p_ind'};
else
    reqOutputs_all = {'Joint', 'PID_atoms', 'Union'};
end
reqOutputs = {'Joint', 'PID_atoms', 'Union'};

if isalphasweep
    PID_v_plugin    = zeros(length(trial_categories),8, numRep, length(alpha));
    PID_v_shuff     = zeros(length(trial_categories),8, numRep, length(alpha));
    PID_v_shuffled  = zeros(length(trial_categories),8, numRep, length(alpha));
    PID_v_qe        = zeros(length(trial_categories),8, numRep, length(alpha));
    PID_v_qeShuff   = zeros(length(trial_categories),8, numRep, length(alpha));
    PID_v_infoCorr  = zeros(length(trial_categories),8, numRep, length(alpha));
    for alpha_idx = alpha
        eps = 1;
        rate1      = (B + alpha_idx*beta*[1-eps,1-eps,eps,eps] + alpha_idx*(1-beta)*[1-eps,eps,1-eps,eps]);
        rate2      = (B + alpha_idx*beta*[1-eps,eps,1-eps,eps] + alpha_idx*(1-beta)*[1-eps,1-eps,eps,eps]);
        rateShared = (gamma_flat*[1,1,1,1]);
        num_stimuli = length(rate1);

        % Compute Bin Edges
        num_trials = 2048;
        X1=[poissrnd(rate1(1),1,num_trials) poissrnd(rate1(2),1,num_trials)  poissrnd(rate1(3),1,num_trials) poissrnd(rate1(4),1,num_trials)];
        X2=[poissrnd(rate2(1),1,num_trials) poissrnd(rate2(2),1,num_trials)  poissrnd(rate2(3),1,num_trials) poissrnd(rate2(4),1,num_trials)];
        Xshared=[poissrnd(rateShared(1),1,num_trials) poissrnd(rateShared(2),1,num_trials)  poissrnd(rateShared(3),1,num_trials) poissrnd(rateShared(4),1,num_trials)];
        X1 = X1+Xshared;
        X2 = X2+Xshared;
        bin_opts.bin_method = {'eqpop'};
        bin_opts.n_bins = {nRbins};
        [~, edges_spikes1] = binning({X1}, bin_opts);
        [~, edges_spikes2] = binning({X2}, bin_opts);
        edges_spikes1(1) = 0;
        edges_spikes2(1) = 0;
        x = nRbins+1;
        edges_spikes1(x) = Inf;
        edges_spikes2(x) = Inf;

        PID_v_plugin_a    = zeros(length(trial_categories),8, numRep);
        PID_v_shuff_a     = zeros(length(trial_categories),8, numRep);
        PID_v_shuffled_a  = zeros(length(trial_categories),8, numRep);
        PID_v_qe_a        = zeros(length(trial_categories),8, numRep);
        PID_v_qeShuff_a   = zeros(length(trial_categories),8, numRep);
        PID_v_infoCorr_a  = zeros(length(trial_categories),8, numRep);
        for rep = 1:numRep
            PID_v_tmp_plugin    = zeros(length(trial_categories),8);
            PID_v_tmp_shuff     = zeros(length(trial_categories),8);
            PID_v_tmp_shuffled  = zeros(length(trial_categories),8);
            PID_v_tmp_qe        = zeros(length(trial_categories),8);
            PID_v_tmp_qeShuff   = zeros(length(trial_categories),8);
            PID_v_tmp_infoCorr  = zeros(length(trial_categories),8);
            for trials=1:length(trial_categories)
                num_trials = trial_categories(trials);
                X1=[poissrnd(rate1(1),1,num_trials) poissrnd(rate1(2),1,num_trials)  poissrnd(rate1(3),1,num_trials) poissrnd(rate1(4),1,num_trials)];
                X2=[poissrnd(rate2(1),1,num_trials) poissrnd(rate2(2),1,num_trials)  poissrnd(rate2(3),1,num_trials) poissrnd(rate2(4),1,num_trials)];
                Xshared=[poissrnd(rateShared(1),1,num_trials) poissrnd(rateShared(2),1,num_trials)  poissrnd(rateShared(3),1,num_trials) poissrnd(rateShared(4),1,num_trials)];
                X1 = X1+Xshared;
                X2 = X2+Xshared;
                Y=[ones(1,num_trials) 2*ones(1,num_trials) 3*ones(1,num_trials) 4*ones(1,num_trials)];
                Yrnd = hShuffle({Y}, {'A'});
                X1_binned = (discretize(X1, edges_spikes1));
                X2_binned = (discretize(X2, edges_spikes2));
                input = {X1_binned, X2_binned, Y};


                % plugin
                PID_res_plugin = PID(input,reqOutputs_all, opts_PID_plugin);
                if strcmp(redundancyMeasure, 'I_BROJA')
                    Pind = mutualInformationXYZ(PID_res_plugin{8});
                    qDist = mutualInformationXYZ(PID_res_plugin{7});
                else
                    Pind = 0;
                    qDist = 0;
                end
                % shuffled
                PID_res_shuffled = PID({X1_binned, X2_binned, Yrnd{1}},reqOutputs_all, opts_PID_plugin);
                if strcmp(redundancyMeasure, 'I_BROJA')
                    Pind_sh = mutualInformationXYZ(PID_res_shuffled{8});
                    qDist_sh = mutualInformationXYZ(PID_res_shuffled{7});
                else
                    Pind_sh = 0;
                    qDist_sh = 0;
                end
                % shuff Sub and plugin
                PID_res_shuff = PID(input,reqOutputs, opts_PID_shuff);
                % qe
                PID_res_qe = PID(input,reqOutputs, opts_PID_qe);
                % infoCorr
                PID_res_infoCorr = PID(input,reqOutputs, opts_PID_infoCorr);

                PID_v_tmp_plugin(trials, :)   = [PID_res_plugin{1}, PID_res_plugin{2}, PID_res_plugin{3}, PID_res_plugin{4}, PID_res_plugin{5}, PID_res_plugin{6}, Pind, qDist];
                PID_v_tmp_shuffled(trials, :) = [PID_res_shuffled{1}, PID_res_shuffled{2}, PID_res_shuffled{3}, PID_res_shuffled{4}, PID_res_shuffled{5}, PID_res_shuffled{6}, Pind_sh, qDist_sh];
                PID_v_tmp_infoCorr(trials, :) = [PID_res_infoCorr{1}, PID_res_infoCorr{2}, PID_res_infoCorr{3}, PID_res_infoCorr{4}, PID_res_infoCorr{5}, PID_res_infoCorr{6}, 0, 0];
                PID_v_tmp_shuff(trials, :)    = [PID_res_shuff{1}, PID_res_shuff{2}, PID_res_shuff{3}, PID_res_shuff{4}, PID_res_shuff{5}, PID_res_shuff{6}, 0, 0];
                PID_v_tmp_qe(trials, :)       = [PID_res_qe{1}, PID_res_qe{2}, PID_res_qe{3}, PID_res_qe{4}, PID_res_qe{5}, PID_res_qe{6}, 0, 0];
                PID_v_tmp_qeShuff(trials,:)   = mean([PID_res_shuff{1}, PID_res_shuff{2}, PID_res_shuff{3}, PID_res_shuff{4}, PID_res_shuff{5}, PID_res_shuff{6}, 0, 0; ...
                    PID_res_qe{1},    PID_res_qe{2},    PID_res_qe{3},    PID_res_qe{4},    PID_res_qe{5},    PID_res_qe{6},    0, 0], 1);
            end
            PID_v_plugin_a(:, :, rep)     = PID_v_tmp_plugin;
            PID_v_shuff_a(:, :, rep)      = PID_v_tmp_shuff;
            PID_v_shuffled_a(:, :, rep)   = PID_v_tmp_shuffled;
            PID_v_qe_a(:, :, rep)         = PID_v_tmp_qe;
            PID_v_qeShuff_a(:, :, rep)    = PID_v_tmp_qeShuff;
            PID_v_infoCorr_a(:, :, rep)   = PID_v_tmp_infoCorr;
        end
        PID_v_plugin(:,:,:,alpha_idx) = PID_v_plugin_a;
        PID_v_shuff(:,:,:,alpha_idx) = PID_v_shuff_a;
        PID_v_shuffled(:,:,:,alpha_idx) = PID_v_shuffled_a;
        PID_v_qe(:,:,:,alpha_idx) = PID_v_qe_a;
        PID_v_qeShuff(:,:,:,alpha_idx) = PID_v_qeShuff_a;
        PID_v_infoCorr(:,:,:,alpha_idx) = PID_v_infoCorr_a;
    end
elseif isBinSweep

    PID_v_plugin    = zeros(length(trial_categories),8, numRep, length(nRbins));
    PID_v_shuff     = zeros(length(trial_categories),8, numRep, length(nRbins));
    PID_v_qe        = zeros(length(trial_categories),8, numRep, length(nRbins));
    PID_v_qeShuff   = zeros(length(trial_categories),8, numRep, length(nRbins));
    PID_v_infoCorr  = zeros(length(trial_categories),8, numRep, length(nRbins));


    rate1_stim = B + alpha*beta*[0,0,0,1,1,1,2,2,2] + alpha*(1-beta)*[0,1,2,0,1,2,0,1,2];
    rate2_stim = B + alpha*beta*[0,1,2,0,1,2,0,1,2] + alpha*(1-beta)*[0,0,0,1,1,1,2,2,2];
    rateShare  = gamma_flat*[1 1 1 1 1 1 1 1 1];

    for binIdx = 1:length(nRbins)
        binIdx

        IdxRate = nRbins(binIdx);
        rateShared = rateShare(1:IdxRate);
        rate1 = rate1_stim(1:IdxRate);
        rate2 = rate2_stim(1:IdxRate);
        num_stimuli = IdxRate;
        num_trials = 10000;
        X1 = [];
        X2 = [];
        Xshared = [];
        for stimIdx = 1:num_stimuli
            X1=[X1, poissrnd(rate1(stimIdx),1,num_trials)];
            X2=[X2, poissrnd(rate2(stimIdx),1,num_trials)];
            Xshared=[Xshared, poissrnd(rateShared(stimIdx),1,num_trials)];
        end
        X1 = X1+Xshared;
        X2 = X2+Xshared;
        bin_opts.bin_method = {'eqpop'};
        bin_opts.n_bins = {nRbins(binIdx)};
        [~, edges_spikes1] = binning({X1}, bin_opts);
        [~, edges_spikes2] = binning({X2}, bin_opts);
        edges_spikes1(1) = 0;
        edges_spikes2(1) = 0;
        x = nRbins(binIdx)+1;
        edges_spikes1(x) = Inf;
        edges_spikes2(x) = Inf;

        PID_v_plugin_b    = zeros(length(trial_categories),8, numRep);
        PID_v_shuff_b     = zeros(length(trial_categories),8, numRep);
        PID_v_qe_b        = zeros(length(trial_categories),8, numRep);
        PID_v_qeShuff_b   = zeros(length(trial_categories),8, numRep);
        PID_v_infoCorr_b  = zeros(length(trial_categories),8, numRep);

        parfor rep = 1:numRep
            PID_v_tmp_plugin    = zeros(length(trial_categories),8);
            PID_v_tmp_shuff     = zeros(length(trial_categories),8);
            PID_v_tmp_qe        = zeros(length(trial_categories),8);
            PID_v_tmp_qeShuff   = zeros(length(trial_categories),8);
            PID_v_tmp_infoCorr  = zeros(length(trial_categories),8);

            for trials=1:length(trial_categories)
                num_trials_total = trial_categories(trials);

                base_trials = floor(num_trials_total / num_stimuli);
                remaining = num_trials_total - base_trials * num_stimuli;
                trial_counts = base_trials * ones(1, num_stimuli);
                trial_counts(randperm(num_stimuli, remaining)) = trial_counts(randperm(num_stimuli, remaining)) + 1;

                X1 = [];
                X2 = [];
                Xshared = [];
                Y = [];

                for stimIdx = 1:num_stimuli
                    n_trials = trial_counts(stimIdx);
                    X1 = [X1, poissrnd(rate1(stimIdx), 1, n_trials)];
                    X2 = [X2, poissrnd(rate2(stimIdx), 1, n_trials)];
                    Xshared = [Xshared, poissrnd(rateShared(stimIdx), 1, n_trials)];
                    Y = [Y, stimIdx * ones(1, n_trials)];
                end

                X1 = X1 + Xshared;
                X2 = X2 + Xshared;
                Yrnd = hShuffle({Y}, {'A'});
                X1_binned = (discretize(X1, edges_spikes1));
                X2_binned = (discretize(X2, edges_spikes2));
                input = {X1_binned, X2_binned, Y};
                % plugin
                PID_res_plugin = PID(input,reqOutputs_all, opts_PID_plugin);
                if strcmp(redundancyMeasure, 'I_BROJA')
                    Pind = mutualInformationXYZ(PID_res_plugin{8});
                    qDist = mutualInformationXYZ(PID_res_plugin{7});
                else
                    Pind = 0;
                    qDist = 0;
                end
                % shuffled
                PID_res_shuffled = PID({X1_binned, X2_binned, Yrnd{1}},reqOutputs_all, opts_PID_plugin);
                if strcmp(redundancyMeasure, 'I_BROJA')
                    Pind_sh = mutualInformationXYZ(PID_res_plugin{8});
                    qDist_sh = mutualInformationXYZ(PID_res_plugin{7});
                else
                    Pind_sh = 0;
                    qDist_sh = 0;
                end
                % shuff Sub and plugin
                PID_res_shuff = PID(input,reqOutputs, opts_PID_shuff);
                % qe
                PID_res_qe = PID(input,reqOutputs, opts_PID_qe);
                % infoCorr
                PID_res_infoCorr = PID(input,reqOutputs, opts_PID_infoCorr);

                PID_v_tmp_plugin(trials, :)   = [PID_res_plugin{1}, PID_res_plugin{2}, PID_res_plugin{3}, PID_res_plugin{4}, PID_res_plugin{5}, PID_res_plugin{6}, Pind, qDist];
                PID_v_tmp_infoCorr(trials, :) = [PID_res_infoCorr{1}, PID_res_infoCorr{2}, PID_res_infoCorr{3}, PID_res_infoCorr{4}, PID_res_infoCorr{5}, PID_res_infoCorr{6}, 0, 0];
                PID_v_tmp_shuff(trials, :)    = [PID_res_shuff{1}, PID_res_shuff{2}, PID_res_shuff{3}, PID_res_shuff{4}, PID_res_shuff{5}, PID_res_shuff{6}, 0, 0];
                PID_v_tmp_qe(trials, :)       = [PID_res_qe{1}, PID_res_qe{2}, PID_res_qe{3}, PID_res_qe{4}, PID_res_qe{5}, PID_res_qe{6}, 0, 0];
                PID_v_tmp_qeShuff(trials,:)   = mean([PID_res_shuff{1}, PID_res_shuff{2}, PID_res_shuff{3}, PID_res_shuff{4}, PID_res_shuff{5}, PID_res_shuff{6}, 0, 0; ...
                    PID_res_qe{1},    PID_res_qe{2},    PID_res_qe{3},    PID_res_qe{4},    PID_res_qe{5},    PID_res_qe{6},    0, 0], 1);
            end
            PID_v_plugin_b(:, :, rep)     = PID_v_tmp_plugin;
            PID_v_shuff_b(:, :, rep)      = PID_v_tmp_shuff;
            PID_v_qe_b(:, :, rep)         = PID_v_tmp_qe;
            PID_v_qeShuff_b(:, :, rep)   = PID_v_tmp_qeShuff;
            PID_v_infoCorr_b(:, :, rep)   = PID_v_tmp_infoCorr;
        end
        PID_v_plugin(:,:,:,binIdx) = PID_v_plugin_b;
        PID_v_shuff(:,:,:,binIdx) = PID_v_shuff_b;
        PID_v_qe(:,:,:,binIdx) = PID_v_qe_b;
        PID_v_qeShuff(:,:,:,binIdx) = PID_v_qeShuff_b;
        PID_v_infoCorr(:,:,:,binIdx) = PID_v_infoCorr_b;
    end
else
    eps = 1;
    rate1 = B + alpha*beta*[1-eps,1-eps,eps,eps] + alpha*(1-beta)*[1-eps,eps,1-eps,eps];
    rate2 = B + alpha*beta*[1-eps,eps,1-eps,eps] + alpha*(1-beta)*[1-eps,1-eps,eps,eps];
    rateShared = gamma_flat*[1,1,1,1];
    num_stimuli = length(rate1);

    % Compute Bin Edges
    num_trials = 2048;
    X1=[poissrnd(rate1(1),1,num_trials) poissrnd(rate1(2),1,num_trials)  poissrnd(rate1(3),1,num_trials) poissrnd(rate1(4),1,num_trials)];
    X2=[poissrnd(rate2(1),1,num_trials) poissrnd(rate2(2),1,num_trials)  poissrnd(rate2(3),1,num_trials) poissrnd(rate2(4),1,num_trials)];
    Xshared=[poissrnd(rateShared(1),1,num_trials) poissrnd(rateShared(2),1,num_trials)  poissrnd(rateShared(3),1,num_trials) poissrnd(rateShared(4),1,num_trials)];
    X1 = X1+Xshared;
    X2 = X2+Xshared;
    bin_opts.bin_method = {'eqpop'};
    bin_opts.n_bins = {nRbins};
    [~, edges_spikes1] = binning({X1}, bin_opts);
    [~, edges_spikes2] = binning({X2}, bin_opts);
    edges_spikes1(1) = 0;
    edges_spikes2(1) = 0;
    x = nRbins+1;
    edges_spikes1(x) = Inf;
    edges_spikes2(x) = Inf;

    PID_v_plugin    = zeros(length(trial_categories),8, numRep);
    PID_v_shuff     = zeros(length(trial_categories),8, numRep);
    PID_v_shuffled  = zeros(length(trial_categories),8, numRep);
    PID_v_qe        = zeros(length(trial_categories),8, numRep);
    PID_v_qeShuff   = zeros(length(trial_categories),8, numRep);
    PID_v_infoCorr  = zeros(length(trial_categories),8, numRep);
    parfor rep = 1:numRep
        PID_v_tmp_plugin    = zeros(length(trial_categories),8);
        PID_v_tmp_shuff     = zeros(length(trial_categories),8);
        PID_v_tmp_shuffled  = zeros(length(trial_categories),8);
        PID_v_tmp_qe        = zeros(length(trial_categories),8);
        PID_v_tmp_qeShuff   = zeros(length(trial_categories),8);
        PID_v_tmp_infoCorr  = zeros(length(trial_categories),8);
        for trials=1:length(trial_categories)
            num_trials = trial_categories(trials);
            X1=[poissrnd(rate1(1),1,num_trials) poissrnd(rate1(2),1,num_trials)  poissrnd(rate1(3),1,num_trials) poissrnd(rate1(4),1,num_trials)];
            X2=[poissrnd(rate2(1),1,num_trials) poissrnd(rate2(2),1,num_trials)  poissrnd(rate2(3),1,num_trials) poissrnd(rate2(4),1,num_trials)];
            Xshared=[poissrnd(rateShared(1),1,num_trials) poissrnd(rateShared(2),1,num_trials)  poissrnd(rateShared(3),1,num_trials) poissrnd(rateShared(4),1,num_trials)];
            X1 = X1+Xshared;
            X2 = X2+Xshared;
            Y=[ones(1,num_trials) 2*ones(1,num_trials) 3*ones(1,num_trials) 4*ones(1,num_trials)];
            Yrnd = Y(randperm(length(Y))); %hShuffle({Y}, {'A'});
            X1_binned = (discretize(X1, edges_spikes1));
            X2_binned = (discretize(X2, edges_spikes2));
            input = {X1_binned, X2_binned, Y};


            % plugin
            PID_res_plugin = PID(input,reqOutputs_all, opts_PID_plugin);
            % Single1_plugin = PID({X1_binned, Y},{'I(A;B)'}, opts_PID_plugin);
            % Single2_plugin = PID({X2_binned, Y},{'I(A;B)'}, opts_PID_plugin);
            if strcmp(redundancyMeasure, 'I_BROJA')
                Pind = mutualInformationXYZ(PID_res_plugin{8});
                qDist = mutualInformationXYZ(PID_res_plugin{7});
            else
                Pind = 0;
                qDist = 0;
            end
            % shuffled
            PID_res_shuffled = PID({X1_binned, X2_binned, Yrnd},reqOutputs_all, opts_PID_plugin);
            % Single1_shuffled = PID({X1_binned, Yrnd},{'I(A;B)'}, opts_PID_plugin);
            % Single2_shuffled = PID({X2_binned, Yrnd},{'I(A;B)'}, opts_PID_plugin);
            if strcmp(redundancyMeasure, 'I_BROJA')
                Pind_sh = mutualInformationXYZ(PID_res_shuffled{8});
                qDist_sh = mutualInformationXYZ(PID_res_shuffled{7});
            else
                Pind_sh = 0;
                qDist_sh = 0;
            end
            % shuff Sub and plugin
            PID_res_shuff = PID(input,reqOutputs, opts_PID_shuff);
            % qe
            PID_res_qe    = PID(input,reqOutputs, opts_PID_qe);
            % infoCorr
            PID_res_infoCorr = PID(input,reqOutputs, opts_PID_infoCorr);

            PID_v_tmp_plugin(trials, :)   = [PID_res_plugin{1}, PID_res_plugin{2}, PID_res_plugin{3}, PID_res_plugin{4}, PID_res_plugin{5}, PID_res_plugin{6}, Pind, qDist];
            PID_v_tmp_shuffled(trials, :) = [PID_res_shuffled{1}, PID_res_shuffled{2}, PID_res_shuffled{3}, PID_res_shuffled{4}, PID_res_shuffled{5}, PID_res_shuffled{6}, Pind_sh, qDist_sh];
            PID_v_tmp_infoCorr(trials, :) = [PID_res_infoCorr{1}, PID_res_infoCorr{2}, PID_res_infoCorr{3}, PID_res_infoCorr{4}, PID_res_infoCorr{5}, PID_res_infoCorr{6}, 0, 0];
            PID_v_tmp_shuff(trials, :)    = [PID_res_shuff{1}, PID_res_shuff{2}, PID_res_shuff{3}, PID_res_shuff{4}, PID_res_shuff{5}, PID_res_shuff{6}, 0, 0];
            PID_v_tmp_qe(trials, :)       = [PID_res_qe{1}, PID_res_qe{2}, PID_res_qe{3}, PID_res_qe{4}, PID_res_qe{5}, PID_res_qe{6}, 0, 0];
            PID_v_tmp_qeShuff(trials,:)   = mean([PID_res_shuff{1}, PID_res_shuff{2}, PID_res_shuff{3}, PID_res_shuff{4}, PID_res_shuff{5}, PID_res_shuff{6}, 0, 0; ...
                PID_res_qe{1},    PID_res_qe{2},    PID_res_qe{3},    PID_res_qe{4},    PID_res_qe{5},    PID_res_qe{6}, 0, 0], 1);
        end
        PID_v_plugin(:, :, rep)     = PID_v_tmp_plugin;
        PID_v_shuff(:, :, rep)      = PID_v_tmp_shuff;
        PID_v_shuffled(:, :, rep)   = PID_v_tmp_shuffled;
        PID_v_qe(:, :, rep)         = PID_v_tmp_qe;
        PID_v_qeShuff (:, :, rep)   = PID_v_tmp_qeShuff;
        PID_v_infoCorr(:, :, rep)   = PID_v_tmp_infoCorr;
    end
end

% Step 3: Save results
if strcmp(opts_PID.redundancy_measure, 'I_BROJA')
    if isalphasweep
        folder_name = 'Alpha_sweep';
    elseif isBinSweep
        folder_name = 'Bin_sweep';
    else
        folder_name = 'Broja';
    end
elseif strcmp(opts_PID.redundancy_measure, 'I_min')
    folder_name = 'Imin';
elseif strcmp(opts_PID.redundancy_measure, 'I_MMI')
    folder_name = 'IMMI';
end

if ~isfolder(['Results/', folder_name])
    mkdir(['Results/', folder_name]);
end

if ~alpha_sweep & strcmp(info_amount, 'alpha_sweep')
    save(['Results/', folder_name, '/', 'Simuldata_', char(simul_case), '_', char(info_amount), '_', char(opts_PID.redundancy_measure), '_', char(alpha_values), '.mat']);
else
    save(['Results/', folder_name, '/', 'Simuldata_', char(simul_case), '_', char(info_amount), '_', char(opts_PID.redundancy_measure), '.mat']);
end

end