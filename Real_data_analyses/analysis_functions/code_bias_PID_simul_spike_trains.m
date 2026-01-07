function code_bias_PID_simul_spike_trains(nRbins, trial_categories, n_iter, simul_case, bias_correction_ALL,info_amount,redundancy_measure, alpha, isalphasweep)
n_neurons = 2;
sharedRate_Stim = [];
rateStim = cell(1,n_neurons);
switch simul_case
    case 'uncorr_unique'
        dir_to_save = ['Figures/',simul_case,'/']; mkdir(dir_to_save);
        B = 5; beta = 0;
        gamma_flat = 0;
        [rateStim{1}, rateStim{2}, sharedRate_Stim] = params_simulations(B, alpha, beta, gamma_flat);

    case 'unique_high_red'
        dir_to_save = ['Figures/',simul_case,'/']; mkdir(dir_to_save);
        B = 5; beta = 0.4;
        gamma_flat = 2;
        [rateStim{1}, rateStim{2}, sharedRate_Stim] = params_simulations(B, alpha, beta, gamma_flat);

    case 'unique_high_syn'
        dir_to_save = ['Figures/',simul_case,'/']; mkdir(dir_to_save);
        B = 5; beta = 0.1;
        gamma_flat = 20;
        [rateStim{1}, rateStim{2}, sharedRate_Stim] = params_simulations(B, alpha, beta, gamma_flat);
end
%% Simulate the time series and compute the Infobreakdown and PID
n_stimuli = length(rateStim{1});
meanRates_trials_iter = cell(1, n_iter);
samples_X1_iter = cell(n_iter,length(trial_categories));
samples_X2_iter = cell(n_iter,length(trial_categories));
samples_states_Y_iter = cell(n_iter,length(trial_categories));
I_iter = zeros(n_iter,length(trial_categories));
ILIN_iter = zeros(n_iter,length(trial_categories));
ISS_iter = zeros(n_iter,length(trial_categories));
ICI_iter = zeros(n_iter,length(trial_categories));
ICD_iter = zeros(n_iter,length(trial_categories));
Ish_iter = zeros(n_iter,length(trial_categories));
Joint_SI_PID_iter = zeros(n_iter,length(trial_categories));
PID_shared_iter = zeros(n_iter,length(trial_categories));
PID_complementary_iter = zeros(n_iter,length(trial_categories));
PID_uniqueX1_iter = zeros(n_iter,length(trial_categories));
PID_uniqueX2_iter = zeros(n_iter,length(trial_categories));
Joint_prob_ind_iter = zeros(n_iter,length(trial_categories));
Joint_minQ_iter = zeros(n_iter,length(trial_categories));
[edges_spikes1uc, edges_spikes2uc] = simulate_edges(nRbins, n_neurons, n_stimuli, rateStim, 0,  2048, false);
[edges_spikes1c, edges_spikes2c] = simulate_edges(nRbins, n_neurons, n_stimuli, rateStim,sharedRate_Stim, 2048, true);

parfor iter = 1:n_iter
    meanRates_trials = cell(1, length(trial_categories));
    samples_X1 = cell(1, length(trial_categories));
    samples_X2 = cell(1, length(trial_categories));
    samples_states_Y = cell(1, length(trial_categories));
    I = zeros(1, length(trial_categories));
    ILIN = zeros(1, length(trial_categories));
    ISS = zeros(1, length(trial_categories));
    ICI = zeros(1, length(trial_categories));
    ICD = zeros(1, length(trial_categories));
    Ish = zeros(1, length(trial_categories));
    Joint_SI_PID = zeros(1, length(trial_categories));
    Joint_prob_ind = zeros(1, length(trial_categories));
    PID_shared = zeros(1, length(trial_categories));
    PID_complementary = zeros(1, length(trial_categories));
    PID_uniqueX1 = zeros(1, length(trial_categories));
    PID_uniqueX2 = zeros(1, length(trial_categories));
    Joint_minQ = zeros(1, length(trial_categories));

    for t = length(trial_categories):-1:1
        % if t== length(trial_categories) && iter >1
        %     continue
        % else
        n_trials = trial_categories(t);
        switch simul_case
            case {'uncorr_unique'}
                [~, ~, meanRates_trials{t}, ~, I(t), ILIN(t), ISS(t), ICI(t), ICD(t), Ish(t), ...
                    Joint_SI_PID(t), PID_shared(t), PID_complementary(t), PID_uniqueX1(t), PID_uniqueX2(t), Joint_prob_ind(t), Joint_minQ(t)] = ...
                    uncorrelated_analyses_bias_PID(info_amount,nRbins, bias_correction_ALL, n_neurons, n_stimuli, rateStim,n_trials,redundancy_measure, edges_spikes1uc, edges_spikes2uc);
            case {'unique_high_red','unique_high_syn'}
                [~, ~, meanRates_trials{t}, ~, I(t), ILIN(t), ISS(t), ICI(t), ICD(t), Ish(t), ...
                    Joint_SI_PID(t), PID_shared(t), PID_complementary(t), PID_uniqueX1(t), PID_uniqueX2(t), Joint_prob_ind(t), Joint_minQ(t)] = ...
                    correlated_analyses_bias_PID(info_amount,nRbins, bias_correction_ALL, n_neurons, n_stimuli, rateStim,sharedRate_Stim,n_trials,redundancy_measure, edges_spikes1c, edges_spikes2c);
        end
    end

    meanRates_trials_iter{iter} = meanRates_trials;
    samples_X1_iter{iter} = samples_X1;
    samples_X2_iter{iter} = samples_X2;
    samples_states_Y_iter{iter} = samples_states_Y;
    I_iter(iter,:) = I;
    ILIN_iter(iter,:) = ILIN;
    ISS_iter(iter,:) = ISS;
    ICI_iter(iter,:) = ICI;
    ICD_iter(iter,:) = ICD;
    Ish_iter(iter,:) = Ish;
    Joint_SI_PID_iter(iter,:) = Joint_SI_PID;
    PID_shared_iter(iter,:) = PID_shared;
    PID_complementary_iter(iter,:) = PID_complementary;
    PID_uniqueX1_iter(iter,:) = PID_uniqueX1;
    PID_uniqueX2_iter(iter,:) = PID_uniqueX2;
    Joint_prob_ind_iter(iter,:) = Joint_prob_ind;
    Joint_minQ_iter(iter,:) = Joint_minQ;
end

if strcmp(redundancy_measure, 'I_BROJA')
    if isalphasweep
        folder_name = 'Alpha_sweep';
    else
        if nRbins == 4
            folder_name = 'Broja';
        else
            folder_name = ['Broja_', num2str(nRbins), 'Bins'];
        end
    end
elseif strcmp(redundancy_measure, 'I_min')
    folder_name = 'Imin';
elseif strcmp(redundancy_measure, 'I_MMI')
    folder_name = 'IMMI';
end

if ~isfolder(['Results/', folder_name])
    mkdir(['Results/', folder_name]);
end

if isalphasweep
    save(['Results/', folder_name, '/', 'Simuldata_', char(simul_case), '_', char(bias_correction_ALL), '_', char(info_amount), '_', char(redundancy_measure),'_', num2str(alpha), '.mat']);
else
    save(['Results/', folder_name, '/', 'Simuldata_', char(simul_case), '_', char(bias_correction_ALL), '_', char(info_amount), '_', char(redundancy_measure), '.mat']);
end

end
