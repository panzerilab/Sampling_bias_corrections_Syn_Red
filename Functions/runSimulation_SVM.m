function runSimulation_SVM(nRbins, trial_categories, numRep, simul_case, alpha,info_amount, popSize, redundancyMeasure)
n_neurons = 2;
% Compute Rates
switch simul_case
    case 'uncorr_unique'
        B = 10; beta = 0;gamma_flat = 0;
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
opts_PID.pid_constrained = true;
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
opts_PID_infoCorr.bias = 'infoCorr';


if strcmp(redundancyMeasure, 'I_BROJA')
    reqOutputs_all = {'Joint', 'PID_atoms', 'Union', 'q_dist', 'p_ind'};
else
    reqOutputs_all = {'Joint', 'PID_atoms', 'Union'};
end
reqOutputs = {'Joint', 'PID_atoms', 'Union'};


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
PID_v_qe        = zeros(length(trial_categories),8, numRep);
PID_v_qeShuff   = zeros(length(trial_categories),8, numRep);
PID_v_infoCorr  = zeros(length(trial_categories),8, numRep);


SVM_opts.libsvm = false;
SVM_opts.svm_family = 'linear';
reqOutputs_SVM = {'labels', 'testIdx'};
parfor rep = 1:numRep
    PID_v_tmp_plugin    = zeros(length(trial_categories),8);
    PID_v_tmp_shuff     = zeros(length(trial_categories),8);
    PID_v_tmp_shuffled  = zeros(length(trial_categories),8);
    PID_v_tmp_qe        = zeros(length(trial_categories),8);
    PID_v_tmp_qeShuff   = zeros(length(trial_categories),8);
    PID_v_tmp_infoCorr  = zeros(length(trial_categories),8);
    for trials=1:length(trial_categories)
        num_trials = trial_categories(trials);
        X1 = zeros(popSize, num_trials*4);
        X2 = zeros(popSize, num_trials*4);
        for neuronIdx = 1:popSize
            X1_tmp=[poissrnd(rate1(1),1,num_trials) poissrnd(rate1(2),1,num_trials)  poissrnd(rate1(3),1,num_trials) poissrnd(rate1(4),1,num_trials)];
            X2_tmp=[poissrnd(rate2(1),1,num_trials) poissrnd(rate2(2),1,num_trials)  poissrnd(rate2(3),1,num_trials) poissrnd(rate2(4),1,num_trials)];
            Xshared=[poissrnd(rateShared(1),1,num_trials) poissrnd(rateShared(2),1,num_trials)  poissrnd(rateShared(3),1,num_trials) poissrnd(rateShared(4),1,num_trials)];
            X1 = X1+Xshared;
            X2 = X2+Xshared;
            X1(neuronIdx,:) = X1_tmp;
            X2(neuronIdx,:) = X2_tmp;
        end
        Y=[ones(1,num_trials) 2*ones(1,num_trials) 3*ones(1,num_trials) 4*ones(1,num_trials)];
        Yrnd = hShuffle({Y}, {'A'});
        Yrnd = Yrnd{1};

        testIdx = evenStimulusKFold(Y);

        hola = SVM_opts;
        hola.cv = {'selectedTestIdx', testIdx};

        SVM_out1 = svm_wrapper({X1, Y}, reqOutputs_SVM, hola);
        SVM_out2 = svm_wrapper({X2, Y}, reqOutputs_SVM, hola);

        Y_pred1 = SVM_out1{1};
        Y_pred2 = SVM_out2{1};

        PID_v_plugin_fold    = zeros(size(testIdx,2),8);
        PID_v_shuffled_fold     = zeros(size(testIdx,2),8);
        PID_v_shuff_fold     = zeros(size(testIdx,2),8);
        PID_v_qe_fold        = zeros(size(testIdx,2),8);
        PID_v_qeShuff_fold   = zeros(size(testIdx,2),8);
        PID_v_infoCorr_fold  = zeros(size(testIdx,2),8);



        input = {Y_pred1, Y_pred2, Y};

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
        PID_res_shuffled = PID({Y_pred1, Y_pred2, Yrnd},reqOutputs_all, opts_PID_plugin);
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
    PID_v_plugin(:, :, rep)     = PID_v_tmp_plugin;
    PID_v_shuff(:, :, rep)      = PID_v_tmp_shuff;
    PID_v_shuffled(:, :, rep)   = PID_v_tmp_shuffled;
    PID_v_qe(:, :, rep)         = PID_v_tmp_qe;
    PID_v_qeShuff (:, :, rep)   = PID_v_tmp_qeShuff;
    PID_v_infoCorr(:, :, rep)   = PID_v_tmp_infoCorr;
end


% Step 3: Save results
folder_name = 'SVM_results';

if ~isfolder(['Results/', folder_name])
    mkdir(['Results/', folder_name]);
end

save(['Results/', folder_name, '/', 'Simuldata_', char(simul_case), '_', char(info_amount), '_', char(opts_PID.redundancy_measure), '_', num2str(popSize), '.mat']);

end

function testIdx = evenStimulusKFold(stimuli)

nTrials = length(stimuli);
testIdx = false(nTrials, 2);

uniqueStimuli = unique(stimuli);

for i = 1:length(uniqueStimuli)
    stim = uniqueStimuli(i);
    stimIdx = find(stimuli == stim);
    stimIdx = stimIdx(randperm(length(stimIdx)));
    nHalf = ceil(length(stimIdx) / 2);
    fold1 = stimIdx(1:nHalf);
    fold2 = stimIdx(nHalf+1:end);
    testIdx(fold1, 1) = true;
    testIdx(fold2, 2) = true;
end
end

