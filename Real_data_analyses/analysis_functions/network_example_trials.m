%rng('default')  % For reproducibility
% GENERATE STIMULUS
nTimesteps = 100;
nNeurons = 6;

PID_opts.bias = 'qe';
PID_opts.xtrp=1;
PID_opts.bin_methodX1 = 'eqpop';
PID_opts.bin_methodX2 = 'eqpop';
PID_opts.bin_methodY  = 'eqpop';
nbins = 4;
PID_opts.n_binsX1 = nbins;             
PID_opts.n_binsX2 = nbins;             
PID_opts.n_binsY  = nbins;       

rho23 = .9;
rho45 = .1;

%% Target neuron 1 and 6. 5 trials
disp('5 trials')
nTrials = 5;

target_neuron1 = 1;
pair_indexes1 = nchoosek(setdiff(1:nNeurons,target_neuron1),2);
PID_naivelowtrialn1 = zeros(nTrials,length(pair_indexes1),4);
PID_qelowtrialn1 = zeros(nTrials,length(pair_indexes1),4);
target_neuron6 = 6;
pair_indexes6 = nchoosek(setdiff(1:nNeurons,target_neuron6),2);
PID_naivelowtrialn6 = zeros(nTrials,length(pair_indexes1),4);
PID_qelowtrialn6 = zeros(nTrials,length(pair_indexes1),4);


for trial=1:nTrials   
    fprintf("Trial %d\n", trial)
    neural_activity = randn(nNeurons,nTimesteps);
    % noise23 = mvnrnd([0,0], [1 rho23; rho23 1], nTimesteps)';
    noise45 = mvnrnd([0,0], [1 rho45; rho45 1], nTimesteps)';
    % neural_activity(2:3, :) = noise23;
    neural_activity(4:5, :) = noise45;
    
    % neural_activity(1, 2:end) = neural_activity(1, 2:end)+ sum(neural_activity(2:3, 1:end-1),1);
    neural_activity(6, 2:end) = neural_activity(6, 2:end)+ sum(neural_activity(4:5, 1:end-1),1);

    PID_naive_perTrial = zeros(length(pair_indexes1),4);
    PID_qe_perTrial = zeros(length(pair_indexes1),4);

    parfor pair_i = 1:length(pair_indexes1)
        n1 = pair_indexes1(pair_i,1);
        n2 = pair_indexes1(pair_i,2);
        [pqe, ~, pnaive] = PID(neural_activity(target_neuron1,2:end),neural_activity(n1,1:end-1),neural_activity(n2,1:end-1),PID_opts);
        PID_naive_perTrial(pair_i,:) = pnaive;
        PID_qe_perTrial(pair_i,:) = pqe;
    end

    PID_naivelowtrialn1(trial,:,:) =  PID_naive_perTrial;
    PID_qelowtrialn1(trial,:,:) = PID_qe_perTrial;

    PID_naive_perTrial = zeros(length(pair_indexes1),4);
    PID_qe_perTrial = zeros(length(pair_indexes1),4);

    parfor pair_i = 1:length(pair_indexes6)
        n1 = pair_indexes6(pair_i,1);
        n2 = pair_indexes6(pair_i,2);
        [pqe, ~, pnaive] = PID(neural_activity(target_neuron6,2:end),neural_activity(n1,1:end-1),neural_activity(n2,1:end-1),PID_opts);
        PID_naive_perTrial(pair_i,:) = pnaive;
        PID_qe_perTrial(pair_i,:) = pqe;
    end
    PID_naivelowtrialn6(trial,:,:) =  PID_naive_perTrial;
    PID_qelowtrialn6(trial,:,:) = PID_qe_perTrial;

end

mean_syn1low = squeeze(mean(PID_naivelowtrialn1(:,:,4), 1));
sem_syn1low = squeeze(std(PID_naivelowtrialn1(:,:,4),[], 2));
mean_synqe1low = squeeze(mean(PID_qelowtrialn1(:,:,4), 1));
sem_synqe1low = squeeze(std(PID_qelowtrialn1(:,:,4), [], 2));
syn_tabletarget1lowtrialn = table(pair_indexes1(:,1), pair_indexes1(:,2), mean_syn1low',mean_synqe1low', 'VariableNames',["Neuron 1","Neuron 2", "Syn plugin", "Syn QE"]);


mean_syn6low = mean(PID_naivelowtrialn6(:,:,4), 1);
sem_syn6low = std(PID_naivelowtrialn6(:,:,4),[], 2);
mean_synqe6low = squeeze(mean(PID_qelowtrialn6(:,:,4), 1));
sem_synqe6low = std(PID_qelowtrialn6(:,:,4),[], 2);
syn_tabletarget6lowtrialn = table(pair_indexes6(:,1), pair_indexes6(:,2), mean_syn6low',mean_synqe6low', 'VariableNames',["Neuron 1","Neuron 2", "Syn plugin", "Syn QE"]);
%% Target neuron 1 and 6. 50 trials
disp('50 trials')
nTrials = 50;

target_neuron1 = 1;
pair_indexes1 = nchoosek(setdiff(1:nNeurons,target_neuron1),2);
PID_naivehightrialn1 = zeros(nTrials,length(pair_indexes1),4);
PID_qehightrialn1 = zeros(nTrials,length(pair_indexes1),4);
target_neuron6 = 6;
pair_indexes6 = nchoosek(setdiff(1:nNeurons,target_neuron6),2);
PID_naivehightrialn6 = zeros(nTrials,length(pair_indexes6),4);
PID_qehightrialn6 = zeros(nTrials,length(pair_indexes6),4);


for trial=1:nTrials   
    fprintf("Trial %d\n", trial)
    neural_activity = randn(nNeurons,nTimesteps);
    % noise23 = mvnrnd([0,0], [1 rho23; rho23 1], nTimesteps)';
    noise45 = mvnrnd([0,0], [1 rho45; rho45 1], nTimesteps)';
    % neural_activity(2:3, :) = noise23;
    neural_activity(4:5, :) = noise45;
    
    % neural_activity(1, 2:end) = neural_activity(1, 2:end)+ sum(neural_activity(2:3, 1:end-1),1);
    neural_activity(6, 2:end) = neural_activity(6, 2:end)+ sum(neural_activity(4:5, 1:end-1),1);

    PID_naive_perTrial = zeros(length(pair_indexes1),4);
    PID_qe_perTrial = zeros(length(pair_indexes1),4);

    parfor pair_i = 1:length(pair_indexes1)
        n1 = pair_indexes1(pair_i,1);
        n2 = pair_indexes1(pair_i,2);
        [pqe, ~, pnaive] = PID(neural_activity(target_neuron1,2:end),neural_activity(n1,1:end-1),neural_activity(n2,1:end-1),PID_opts);
        PID_naive_perTrial(pair_i,:) = pnaive;
        PID_qe_perTrial(pair_i,:) = pqe;
    end

    PID_naivehightrialn1(trial,:,:) =  PID_naive_perTrial;
    PID_qehightrialn1(trial,:,:) = PID_qe_perTrial;

    PID_naive_perTrial = zeros(length(pair_indexes1),4);
    PID_qe_perTrial = zeros(length(pair_indexes1),4);

    parfor pair_i = 1:length(pair_indexes6)
        n1 = pair_indexes6(pair_i,1);
        n2 = pair_indexes6(pair_i,2);
        [pqe, ~, pnaive] = PID(neural_activity(target_neuron6,2:end),neural_activity(n1,1:end-1),neural_activity(n2,1:end-1),PID_opts);
        PID_naive_perTrial(pair_i,:) = pnaive;
        PID_qe_perTrial(pair_i,:) = pqe;
    end
    PID_naivehightrialn6(trial,:,:) =  PID_naive_perTrial;
    PID_qehightrialn6(trial,:,:) = PID_qe_perTrial;

end

mean_syn1high = squeeze(mean(PID_naivehightrialn1(:,:,4), 1));
sem_syn1high = squeeze(std(PID_naivehightrialn1(:,:,4),[], 2));
mean_synqe1high = squeeze(mean(PID_qehightrialn1(:,:,4), 1));
sem_synqe1high = squeeze(std(PID_qehightrialn1(:,:,4), [], 2));
syn_tabletarget1hightrialn = table(pair_indexes1(:,1), pair_indexes1(:,2), mean_syn1high',mean_synqe1high', 'VariableNames',["Neuron 1","Neuron 2", "Syn plugin", "Syn QE"]);


mean_syn6high = mean(PID_naivehightrialn6(:,:,4), 1);
sem_syn6high = std(PID_naivehightrialn6(:,:,4),[], 2);
mean_synqe6high = squeeze(mean(PID_qehightrialn6(:,:,4), 1));
sem_synqe6high = std(PID_qehightrialn6(:,:,4),[], 2);
syn_tabletarget6hightrialn = table(pair_indexes6(:,1), pair_indexes6(:,2), mean_syn6high',mean_synqe6high', 'VariableNames',["Neuron 1","Neuron 2", "Syn plugin", "Syn QE"]);

%%
figure;
max_val = .3; %max([max(mean_syn1low), max(mean_syn6low), max(mean_syn1high), max(mean_syn6high)]);
subplot(2,4,1)
title('Synergy')
heatmap(syn_tabletarget1hightrialn, "Neuron 1","Neuron 2", "ColorVariable", "Syn plugin", 'MissingDataColor', 'white')
clim([0, max_val])

subplot(2,4,2)
title('Synergy corrected')
heatmap(syn_tabletarget1lowtrialn, "Neuron 1","Neuron 2", "ColorVariable", "Syn QE", 'MissingDataColor', 'white')
clim([0, max_val])

subplot(2,4,3)
title('Synergy')
heatmap(syn_tabletarget6lowtrialn, "Neuron 1","Neuron 2", "ColorVariable", "Syn plugin", 'MissingDataColor', 'white')
clim([0, max_val])

subplot(2,4,4)
title('Synergy corrected')
heatmap(syn_tabletarget6lowtrialn, "Neuron 1","Neuron 2", "ColorVariable", "Syn QE", 'MissingDataColor', 'white')
clim([0, max_val])

subplot(2,4,5)
title('Synergy')
heatmap(syn_tabletarget1hightrialn, "Neuron 1","Neuron 2", "ColorVariable", "Syn plugin", 'MissingDataColor', 'white')
clim([0, max_val])

subplot(2,4,6)
title('Synergy corrected')
heatmap(syn_tabletarget1hightrialn, "Neuron 1","Neuron 2", "ColorVariable", "Syn QE", 'MissingDataColor', 'white')
clim([0, max_val])

subplot(2,4,7)
title('Synergy')
heatmap(syn_tabletarget6hightrialn, "Neuron 1","Neuron 2", "ColorVariable", "Syn plugin", 'MissingDataColor', 'white')
clim([0, max_val])

subplot(2,4,8)
title('Synergy corrected')
heatmap(syn_tabletarget6hightrialn, "Neuron 1","Neuron 2", "ColorVariable", "Syn QE", 'MissingDataColor', 'white')
clim([0, max_val])

%% Significance tests

ptestnaive1low  = zeros(3,length(pair_indexes1));
ptestnaive6low  = zeros(3,length(pair_indexes6));

ptestnaive1high = zeros(3,length(pair_indexes1));
ptestnaive6high = zeros(3,length(pair_indexes6));

ptestqe1low     = zeros(3,length(pair_indexes1));
ptestqe6low     = zeros(3,length(pair_indexes6));

ptestqe1high    = zeros(3,length(pair_indexes1));
ptestqe6high    = zeros(3,length(pair_indexes6));

for i=1:length(pair_indexes1)
    ptestnaive1low(1:2,i)  = pair_indexes1(i,:);
    ptestnaive6low(1:2,i)  = pair_indexes6(i,:);
    ptestnaive1high(1:2,i)  = pair_indexes1(i,:);
    ptestnaive6high(1:2,i)  = pair_indexes6(i,:);
    ptestqe1low(1:2,i)  = pair_indexes1(i,:);
    ptestqe6low(1:2,i)  = pair_indexes6(i,:);
    ptestqe1high(1:2,i)  = pair_indexes1(i,:);
    ptestqe6high(1:2,i)  = pair_indexes6(i,:);
    [~, ptestnaive1low(3,:) ]= ttest(squeeze(PID_naivelowtrialn1(:,:,4)),0,"Tail","right");
    [~, ptestnaive6low(3,:) ]= ttest(squeeze(PID_naivelowtrialn6(:,:,4)),0,"Tail","right");
    [~, ptestnaive1high(3,:)] = ttest(squeeze(PID_naivehightrialn1(:,:,4)),0,"Tail","right");
    [~, ptestnaive6high(3,:)] = ttest(squeeze(PID_naivehightrialn6(:,:,4)),0,"Tail","right");
    [~, ptestqe1low(3,:) ]= ttest(squeeze(PID_qelowtrialn1(:,:,4)),0,"Tail","right");
    [~, ptestqe6low(3,:) ]= ttest(squeeze(PID_qelowtrialn6(:,:,4)),0,"Tail","right");
    [~, ptestqe1high(3,:)] = ttest(squeeze(PID_qehightrialn1(:,:,4)),0,"Tail","right");
    [~, ptestqe6high(3,:)] = ttest(squeeze(PID_qehightrialn6(:,:,4)),0,"Tail","right");
end

