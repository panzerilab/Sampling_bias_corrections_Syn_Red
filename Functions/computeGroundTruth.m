function computeGroundTruth(nRbins, simul_case,alpha,info_amount, isalphasweep, isBinSweep, redundancy_measure)
n_neurons = 2;
dir_to_save = ['Figures/',simul_case,'/'];
if ~exist(dir_to_save, 'dir')
    mkdir(dir_to_save);
end

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

rate1_stim = B + alpha*beta*[0,0,0,1,1,1,2,2,2] + alpha*(1-beta)*[0,1,2,0,1,2,0,1,2];
rate2_stim = B + alpha*beta*[0,1,2,0,1,2,0,1,2] + alpha*(1-beta)*[0,0,0,1,1,1,2,2,2];
rateShare  = gamma_flat*[1 1 1 1 1 1 1 1 1];
if length(nRbins) > 1
    
    edges_spikes1_all = cell(1, length(nRbins));
    edges_spikes2_all = cell(1, length(nRbins));
    for binIdx = 1:length(nRbins)
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
        [~, edges_spikes1_all{binIdx}] = binning({X1}, bin_opts);
        [~, edges_spikes2_all{binIdx}] = binning({X2}, bin_opts);
        edges_spikes1_all{binIdx}(1) = 0;
        edges_spikes2_all{binIdx}(1) = 0;
        x = nRbins(binIdx)+1;
        edges_spikes1_all{binIdx}(x) = Inf;
        edges_spikes2_all{binIdx}(x) = Inf;
    end
    edges_spikes1 = edges_spikes1_all;
    edges_spikes2 = edges_spikes2_all;
else     
    bin_opts.n_bins = {nRbins};
    [~, edges_spikes1] = binning({X1}, bin_opts);
    [~, edges_spikes2] = binning({X2}, bin_opts);
    edges_spikes1(1) = 0;
    edges_spikes2(1) = 0;
    x = nRbins+1;
    edges_spikes1(x) = Inf;
    edges_spikes2(x) = Inf;
end

if length(nRbins)>1
    for binidx=1:length(nRbins)
        IdxRate = nRbins(binidx);
        rateShared = rateShare(1:IdxRate);
        rate1 = rate1_stim(1:IdxRate);
        rate2 = rate2_stim(1:IdxRate);
        num_stimuli = IdxRate;
        GroundTruth_value{binidx} = GroundTruth(edges_spikes1{binidx}, edges_spikes2{binidx}, rate1, rate2,rateShared, nRbins(binidx), num_stimuli, 100, redundancy_measure);
    end
else
    GroundTruth_value =  GroundTruth(edges_spikes1, edges_spikes2, rate1, rate2,rateShared, nRbins, num_stimuli, 100, redundancy_measure);
end
if strcmp(redundancy_measure, 'I_BROJA')
    if isalphasweep
        folder_name = 'Alpha_sweep';
    elseif isBinSweep
        folder_name = 'Bin_sweep';
    else
        folder_name = 'Broja';
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
    save(['Results/', folder_name, '/', 'GroundTruth_', char(simul_case), '_', char(info_amount), '_', num2str(alpha), '.mat']);
elseif isBinSweep
    save(['Results/', folder_name, '/', 'GroundTruth_', char(simul_case), '_', char(info_amount),'.mat']);
elseif strcmp(info_amount, 'alpha_sweep')
    save(['Results/', folder_name, '/', 'GroundTruth_', char(simul_case), '_', char(info_amount), '_', num2str(alpha), '.mat']);
else
    save(['Results/', folder_name, '/', 'GroundTruth_', char(simul_case), '_', char(info_amount), '.mat']);
end