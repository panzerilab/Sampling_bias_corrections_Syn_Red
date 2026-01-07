clc, clear close all;
rng(0);
dir_to_save = 'Results';
if ~exist(dir_to_save, 'dir')
    mkdir(dir_to_save);
end
%% Simulation Fig 1 
% The simulation is already done in the plotting function, run that instead

%% Simulation Fig 2A, 3, 4A, S4, S5, S6
simul_cases = {'uncorr_unique','unique_high_red','unique_high_syn', 'bit_of_all'};
redundancy_measures = {'I_BROJA','I_min','I_MMI'};
info_vals = {'low', 'high'};
trial_categories = [16,32,64,128,256,512,1024,2048];
n_iter = 100;
nRbins = 4;

for redMeas = 1:length(redundancy_measures)
    redundancy_measure = redundancy_measures{redMeas};
    for simul = 1:length(simul_cases)
        simul_case = simul_cases{simul};
        for infoLevel = 1:length(info_vals)
            info_amount = info_vals{infoLevel};

            if infoLevel == 1
                alpha = 7;
            elseif infoLevel == 2
                alpha = 10;
            end    
            runSimulation(nRbins, trial_categories, n_iter, simul_case, alpha, info_amount, false, false, redundancy_measure);
            computeGroundTruth(nRbins, simul_case, alpha, info_amount, false, false, redundancy_measure)
        end 
    end
end 

%% Simulation Fig 5A, S2ABC
dir_to_save = 'Results';
if ~exist(dir_to_save, 'dir')
    mkdir(dir_to_save);
end
redundancy_measure = 'I_BROJA';
n_iter = 100;

trial_categories = [256, 512, 1024];
alpha = [7];
info_amount = {'low'};
simul_case = {'bit_of_all'};
nRbins = [2:9];
computeGroundTruth(nRbins, simul_case{1}, alpha, info_amount, false, true, redundancy_measure)
runSimulation(nRbins, trial_categories, n_iter, simul_case{1}, alpha, info_amount, false, true, redundancy_measure);
alpha = [10];
info_amount = {'high'};
simul_case = {'bit_of_all'};
nRbins = [2:9];
computeGroundTruth(nRbins, simul_case{1}, alpha, info_amount, false, true, redundancy_measure)
runSimulation(nRbins, trial_categories, n_iter, simul_case{1}, alpha, info_amount, false, true, redundancy_measure);

%% Simulation Fig S3, S7A, S11
simul_case = {'bit_of_all','uncorr_unique','unique_high_red','unique_high_syn'};
trial_categories = [16,32,64,128,256,512,1024,2048];
alpha_values = [1:15 20];
n_iter = 100;
redundancy_measure = 'I_BROJA';
nRbins = 4;
info_amount = 'alpha_sweep';
for simul = 1:length(simul_case)
        runSimulation(nRbins, trial_categories, n_iter, simul_case{simul},  alpha_values, info_amount, true, false, redundancy_measure);
end
 
for simul = 1:length(simul_case)
    for a = 1:length(alpha_values)
        runSimulation(nRbins, trial_categories, n_iter, simul_case{simul},  alpha_values(a), info_amount, false, false, redundancy_measure);
        computeGroundTruth(nRbins, simul_case, alpha_values(a), info_amount, false, false, redundancy_measure)
    end
end

%% Simulation Fig S16
simul_case = 'uncorr_unique';
trial_categories = [16,32,64,128,256,512,1024,2048];
alpha_values = 1;
n_iter = 100;
redundancy_measure = 'I_BROJA';
nRbins = 4;
alpha=1;
info_amount = {'low'};
runSimulation_SVM(nRbins, trial_categories, n_iter, simul_case,  alpha_values, info_amount, 10, redundancy_measure);

%% 3 Source PID Simulation discrete
run("Functions/three_sources_nonzeroinfo.m")
