clc; clear; close all;
ResultFolder = 'Alpha_Bin_sweep';
simul_cases = {'uncorr_unique','unique_high_red','unique_high_syn', 'bit_of_all'};
bias_corrections = {'qe', 'shuff'};
atom = 'Joint';
nCases = length(simul_cases);
% Load all data
nBins = 8;
PID_qe_cases = cell(nCases,nBins);
PID_shuff_cases = cell(nCases,nBins);
for c = 1:nCases
    for binIdx = 1:nBins
        simul_case = simul_cases{c};
        filename = sprintf('Results/%s/Simuldata_%s_bin_alpha_sweep_I_BROJA.mat', ResultFolder, simul_case);
       
            nameLoad_qe = 'PID_v_qe';
            nameLoad_shuff = 'PID_v_shuff';
            data_qe = load(filename, nameLoad_qe);
            data_shuff = load(filename, nameLoad_shuff);
            data_qe =  data_qe.(nameLoad_qe);
            data_shuff  = data_shuff.(nameLoad_shuff);
            PID_qe_cases{c,binIdx} = data_qe(:,:,:,:,binIdx);
            PID_shuff_cases{c,binIdx} = data_shuff(:,:,:,:,binIdx);
    end
end
% Sizes (assumed consistent across cases)
[nTrialsAvail, ~, nLevels] = size(PID_qe_cases{1,1});
% Prepare sigmoid
sigmoid = @(p, x) 1 ./ (1 + exp(-p(1)*(x - p(2))));
% Prepare outputs
weight_matrix = zeros(nBins, nTrialsAvail, 200);  % fine grid (adjust size if needed)
info_levels_fine = linspace(0, 1, 200);
p_fits = zeros(nBins, nTrialsAvail, 2);
for binIdx = 1:nBins
    max_info = log2(binIdx+1);
    for i = 1:nTrialsAvail
        all_info_fractions = [];
        all_w_qe_emp = [];
        % Collect all dots from all cases
        for c = 1:nCases
            % Case-specific GT and info fraction (normalize by max_info)
            GT_case = squeeze(PID_qe_cases{c,binIdx}(8,1,:));
            info_frac_case = GT_case / max_info;
            for j = 1:nLevels
                est_qe = PID_qe_cases{c,binIdx}(i,1,j);
                est_shuff = PID_shuff_cases{c,binIdx}(i,1,j);
                R_qe = abs(est_qe - GT_case(j));
                R_shuff = abs(est_shuff - GT_case(j));
                w_qe_emp = R_shuff / (R_qe + R_shuff);
                w_qe_emp = min(max(w_qe_emp, 0), 1);
                % Collect for fitting
                all_info_fractions(end+1) = info_frac_case(j);
                all_w_qe_emp(end+1) = w_qe_emp;
            end
        end
        % Fit sigmoid on pooled points for this trial row i
        loss = @(p) sum((sigmoid(p, all_info_fractions) - all_w_qe_emp).^2);
        p0 = [10, 0.5];
        p_fit = fminsearch(loss, p0);
        p_fits(binIdx,i,:) = p_fit;
        % Generate weights on fine grid
        weight_matrix(binIdx,i,:) = sigmoid(p_fit, info_levels_fine);
    end
end
%%
save('adaptive_weight_matrix_allcases.mat', 'weight_matrix', 'info_levels_fine', 'p_fits');