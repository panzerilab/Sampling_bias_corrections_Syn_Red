addpath(genpath('Functions'))

%% Figure 1
% Run script zero_info.py to generate the histograms
% run("Functions/PlotFigure1AC.m")
close all;
PlotFigure1AC()
%% Figure 2
close all;
PlotFigure2discr()
PlotFigure2gauss()

%% Figure 3
close all;
PlotFigure3()
%% Figure 4
close all;
%discrete
trial_categories = [16 32 64 128 256 512 1024 2048];
PlotFigure4(trial_categories, 'bit_of_all', 'low', ...
            'I_BROJA', 'Broja', [], 'discr');
set(gcf, 'Renderer', 'painters');
outDir = 'Figures_mat';
if ~exist(outDir, 'dir')
    mkdir(outDir);
end
exportgraphics(gcf, fullfile(outDir, 'Figure_4_discrete.svg'), 'ContentType','vector');
%gauss
info_amount = 3;
PlotFigure4(trial_categories, 'bit_of_all', info_amount, ...
            'I_BROJA', 'Gauss', [], 'gauss');
set(gcf, 'Renderer', 'painters');
exportgraphics(gcf, fullfile(outDir, 'Figure_4_gauss.svg'), 'ContentType','vector');
%% Figure 5
close all;
outDir = 'Figures_mat';
figure_handle = figure('Units', 'centimeters', 'Position', [1, 1, 18, 12]); %'Position',[1 1 18 6]
Gauss_Figure = tiledlayout(figure_handle, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
worst2 = [];
worst_suff2 = [];
subplot_idx = 1;
for idxRow = 1:2
    if idxRow == 1
        % ---------- Discrete (Alpha_Bin_sweep) ----------
        ResultFolder = 'Bin_sweep';
        categories = 2:9;
        bias_correction = {'qe','shuff', 'weighted'}; %qeShuff weighted 'plugin',
        simul_cases = {'bit_of_all'};
        info_amount = 'low';
        redundancy_measure = 'I_BROJA';
        max_v = 0.3;  min_v = -0.4;
        isGauss = false; isUnion = false;
        trialIdx = 2;  % your original setting
    else
        % ------------- Gaussian (dimension sweep) -------------
        ResultFolder = 'Gauss';
        categories = 4:8:84; %4:2:42;                          % dimensions
        bias_correction = {'resample','shuffle','shuff-resamp','Venkatesh'}; %'plugin',
        simul_cases = {''};
        info_amount = 3;                               % as in your code
        redundancy_measure = '';
        % y-lims: plugin panel gets a different range in your original;
        % but here we use a common range for all lines within each panel:
        max_v = 10;  min_v = -4;
        isGauss = true; isUnion = false;
        trialIdx = 256;
    end
    [w, ws] = PlotFigure5( ...
        trialIdx, categories, simul_cases{1}, bias_correction, info_amount, redundancy_measure, ...
        subplot_idx, Gauss_Figure, 2, 3, ResultFolder, min_v, max_v, isGauss, isUnion,false);
    subplot_idx = subplot_idx + 3;  
end
set(gcf, 'Renderer', 'painters');
outDir = 'Figures_mat';
if ~exist(outDir, 'dir')
    mkdir(outDir);
end
exportgraphics(gcf, fullfile(outDir, 'Figure_5.svg'), 'ContentType','vector');



%% Figure 6
run("Real_data_analyses/PlotFigure6.m")

%% Figure 7
run("Functions/PlotFigure7.m")

%% Figure 8 
run("Functions/PlotFigure8.m")

%% Supplementary figures
%% Figure S2
PlotFigureS2discr()
PlotFigureS2gauss()
%% Figure S3
% AlphaSweep
close all, clc, clear
ResultFolder = 'Alpha_sweep';
dir_to_save = 'Figures/';
simul_case = {'bit_of_all'}; %unique_high_red, unique_high_syn
bias_correction = {'plugin','qe','shuff', 'qeShuff','infoCorr'};
Atoms = {'Joint', 'Synergy', 'Redundancy'};
trial_categories = [16,32,64,128,256,512,1024,2048];
alpha_values = 1:15;
redundancy_measure = 'I_BROJA';
nRbins = 4;
info_amount = {'alpha_sweep'};

figure_handle = figure('Units', 'inches', 'Position', [1, 1,15, 8]);
DiscreteInfosweep = tiledlayout(figure_handle, 3, 5, 'TileSpacing', 'compact', 'Padding', 'compact');
subplot_idx = 1;
for atomIdx = 1:length(Atoms)
    atom_v = Atoms{atomIdx};
    if atomIdx == 1
        min_v = 0;
        max_v = 1.4;
    elseif atomIdx == 2
        min_v = -0.15;
        max_v = 0.4;
    elseif atomIdx == 3
        min_v = -0.05;
        max_v = 0.55;
    elseif atomIdx == 4
        min_v = -0.05;
        max_v = 0.4;
    end

    PlotFigureS3(trial_categories, 0, atom_v, simul_case{1}, ...
        bias_correction, info_amount, redundancy_measure, ...
        subplot_idx, DiscreteInfosweep, ResultFolder, max_v, min_v)
    subplot_idx = subplot_idx + 5;
end
outDir = 'Figures_mat';
if ~exist(outDir, 'dir')
    mkdir(outDir);
end
exportgraphics(gcf, fullfile(outDir, 'Figure_S3.svg'), 'ContentType','vector');
%% Figure S4
cfg = struct;
cfg.ResultFolder = 'Broja';
cfg.mode = 'discrete';
cfg.simul_cases = {'bit_of_all','uncorr_unique','unique_high_red','unique_high_syn'};
cfg.info_types = {'Joint','Syn','Red'};
cfg.bias_correction = {'plugin','qe','shuff','qeShuff','infoCorr'};
cfg.info_amount = 'low';
cfg.redundancy_measure = 'I_BROJA';
cfg.trial_categories = [16 32 64 128 256 512 1024 2048];
cfg.row_names = {'Bit of all','Uncorr','High Red','High Syn'};
cfg.row_ylim = [-0.1 1.1; -0.1 1.8; -0.1 1.2; -0.1 1.4];
cfg.save_path = 'Figures_mat/Figure_S4_low.svg';
cfg.thumb_rule_x = 256;
pid_plot_grid(cfg);
cfg.info_amount = 'high';
cfg.save_path = 'Figures_mat/Figure_S4_high.svg';
cfg.thumb_rule_x = 256;
pid_plot_grid(cfg);


%% Figure S5
cfg = struct;
cfg.ResultFolder = 'Imin';
cfg.mode = 'discrete';
cfg.simul_cases = {'bit_of_all','uncorr_unique','unique_high_red','unique_high_syn'};
cfg.info_types = {'Joint','Syn','Red'};
cfg.bias_correction = {'plugin','qe','shuff','qeShuff','infoCorr'};
cfg.colors = [0 0.447 0.741; 0.85 0.325 0.098; 0.929 0.694 0.125; 0.494 0.184 0.556; 0.466 0.674 0.188];
cfg.info_amount = 'low';
cfg.redundancy_measure = 'I_min';
cfg.trial_categories = [16 32 64 128 256 512 1024 2048];
cfg.row_names = {'Bit of all','Uncorr','High Red','High Syn'};
cfg.row_ylim = [-0.1 1.1; -0.1 1.8; -0.1 1.2; -0.1 1.4];
cfg.save_path = 'Figures_mat/Figure_S5_low.svg';
cfg.thumb_rule_x = 256;
pid_plot_grid(cfg);
cfg.info_amount = 'high';
cfg.save_path = 'Figures_mat/Figure_S5_high.svg';
cfg.thumb_rule_x = 256;
pid_plot_grid(cfg);


%% Figure S6
cfg = struct;
cfg.ResultFolder = 'IMMI';
cfg.mode = 'discrete';
cfg.simul_cases = {'bit_of_all','uncorr_unique','unique_high_red','unique_high_syn'};
cfg.info_types = {'Joint','Syn','Red'};
cfg.bias_correction = {'plugin','qe','shuff','qeShuff','infoCorr'};
cfg.info_amount = 'low';
cfg.redundancy_measure = 'I_MMI';
cfg.trial_categories = [16 32 64 128 256 512 1024 2048];
cfg.row_names = {'Bit of all','Uncorr','High Red','High Syn'};
cfg.row_ylim = [-0.1 1.1; -0.1 1.8; -0.1 1.2; -0.1 1.4];
cfg.save_path = 'Figures_mat/Figure_S6_low.svg';
cfg.thumb_rule_x = 256;
pid_plot_grid(cfg);
cfg.info_amount = 'high';
cfg.save_path = 'Figures_mat/Figure_S6_high.svg';
cfg.thumb_rule_x = 256;
pid_plot_grid(cfg);


%% Figure S7
PlotFigureS7('discr', 'Alphas', 1:15, 'ResultFolder','Broja', 'SavePrefix','Figures_mat/Figure_S7_discrete');
PlotFigureS7('gauss', 'Alphas', [1 1.1 1.2 1.5 2 10000], 'ResultFolder','Gauss', 'SavePrefix','Figures_mat/Figure_S7_gauss');

%% Figure S8
% Figure InfoSweep Gaussian across dim
% AlphaSweep
close all, clc, clear
ResultFolder = 'Gauss';
simul_case = '';
bias_correction = {'plugin','resample','shuffle', 'Venkatesh', 'shuff-resamp'};
Atoms = {'Joint', 'Synergy', 'Redundancy'};
trial_categories = [256];
redundancy_measure = '';
nRbins = 4;
info_amount = {'alpha_sweep'};
%categories = 4:2:40;
categories = 4:8:80;
figure_handle = figure('Units', 'inches', 'Position', [1, 1, 15, 12]);
Fig = tiledlayout(figure_handle, 3, 5, 'TileSpacing', 'compact', 'Padding', 'compact');
subplot_idx = 1;
for atomIdx = 1:length(Atoms)
    atom_v = Atoms{atomIdx};
    if atomIdx == 1
        min_v = 0;
        max_v = 45;
    elseif atomIdx == 2
        min_v = 0;
        max_v = 42;
    elseif atomIdx == 3
        min_v = 0;
        max_v = 15;
    elseif atomIdx == 4
        min_v = 0;
        max_v = 15;
    end
    PlotFigureS8(trial_categories, categories, atom_v, simul_case, bias_correction, info_amount, redundancy_measure, subplot_idx, Fig, ResultFolder, max_v, min_v, true)
    subplot_idx = subplot_idx + 5;
end
dir_to_save = 'Figures_mat';
fileName = 'Figure_S8.svg';
filePath = fullfile(dir_to_save, fileName);
saveas(Fig, filePath, 'svg');

%% Figure S9
% Figure InfoSweep Gaussian across trials
% AlphaSweep
clc, clear, close all
ResultFolder = 'Gauss';
dir_to_save = 'Figures/';
simul_case = '';
bias_correction = {'plugin','resample','shuffle', 'Venkatesh', 'shuff-resample'};
Atoms = {'Joint', 'Synergy', 'Redundancy'};
trial_categories = [64,128,256,512,1024,2058];
redundancy_measure = '';
nRbins = 4;
info_amount = {'alpha_sweep'};

figure_handle = figure('Units', 'inches', 'Position', [1, 1, 14, 6]);
Fig = tiledlayout(figure_handle, 3, 5, 'TileSpacing', 'compact', 'Padding', 'compact');
subplot_idx = 1;

for atomIdx = 1:numel(Atoms)
    atom_v = Atoms{atomIdx};

    if atomIdx == 1
        min_v = 0;    max_v = 30;
    elseif atomIdx == 2
        min_v = -1;   max_v = 15;
    elseif atomIdx == 3
        min_v = -1;   max_v = 7;
    elseif atomIdx == 4
        min_v = 0;    max_v = 6;
    end

    PlotFigureS9( ...
        trial_categories, ...
        20, ...
        atom_v, ...
        bias_correction, ...
        subplot_idx, ...
        Fig, ...
        ResultFolder, ...
        max_v, ...
        min_v);

    subplot_idx = subplot_idx + 5;
end

dir_to_save = 'Figures_mat';
fileName = 'Figure_S9.svg';
filePath = fullfile(dir_to_save, fileName);
saveas(Fig, filePath, 'svg');

%% Figure S10
cfg = struct;
cfg.ResultFolder = 'Gauss';
cfg.mode = 'gauss';
cfg.redundancy_measure = 'gauss';
cfg.simul_cases = {'bit_of_all','high_syn','zero_syn','both_unq'};
cfg.info_types = {'Joint','Syn','Red'};
cfg.bias_correction = {'plugin','resample','shuff','shuff-resample','venkatesh'};
cfg.colors = [0 0.447 0.741; 0.85 0.325 0.098; 0.929 0.694 0.125; 0.466 0.674 0.188; 0.494 0.184 0.556];
cfg.info_amount = 3;       % i index for the info dimension in Gauss files
cfg.gauss_dim = 20;        % choose M
cfg.trial_categories = [16 32 64 128 256 512 1024 2048];
cfg.row_names = {'Bit of all','High Syn','Zero Syn','Both Unique'};
cfg.row_ylim = [0 15; 0 15; 0 15; 0 15];
cfg.thumb_rule_x = 128;     % optional vertical guide
cfg.save_path = 'Figures_mat/Figure_S10_low.svg';
pid_plot_grid(cfg);

cfg = struct;
cfg.ResultFolder = 'Gauss';
cfg.redundancy_measure = 'gauss';
cfg.mode = 'gauss';
cfg.simul_cases = {'bit_of_all','high_syn','zero_syn','both_unq'};
cfg.info_types = {'Joint','Syn','Red'};
cfg.bias_correction = {'plugin','resample','shuff','shuff-resample','venkatesh'};
cfg.colors = [0 0.447 0.741; 0.85 0.325 0.098; 0.929 0.694 0.125; 0.466 0.674 0.188; 0.494 0.184 0.556];
cfg.info_amount = 1;       % i index for the info dimension in Gauss files
cfg.gauss_dim = 20;        % choose M
cfg.trial_categories = [16 32 64 128 256 512 1024 2048];
cfg.row_names = {'Bit of all','High Syn','Zero Syn','Both Unique'};
cfg.row_ylim = [0 20; 0 30; 0 20; 0 20];
cfg.thumb_rule_x = 128;     % optional vertical guide
cfg.save_path = 'Figures_mat/Figure_S10_high.svg';
pid_plot_grid(cfg);

%% Figure S11
clc, clear, close all;
run("Functions/PlotFigureS11.m")

%% Figure S12
clc, clear, close all;
run('Figure_S12/panel_AB/panelAB.m')
run('Figure_S12/pnel_CDE/Simulation_DiscretVsGaussian_2Stim.m')
run('Figure_S12/panel_CDE/Simulation_DiscretVsGaussian_4Stim.m')
run('Figure_S12/panel_CDE/Simulation_DiscretVsGaussian_12Stim.m')
run('Figure_S12/panel_F/Plot_GPIDHighInfo.m')
%% Figure S13
clc, clear, close all;
run("Functions/PlotFigureS13.m")

%% Figure S14
clc, clear, close all;
PlotFigureS14()

%% Figure S15
clc, clear, close all;
run("Real_data_analyses/PlotFigureS15.m")

%% Figure S16
clc, clear, close all;
PlotFigureS16('Full')
PlotFigureS16('256tp')

%% Figure S17
clc, clear, close all;
run("Functions/PlotFigureS17.m")
