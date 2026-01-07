clear, close all
%% Groundtruth Simulation
clc, clear, close all
% Groundtruth Simulation
rng(0);
num_stimuli=2;
alphasweep= 1:10;
Bsweep=     1:10;
trialsSweep =  10000;
num_trials = trialsSweep(1);
gamma_flat = 0;
X1 = zeros(length(alphasweep), length(Bsweep), num_stimuli*num_trials);
X2 = zeros(length(alphasweep), length(Bsweep), num_stimuli*num_trials);
Y  = zeros(1,num_stimuli*num_trials);
I_gt = zeros(length(alphasweep), length(Bsweep));
for a=1:length(alphasweep)
    alpha=alphasweep(a);
    disp(a)
    parfor b=1:length(Bsweep)
        B= Bsweep(b);
        beta = 0;
        
        rate1 = B + alpha*[1,0];
        rate2 = B + alpha*[0,1];
        X1_tmp=[];
        X2_tmp=[];
        for rateIdx = 1:length(rate1)
            X1_tmp =[X1_tmp, poissrnd(rate1(rateIdx),1,num_trials)];
            X2_tmp =[X2_tmp, poissrnd(rate2(rateIdx),1,num_trials)];
        end
        X1(a,b,:) = X1_tmp;
        X2(a,b,:) = X2_tmp;
        I_gt(a,b) = groundtruthinformation(rate1, rate2, 1/num_stimuli, num_stimuli, 100);
    end
end
Y =[ones(1,num_trials) 2*ones(1,num_trials)];

datapath = fullfile('Discrete','SimulationR1','R1A_B');
if isfolder(datapath) == false
    % Create directory
    mkdir datapath % This part needs modification
end

save(fullfile(datapath,'I_gt_2Stim.mat'), "I_gt" );
save(fullfile(datapath,'data_2Stim.mat'), 'X1', 'X2', 'Y', 'Bsweep', 'alphasweep');
%% Analysis Simulation
load('data_2stim.mat');
opts_MI.bias='plugin';
opts_MI.bin_method={'eqpop', 'none'};
opts_MI.n_bins={4};
opts_MI.supressWarnings = true;
opts_MI.computeNulldist = false;
I_discrete = zeros(length(alphasweep), length(Bsweep));
for a = 1:length(alphasweep)
    disp(a) 
    parfor b = 1:length(Bsweep)
        X1_ab = (squeeze(X1(a,b,:)))';
        X2_ab = (squeeze(X2(a,b,:)))';
        MI_v = MI({[X1_ab;X2_ab], Y}, {'I(A;B)'}, opts_MI);
        I_discrete(a,b)=MI_v{1};
    end
end
% save('Discrete/SimulationR1/R1A_B/I_discrete_2S.mat', 'I_discrete');


datapath = fullfile('Discrete','SimulationR1','R1A_B');
if isfolder(datapath) == false
    % Create directory
    mkdir datapath % This part needs modification
end

save(fullfile(datapath,'I_discrete_2S.mat'), "I_discrete" );
%% Plot
clc, clear, close all
greyColors = [ ...
    255, 255, 255;  
    220, 220, 220;
    190, 190, 200;
    175, 175, 190;
    160, 160, 180;
    145, 145, 170;
    130, 130, 160;
    115, 115, 150;
    100, 100, 140;
    85, 85, 130;
    70, 70, 120;
    55, 55, 110;
    40, 40, 100;
] / 255;
numColors = 256; 
greyColormap = interp1(linspace(0, 1, size(greyColors, 1)), greyColors, linspace(0, 1, numColors));

blueColors = [ ...
    255, 255, 255;  
    240, 245, 250;  
    220, 230, 245;  
    200, 215, 240;  
    180, 200, 235;  
    160, 185, 230;  
    140, 170, 225;  
    120, 155, 220;  
    100, 140, 215;  
    80, 125, 210;   
    60, 110, 205; 
    40, 95, 200;
] / 255;
numColors = 256; 
blueColormap = interp1(linspace(0, 1, size(blueColors, 1)), blueColors, linspace(0, 1, numColors));


customColormap1 = blueColormap;
customColormap2 = greyColormap;

% Plot
load('Discrete/SimulationR1/R1A_B/I_discrete_2S.mat', 'I_discrete');
info_venk = load('info_Gpid_2Stim.mat');
load('I_gt_2Stim.mat');
load('data_2Stim.mat');

fontSize = 19;
min_value = 0;
max_value = 1;

min_value_1 = 0;
max_value_1 = 1;

alphasweep = 1:10;
B = 1:10;
alpha = 1:length(alphasweep);

I_noBinning_ratio = I_discrete ./ I_gt;
I_venka_ratio = info_venk.info_FullGauss ./ I_gt;

figure_handle = figure('Units', 'inches', 'Position', [1, 1, 24, 4]);
fontname("arial");
tcl = tiledlayout(1, 5, 'TileSpacing', 'tight');

ax = nexttile(tcl, 1);
h(1) = heatmap(Bsweep(B), alphasweep(alpha), I_gt(alpha ,B), 'Colormap',customColormap2);
clim([min_value_1, max_value_1]);
h(1).YDisplayData = flipud(h(1).YDisplayData);
title('Ground Truth');
ylabel('alpha' );
xlabel('B');
h(1).FontSize = fontSize ;
set(struct(h(1)).Axes, 'FontSize', fontSize );
set(struct(h(1)).Axes.Title, 'FontSize', fontSize + 2);
set(struct(h(1)).Axes.XLabel, 'FontSize', fontSize + 2);
set(struct(h(1)).Axes.YLabel, 'FontSize', fontSize + 2);
set(struct(h(1)).Axes.Title, 'FontName', 'Arial');
set(struct(h(1)).Axes.XLabel, 'FontName', 'Arial');
set(struct(h(1)).Axes.YLabel, 'FontName', 'Arial');
xLabels = strings(size(B));
xLabels(1) = "1";
xLabels(5) = "5";
xLabels(10) = "10";

yLabels = strings(size(alpha));
yLabels(1) = "10";
yLabels(5) = "5";
yLabels(10) = "1";
h(1).XDisplayLabels = xLabels;
h(1).YDisplayLabels = yLabels;

ax = nexttile(tcl, 2);
h(2) = heatmap(Bsweep(B), alphasweep(alpha), I_discrete(alpha,B), 'Colormap', customColormap2);
clim([min_value_1, max_value_1]);
h(2).YDisplayData = flipud(h(2).YDisplayData);
title('Discrete');
xlabel('B');
h(2).FontSize = fontSize ;
set(struct(h(2)).Axes, 'FontSize', fontSize );
set(struct(h(2)).Axes.Title, 'FontSize', fontSize + 2);
set(struct(h(2)).Axes.XLabel, 'FontSize', fontSize + 2);
set(struct(h(2)).Axes.YLabel, 'FontSize', fontSize + 2);
set(struct(h(2)).Axes.Title, 'FontName', 'Arial');
set(struct(h(2)).Axes.XLabel, 'FontName', 'Arial');
set(struct(h(2)).Axes.YLabel, 'FontName', 'Arial');
h(2).XDisplayLabels = xLabels;
h(2).YDisplayLabels = yLabels;


ax = nexttile(tcl, 3);
h(3) = heatmap(Bsweep(B), alphasweep(alpha), info_venk.info_FullGauss(alpha,B), 'Colormap', customColormap2);
clim([min_value_1, max_value_1]);

h(3).YDisplayData = flipud(h(3).YDisplayData);
title('Full Gaussian');
xlabel('B');
h(3).FontSize = fontSize ;
set(struct(h(3)).Axes, 'FontSize', fontSize );
set(struct(h(3)).Axes.Title, 'FontSize', fontSize + 2);
set(struct(h(3)).Axes.XLabel, 'FontSize', fontSize + 2);
set(struct(h(3)).Axes.YLabel, 'FontSize', fontSize + 2);
set(struct(h(3)).Axes.Title, 'FontName', 'Arial');
set(struct(h(3)).Axes.XLabel, 'FontName', 'Arial');
set(struct(h(3)).Axes.YLabel, 'FontName', 'Arial');
h(3).XDisplayLabels = xLabels;
h(3).YDisplayLabels = yLabels;

ax = nexttile(tcl, 4);
h(4) = heatmap(Bsweep(B), alphasweep(alpha), I_noBinning_ratio(alpha,B), 'Colormap', customColormap1);
clim([min_value, max_value]);
h(4).YDisplayData = flipud(h(4).YDisplayData);
title(sprintf('Discrete'));
xlabel('B');
h(4).FontSize = fontSize ;
set(struct(h(4)).Axes, 'FontSize', fontSize );
set(struct(h(4)).Axes.Title, 'FontSize', fontSize + 2);
set(struct(h(4)).Axes.XLabel, 'FontSize', fontSize + 2);
set(struct(h(4)).Axes.YLabel, 'FontSize', fontSize + 2);
set(struct(h(4)).Axes.Title, 'FontName', 'Arial');
set(struct(h(4)).Axes.XLabel, 'FontName', 'Arial');
set(struct(h(4)).Axes.YLabel, 'FontName', 'Arial');
h(4).XDisplayLabels = xLabels;
h(4).YDisplayLabels = yLabels;

ax = nexttile(tcl, 5);
h(5) = heatmap(Bsweep(B), alphasweep(alpha), I_venka_ratio(alpha,B), 'Colormap', customColormap1);
clim([min_value, max_value]);
h(5).YDisplayData = flipud(h(5).YDisplayData);
title(sprintf('Full Gaussian'));
xlabel('B');
h(5).FontSize = fontSize ;
set(struct(h(5)).Axes, 'FontSize', fontSize );
set(struct(h(5)).Axes.Title, 'FontSize', fontSize + 2);
set(struct(h(5)).Axes.XLabel, 'FontSize', fontSize + 2);
set(struct(h(5)).Axes.YLabel, 'FontSize', fontSize + 2);
set(struct(h(5)).Axes.Title, 'FontName', 'Arial');
set(struct(h(5)).Axes.XLabel, 'FontName', 'Arial');
set(struct(h(5)).Axes.YLabel, 'FontName', 'Arial');
h(5).XDisplayLabels = xLabels;
h(5).YDisplayLabels = yLabels;

saveas(figure_handle, 'Combined_Info_Map_2Stim.svg');

%% Helper Function
function I = groundtruthinformation(rate1, rate2, p_s, num_stimuli, numX)
p_x1x2_s = zeros(numX, numX, num_stimuli);
p_x1x2 = zeros(numX, numX);
for s = 1:num_stimuli
    for x1 = 0:numX-1
        for x2 = 0:numX-1
            p_x1x2_s(x1+1, x2+1, s) = poisspdf(x1, rate1(s)) * poisspdf(x2, rate2(s));
        end
    end
end
for s = 1:num_stimuli
    p_x1x2 = p_x1x2 + p_x1x2_s(:, :, s) * p_s;
end
I = 0;
for s = 1:num_stimuli
    for x1 = 0:numX-1
        for x2 = 0:numX-1
            p_cond = p_x1x2_s(x1+1, x2+1, s);
            p_joint = p_x1x2(x1+1, x2+1);
            if p_cond > 0 && p_joint > 0
                I = I + p_s * p_cond * log2(p_cond / p_joint);
            end
        end
    end
end
end