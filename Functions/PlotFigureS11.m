clc; clear; close all;

folder = 'Alpha_sweep';
simul_case = 'bit_of_all';
redundancy_measure = 'I_BROJA';
datafile = sprintf('Results/%s/Simuldata_%s_alpha_sweep_%s100.mat', ...
                   folder, simul_case, redundancy_measure);

%% --- Load PID variants from the big .mat file ---
variants = {'qeShuff','shuff','qe'};
S = struct();
for v = 1:numel(variants)
    nameLoad = sprintf('PID_v_%s', variants{v});
    tmp = load(datafile, nameLoad);
    S.(nameLoad) = tmp.(nameLoad);
end
PID_v_qeShuff = S.PID_v_qeShuff;
PID_v_shuff   = S.PID_v_shuff;
PID_v_qe      = S.PID_v_qe;

%% --- Trial categories ---
trial_categories = [16, 32, 64, 128, 256, 512, 1024, 2048] * 4;

max_info = log2(4);
infoLevels = size(PID_v_shuff,4);

%% --- Build optimally weighted PID and load ground truth ---
PID_v_weighted = zeros(size(PID_v_shuff));

GT_Joint = zeros(1,infoLevels);
GT_Syn   = zeros(1,infoLevels);
GT_Red   = zeros(1,infoLevels);

% Load weighting matrix once
W = load('adaptive_weight_matrix_allcases.mat', 'weight_matrix', 'info_levels_fine');
weight_matrix    = W.weight_matrix;
info_levels_fine = W.info_levels_fine;

for infolevel = 1:infoLevels
    PID_shuff = squeeze(PID_v_shuff(:,:,:,infolevel));
    PID_qe    = squeeze(PID_v_qe(:,:,:,infolevel));

    PID_qeshuff = (PID_qe(:,1,:) + PID_shuff(:,1,:)) / 2;
    infoval_Joint = squeeze(mean(PID_qeshuff,3));
    info_frac_matrix = infoval_Joint / max_info;

    PID_weighted = zeros(size(PID_shuff));
    for t = 1:size(PID_shuff,1)
        val = info_frac_matrix(t);
        [~, idx] = min(abs(info_levels_fine - val));
        w = weight_matrix(3,t,idx);

        QE_slice    = squeeze(PID_qe(t,:,:));
        Shuff_slice = squeeze(PID_shuff(t,:,:));
        PID_weighted(t,:,:) = w * QE_slice + (1 - w) * Shuff_slice;
    end
    PID_v_weighted(:,:,:,infolevel) = PID_weighted;

    gtfile = sprintf('Results/Broja/GroundTruth_%s_%d.mat', simul_case, infolevel);
    GT = load(gtfile, 'GroundTruth_value');
    GT = GT.GroundTruth_value;

    GT_Joint(infolevel) = sum(GT(1:4));
    GT_Syn(infolevel)   = GT(4);
    GT_Red(infolevel)   = GT(1);
end

%% --- Extract components: [nTrials x nReps x nInfoLevels]
Joint_equal = squeeze(PID_v_qeShuff(:,1,:,:));
Syn_equal   = squeeze(PID_v_qeShuff(:,2,:,:));
Red_equal   = squeeze(PID_v_qeShuff(:,3,:,:));

Joint_opt = squeeze(PID_v_weighted(:,1,:,:));
Syn_opt   = squeeze(PID_v_weighted(:,2,:,:));
Red_opt   = squeeze(PID_v_weighted(:,3,:,:));

nTrials     = size(Joint_equal,1);
nInfoLevels = size(Joint_equal,3);

%% --- Residual bias: subtract GT (broadcast over trials,reps) ---
GT_Joint_mat = reshape(GT_Joint, [1 1 nInfoLevels]);
GT_Syn_mat   = reshape(GT_Syn,   [1 1 nInfoLevels]);
GT_Red_mat   = reshape(GT_Red,   [1 1 nInfoLevels]);
GT_Joint(1);
GT_Joint(15);


resJoint_equal = Joint_equal - GT_Joint_mat;
resSyn_equal   = Syn_equal   - GT_Syn_mat;
resRed_equal   = Red_equal   - GT_Red_mat;

resJoint_opt = Joint_opt - GT_Joint_mat;
resSyn_opt   = Syn_opt   - GT_Syn_mat;
resRed_opt   = Red_opt   - GT_Red_mat;

% Mean across repetitions (dim 2)
mean_resJoint_equal = squeeze(mean(resJoint_equal, 2, 'omitnan')); % [nTrials x nInfoLevels]
mean_resSyn_equal   = squeeze(mean(resSyn_equal,   2, 'omitnan'));
mean_resRed_equal   = squeeze(mean(resRed_equal,   2, 'omitnan'));

mean_resJoint_opt = squeeze(mean(resJoint_opt, 2, 'omitnan'));
mean_resSyn_opt   = squeeze(mean(resSyn_opt,   2, 'omitnan'));
mean_resRed_opt   = squeeze(mean(resRed_opt,   2, 'omitnan'));

%% --- DEFINE COLORMAPS ---
numColors = 256;

yellowColors = [254 227 145; 254 196 79; 254 153 41; 217 95 14; 153 52 4] / 255;
JointColormap = interp1(linspace(0,1,size(yellowColors,1)), yellowColors, linspace(0,1,numColors));

greenColors = [229 245 224; 161 217 155; 116 196 118; 49 163 84; 0 109 44] / 255;
SynColormap = interp1(linspace(0,1,size(greenColors,1)), greenColors, linspace(0,1,numColors));

blueColors = [198 219 239; 158 202 225; 107 174 214; 66 146 198; 33 113 181; 8 69 148] / 255;
RedColormap = interp1(linspace(0,1,size(blueColors,1)), blueColors, linspace(0,1,numColors));

%% --- Plot config ---
metrics    = {'Joint', 'Synergy', 'Redundancy'};
equal_data = {mean_resJoint_equal, mean_resSyn_equal, mean_resRed_equal};
opt_data   = {mean_resJoint_opt,   mean_resSyn_opt,   mean_resRed_opt};
colormaps  = {JointColormap, SynColormap, RedColormap};

figure('Position',[100 100 700 900], 'Color', 'w');

% Global y limits across all panels
all_vals = [];
for m = 1:3
    all_vals = [all_vals; equal_data{m}(:); opt_data{m}(:)];
end
ymax  = max(abs(all_vals), [], 'omitnan');
ylims = [-ymax ymax];

% --- Log ticks & minor ticks (between decades) ---
xmin = min(trial_categories);
xmax = max(trial_categories);

pmin = floor(log10(xmin));
pmax = ceil(log10(xmax));

% Major ticks: 10^p
xticks_pow = 10.^(pmin:pmax);
xticks_pow = xticks_pow(xticks_pow >= xmin & xticks_pow <= xmax);
xticklabels_pow = arrayfun(@(x) sprintf('10^{%d}', round(log10(x))), xticks_pow, 'UniformOutput', false);

% Minor ticks: 2..9 * 10^p (unlabeled)
minorVals = [];
for p = pmin:(pmax-1)
    minorVals = [minorVals, (2:9) * 10^p];
end
minorVals = minorVals(minorVals >= xmin & minorVals <= xmax);

%% --- Plotting ---
for m = 1:3
    data_equal = equal_data{m};
    data_opt   = opt_data{m};
    cmap       = colormaps{m};

    % infolevel index increases with info: 1=low info -> light, end=high info -> dark
    colorIndices = round(linspace(1, numColors, nInfoLevels));  % NO flip

    % ===== Equal =====
    subplot(3,2,(m-1)*2+1); hold on
    for info = 1:nInfoLevels
        plot(trial_categories, data_equal(:,info), ...
            'Color', cmap(colorIndices(info),:), 'LineWidth', 1.5);
    end
    set(gca,'XScale','log','TickLabelInterpreter','tex','FontSize',11)
    ax = gca;
    ax.XTick = xticks_pow;
    ax.XTickLabel = xticklabels_pow;
    ax.XMinorTick = 'on';
    ax.XAxis.MinorTickValues = minorVals;

    xlim([xmin xmax]); ylim(ylims);
    box off; grid off;
    yline(0,'k--','LineWidth',0.8);
    ylabel('Residual bias [bits]','FontWeight','bold')
    title(sprintf('%s – Equal weighting', metrics{m}),'FontWeight','bold')
    if m == 3, xlabel('Trials','FontWeight','bold'); end

    % ===== Optimal =====
    subplot(3,2,(m-1)*2+2); hold on
    for info = 1:nInfoLevels
        plot(trial_categories, data_opt(:,info), ...
            'Color', cmap(colorIndices(info),:), 'LineWidth', 1.5);
    end
    set(gca,'XScale','log','TickLabelInterpreter','tex','FontSize',11)
    ax = gca;
    ax.XTick = xticks_pow;
    ax.XTickLabel = xticklabels_pow;
    ax.XMinorTick = 'on';
    ax.XAxis.MinorTickValues = minorVals;

    xlim([xmin xmax]); ylim(ylims);
    box off; grid off;
    xline(512,'k--','LineWidth',0.8); 
    yline(0,'k--','LineWidth',0.8);
    title(sprintf('%s – Optimal weighting', metrics{m}),'FontWeight','bold')
    if m == 3, xlabel('Trials','FontWeight','bold'); end
end



%% --- Grey legend: light=Low info (bottom), dark=High info (top) ---
axL = axes('Position',[0.93 0.15 0.02 0.7]);
greyColors = [245 245 245; 200 200 200; 150 150 150; 100 100 100; 60 60 60; 20 20 20] / 255;
colormap(axL, greyColors);
cb = colorbar(axL, 'Ticks',[0 1], 'TickLabels',{'0.01','0.98'}, 'FontSize',11);
ylabel(cb, 'Information level', 'FontWeight','bold');
axis(axL,'off')

%% --- Save as SVG ---
outDir = fullfile(fileparts(pwd), 'Figures_mat');  

if ~exist(outDir, 'dir')
    mkdir(outDir);
end

exportgraphics(gcf, fullfile(outDir, 'Figure_S11.svg'), 'ContentType','vector');
