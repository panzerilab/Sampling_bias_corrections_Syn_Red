clc, clear, close all;

set(0, 'DefaultTextFontName', 'Arial');
set(0, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultUicontrolFontName', 'Arial');
set(0, 'DefaultUitableFontName', 'Arial');
set(0, 'DefaultUipanelFontName', 'Arial');
set(0, 'DefaultLegendFontName', 'Arial');

customColors_info_PID = [
    0.9290, 0.6940, 0.1250;  % Joint
    0, 0.4470, 0.7410;       % Red
    0.4660, 0.6740, 0.1880;  % Syn
    0.5, 0.5, 0.5;           % Unq
];

chosenAtom = 'Red';

% -------- weight rule (edit here to change behavior) --------
get_weight = @(f) min(max(f,0),1);  % linear clamp in [0,1]

%% ======================== Prepare A1 data ========================
filename = ['Results/Real_data_analysis/A1_PID.mat'];
load(filename);

plugin_all=[]; qe_all=[]; qeshuffSub_all=[]; shuffSub_all=[];
for dataSet = 1:length(Results_plugin)
    if ~isempty(Results_plugin(dataSet).HalfData)
        plugin_all     = [plugin_all    ; Results_plugin(dataSet).HalfData.Result];
        qe_all         = [qe_all        ; Results_qe(dataSet).HalfData.Result];
        qeshuffSub_all = [qeshuffSub_all; Results_qeshuffSub(dataSet).HalfData.Result];
        shuffSub_all   = [shuffSub_all  ; Results_shuffSub(dataSet).HalfData.Result];
    end
end
HalfData_plugin_all     = squeeze(mean(plugin_all     ,2));
HalfData_qe_all         = squeeze(mean(qe_all         ,2));
HalfData_qeshuffSub_all = squeeze(mean(qeshuffSub_all ,2));
HalfData_shuffSub_all   = squeeze(mean(shuffSub_all   ,2));

plugin_all=[]; qe_all=[]; qeshuffSub_all=[]; shuffSub_all=[];
for dataSet = 1:length(Results_plugin)
    if ~isempty(Results_plugin(dataSet).FullData)
        plugin_all     = [plugin_all    ; Results_plugin(dataSet).FullData.Result];
        qe_all         = [qe_all        ; Results_qe(dataSet).FullData.Result];
        qeshuffSub_all = [qeshuffSub_all; Results_qeshuffSub(dataSet).FullData.Result];
        shuffSub_all   = [shuffSub_all  ; Results_shuffSub(dataSet).FullData.Result];
    end
end
FullData_plugin_all     = plugin_all;
FullData_qe_all         = qe_all;
FullData_qeshuffSub_all = qeshuffSub_all;
FullData_shuffSub_all   = shuffSub_all;

% --------- NEW: build adaptive-weighted blends (A1) ----------
max_info_half_A1 = log2(2); %max( (HalfData_qe_all(:,1) + HalfData_shuffSub_all(:,1))/2 );
max_info_full_A1 = log2(2); %max( (FullData_qe_all(:,1) + FullData_shuffSub_all(:,1))/2 );
[HalfData_weighted_all_A1,  w_half_A1] = apply_adaptive_weighting(HalfData_qe_all,  HalfData_shuffSub_all,  max_info_half_A1, get_weight);
[FullData_weighted_all_A1,  w_full_A1] = apply_adaptive_weighting(FullData_qe_all,  FullData_shuffSub_all,  max_info_full_A1, get_weight);

% t-tests (Full vs Half) for each component order 1→(1,3,2)
for k = 1:3
    if     k==1, i=1;
    elseif k==2, i=3;
    else         i=2;
    end
    [~, P_joint_A1_plugin(k)]   = ttest(FullData_plugin_all(:,i),   HalfData_plugin_all(:,i));
    [~, P_joint_A1_qe(k)]       = ttest(FullData_qe_all(:,i),       HalfData_qe_all(:,i));
    [~, P_joint_A1_shuffSub(k)] = ttest(FullData_shuffSub_all(:,i), HalfData_shuffSub_all(:,i));
    % replaced qeshuffSub with weighted:
    [~, P_joint_A1_weighted(k)] = ttest(FullData_weighted_all_A1(:,i), HalfData_weighted_all_A1(:,i));
end

% -------- A1 --------
meanpluginA1 = [ mean(FullData_plugin_all(:,1)), mean(FullData_plugin_all(:,3)), mean(FullData_plugin_all(:,2)), ...
                 mean(HalfData_plugin_all(:,1)), mean(HalfData_plugin_all(:,3)), mean(HalfData_plugin_all(:,2)) ];
sempluginA1  = [ std(FullData_plugin_all(:,1))/sqrt(size(FullData_plugin_all,1)), ...
                 std(FullData_plugin_all(:,3))/sqrt(size(FullData_plugin_all,1)), ...
                 std(FullData_plugin_all(:,2))/sqrt(size(FullData_plugin_all,1)), ...
                 std(HalfData_plugin_all(:,1))/sqrt(size(HalfData_plugin_all,1)), ...
                 std(HalfData_plugin_all(:,3))/sqrt(size(HalfData_plugin_all,1)), ...
                 std(HalfData_plugin_all(:,2))/sqrt(size(HalfData_plugin_all,1)) ];

meanqeA1 = [ mean(FullData_qe_all(:,1)), mean(FullData_qe_all(:,3)), mean(FullData_qe_all(:,2)), ...
             mean(HalfData_qe_all(:,1)), mean(HalfData_qe_all(:,3)), mean(HalfData_qe_all(:,2)) ];
semqeA1  = [ std(FullData_qe_all(:,1))/sqrt(size(FullData_qe_all,1)), ...
             std(FullData_qe_all(:,3))/sqrt(size(FullData_qe_all,1)), ...
             std(FullData_qe_all(:,2))/sqrt(size(FullData_qe_all,1)), ...
             std(HalfData_qe_all(:,1))/sqrt(size(HalfData_qe_all,1)), ...
             std(HalfData_qe_all(:,3))/sqrt(size(HalfData_qe_all,1)), ...
             std(HalfData_qe_all(:,2))/sqrt(size(HalfData_qe_all,1)) ];

meanshuffSubA1 = [ mean(FullData_shuffSub_all(:,1)), mean(FullData_shuffSub_all(:,3)), mean(FullData_shuffSub_all(:,2)), ...
                   mean(HalfData_shuffSub_all(:,1)), mean(HalfData_shuffSub_all(:,3)), mean(HalfData_shuffSub_all(:,2)) ];
semshuffSubA1  = [ std(FullData_shuffSub_all(:,1))/sqrt(size(FullData_shuffSub_all,1)), ...
                   std(FullData_shuffSub_all(:,3))/sqrt(size(FullData_shuffSub_all,1)), ...
                   std(FullData_shuffSub_all(:,2))/sqrt(size(FullData_shuffSub_all,1)), ...
                   std(HalfData_shuffSub_all(:,1))/sqrt(size(HalfData_shuffSub_all,1)), ...
                   std(HalfData_shuffSub_all(:,3))/sqrt(size(HalfData_shuffSub_all,1)), ...
                   std(HalfData_shuffSub_all(:,2))/sqrt(size(HalfData_shuffSub_all,1)) ];

meanWeightedA1 = [ mean(FullData_weighted_all_A1(:,1)), mean(FullData_weighted_all_A1(:,3)), mean(FullData_weighted_all_A1(:,2)), ...
                   mean(HalfData_weighted_all_A1(:,1)), mean(HalfData_weighted_all_A1(:,3)), mean(HalfData_weighted_all_A1(:,2)) ];
semWeightedA1  = [ std(FullData_weighted_all_A1(:,1))/sqrt(size(FullData_weighted_all_A1,1)), ...
                   std(FullData_weighted_all_A1(:,3))/sqrt(size(FullData_weighted_all_A1,1)), ...
                   std(FullData_weighted_all_A1(:,2))/sqrt(size(FullData_weighted_all_A1,1)), ...
                   std(HalfData_weighted_all_A1(:,1))/sqrt(size(HalfData_weighted_all_A1,1)), ...
                   std(HalfData_weighted_all_A1(:,3))/sqrt(size(HalfData_weighted_all_A1,1)), ...
                   std(HalfData_weighted_all_A1(:,2))/sqrt(size(HalfData_weighted_all_A1,1)) ];


%% ======================== Prepare CA1 data ========================
filename = ['Results/Real_data_analysis/CA1_PID.mat'];
load(filename);

plugin_all=[]; qe_all=[]; qeshuffSub_all=[]; shuffSub_all=[];
for dataSet = 1:length(Results_plugin)
    if ~isempty(Results_plugin(dataSet).HalfData)
        plugin_all     = [plugin_all    ; Results_plugin(dataSet).HalfData.Result];
        qe_all         = [qe_all        ; Results_qe(dataSet).HalfData.Result];
        qeshuffSub_all = [qeshuffSub_all; Results_qeshuffSub(dataSet).HalfData.Result];
        shuffSub_all   = [shuffSub_all  ; Results_shuffSub(dataSet).HalfData.Result];
    end
end
HalfData_plugin_all     = squeeze(mean(plugin_all     ,2));
HalfData_qe_all         = squeeze(mean(qe_all         ,2));
HalfData_qeshuffSub_all = squeeze(mean(qeshuffSub_all ,2));
HalfData_shuffSub_all   = squeeze(mean(shuffSub_all   ,2));

plugin_all=[]; qe_all=[]; qeshuffSub_all=[]; shuffSub_all=[];
for dataSet = 1:length(Results_plugin)
    if ~isempty(Results_plugin(dataSet).FullData)
        plugin_all     = [plugin_all    ; Results_plugin(dataSet).FullData.Result];
        qe_all         = [qe_all        ; Results_qe(dataSet).FullData.Result];
        qeshuffSub_all = [qeshuffSub_all; Results_qeshuffSub(dataSet).FullData.Result];
        shuffSub_all   = [shuffSub_all  ; Results_shuffSub(dataSet).FullData.Result];
    end
end
FullData_plugin_all     = plugin_all;
FullData_qe_all         = qe_all;
FullData_qeshuffSub_all = qeshuffSub_all;
FullData_shuffSub_all   = shuffSub_all;

% --------- NEW: adaptive-weighted (CA1) ----------
max_info_half_CA1 = log2(4); %max( (HalfData_qe_all(:,1) + HalfData_shuffSub_all(:,1))/2 );
max_info_full_CA1 = log2(4); %max( (FullData_qe_all(:,1) + FullData_shuffSub_all(:,1))/2 );
[HalfData_weighted_all_CA1,  w_half_CA1] = apply_adaptive_weighting(HalfData_qe_all,  HalfData_shuffSub_all,  max_info_half_CA1, get_weight);
[FullData_weighted_all_CA1,  w_full_CA1] = apply_adaptive_weighting(FullData_qe_all,  FullData_shuffSub_all,  max_info_full_CA1, get_weight);

for k = 1:3
    if     k==1, i=1;
    elseif k==2, i=3;
    else         i=2;
    end
    [~, P_joint_CA1_plugin(k)]   = ttest(FullData_plugin_all(:,i),   HalfData_plugin_all(:,i));
    [~, P_joint_CA1_qe(k)]       = ttest(FullData_qe_all(:,i),       HalfData_qe_all(:,i));
    [~, P_joint_CA1_shuffSub(k)] = ttest(FullData_shuffSub_all(:,i), HalfData_shuffSub_all(:,i));
    [~, P_joint_CA1_weighted(k)] = ttest(FullData_weighted_all_CA1(:,i), HalfData_weighted_all_CA1(:,i));
end

% -------- CA1 --------
meanpluginCA1 = [ mean(FullData_plugin_all(:,1)), mean(FullData_plugin_all(:,3)), mean(FullData_plugin_all(:,2)), ...
                  mean(HalfData_plugin_all(:,1)), mean(HalfData_plugin_all(:,3)), mean(HalfData_plugin_all(:,2)) ];
sempluginCA1  = [ std(FullData_plugin_all(:,1))/sqrt(size(FullData_plugin_all,1)), ...
                  std(FullData_plugin_all(:,3))/sqrt(size(FullData_plugin_all,1)), ...
                  std(FullData_plugin_all(:,2))/sqrt(size(FullData_plugin_all,1)), ...
                  std(HalfData_plugin_all(:,1))/sqrt(size(HalfData_plugin_all,1)), ...
                  std(HalfData_plugin_all(:,3))/sqrt(size(HalfData_plugin_all,1)), ...
                  std(HalfData_plugin_all(:,2))/sqrt(size(HalfData_plugin_all,1)) ];

meanqeCA1 = [ mean(FullData_qe_all(:,1)), mean(FullData_qe_all(:,3)), mean(FullData_qe_all(:,2)), ...
              mean(HalfData_qe_all(:,1)), mean(HalfData_qe_all(:,3)), mean(HalfData_qe_all(:,2)) ];
semqeCA1  = [ std(FullData_qe_all(:,1))/sqrt(size(FullData_qe_all,1)), ...
              std(FullData_qe_all(:,3))/sqrt(size(FullData_qe_all,1)), ...
              std(FullData_qe_all(:,2))/sqrt(size(FullData_qe_all,1)), ...
              std(HalfData_qe_all(:,1))/sqrt(size(HalfData_qe_all,1)), ...
              std(HalfData_qe_all(:,3))/sqrt(size(HalfData_qe_all,1)), ...
              std(HalfData_qe_all(:,2))/sqrt(size(HalfData_qe_all,1)) ];

meanshuffSubCA1 = [ mean(FullData_shuffSub_all(:,1)), mean(FullData_shuffSub_all(:,3)), mean(FullData_shuffSub_all(:,2)), ...
                    mean(HalfData_shuffSub_all(:,1)), mean(HalfData_shuffSub_all(:,3)), mean(HalfData_shuffSub_all(:,2)) ];
semshuffSubCA1  = [ std(FullData_shuffSub_all(:,1))/sqrt(size(FullData_shuffSub_all,1)), ...
                    std(FullData_shuffSub_all(:,3))/sqrt(size(FullData_shuffSub_all,1)), ...
                    std(FullData_shuffSub_all(:,2))/sqrt(size(FullData_shuffSub_all,1)), ...
                    std(HalfData_shuffSub_all(:,1))/sqrt(size(HalfData_shuffSub_all,1)), ...
                    std(HalfData_shuffSub_all(:,3))/sqrt(size(HalfData_shuffSub_all,1)), ...
                    std(HalfData_shuffSub_all(:,2))/sqrt(size(HalfData_shuffSub_all,1)) ];

meanWeightedCA1 = [ mean(FullData_weighted_all_CA1(:,1)), mean(FullData_weighted_all_CA1(:,3)), mean(FullData_weighted_all_CA1(:,2)), ...
                    mean(HalfData_weighted_all_CA1(:,1)), mean(HalfData_weighted_all_CA1(:,3)), mean(HalfData_weighted_all_CA1(:,2)) ];
semWeightedCA1  = [ std(FullData_weighted_all_CA1(:,1))/sqrt(size(FullData_weighted_all_CA1,1)), ...
                    std(FullData_weighted_all_CA1(:,3))/sqrt(size(FullData_weighted_all_CA1,1)), ...
                    std(FullData_weighted_all_CA1(:,2))/sqrt(size(FullData_weighted_all_CA1,1)), ...
                    std(HalfData_weighted_all_CA1(:,1))/sqrt(size(HalfData_weighted_all_CA1,1)), ...
                    std(HalfData_weighted_all_CA1(:,3))/sqrt(size(HalfData_weighted_all_CA1,1)), ...
                    std(HalfData_weighted_all_CA1(:,2))/sqrt(size(HalfData_weighted_all_CA1,1)) ];


%% ==================== Prepare A1_monkey data =====================
filename = ['Results/Real_data_analysis/KayserA1_PID.mat'];
load(filename);
plugin_all=[]; qe_all=[]; qeshuffSub_all=[]; shuffSub_all=[];
for dataSet = 1:length(Results_plugin)
    if ~isempty(Results_plugin(dataSet).HalfData)
        plugin_all     = [plugin_all    ; Results_plugin(dataSet).HalfData.Result];
        qe_all         = [qe_all        ; Results_qe(dataSet).HalfData.Result];
        qeshuffSub_all = [qeshuffSub_all; Results_qeshuffSub(dataSet).HalfData.Result];
        shuffSub_all   = [shuffSub_all  ; Results_shuffSub(dataSet).HalfData.Result];
    end
end
HalfData_plugin_all     = squeeze(mean(plugin_all     ,2));
HalfData_qe_all         = squeeze(mean(qe_all         ,2));
HalfData_qeshuffSub_all = squeeze(mean(qeshuffSub_all ,2));
HalfData_shuffSub_all   = squeeze(mean(shuffSub_all   ,2));

load(filename);
plugin_all=[]; qe_all=[]; qeshuffSub_all=[]; shuffSub_all=[];
for dataSet = 1:length(Results_plugin)
    if ~isempty(Results_plugin(dataSet).FullData)
        plugin_all     = [plugin_all    ; Results_plugin(dataSet).FullData.Result];
        qe_all         = [qe_all        ; Results_qe(dataSet).FullData.Result];
        qeshuffSub_all = [qeshuffSub_all; Results_qeshuffSub(dataSet).FullData.Result];
        shuffSub_all   = [shuffSub_all  ; Results_shuffSub(dataSet).FullData.Result];
    end
end
FullData_plugin_all     = plugin_all;
FullData_qe_all         = qe_all;
FullData_qeshuffSub_all = qeshuffSub_all;
FullData_shuffSub_all   = shuffSub_all;

% --------- NEW: adaptive-weighted (A1_monkey) ----------
max_info_half_A1m = log2(2); %max( (HalfData_qe_all(:,1) + HalfData_shuffSub_all(:,1))/2 );
max_info_full_A1m = log2(2); %max( (FullData_qe_all(:,1) + FullData_shuffSub_all(:,1))/2 );
[HalfData_weighted_all_A1m,  w_half_A1m] = apply_adaptive_weighting(HalfData_qe_all,  HalfData_shuffSub_all,  max_info_half_A1m, get_weight);
[FullData_weighted_all_A1m,  w_full_A1m] = apply_adaptive_weighting(FullData_qe_all,  FullData_shuffSub_all,  max_info_full_A1m, get_weight);

for k = 1:3
    if     k==1, i=1;
    elseif k==2, i=3;
    else         i=2;
    end
    [~, P_joint_A1_m_plugin(k)]   = ttest(FullData_plugin_all(:,i),   HalfData_plugin_all(:,i));
    [~, P_joint_A1_m_qe(k)]       = ttest(FullData_qe_all(:,i),       HalfData_qe_all(:,i));
    [~, P_joint_A1_m_shuffSub(k)] = ttest(FullData_shuffSub_all(:,i), HalfData_shuffSub_all(:,i));
    [~, P_joint_A1_m_weighted(k)] = ttest(FullData_weighted_all_A1m(:,i), HalfData_weighted_all_A1m(:,i));
end

% -------- A1_monkey --------
meanpluginA1_monkey = [ mean(FullData_plugin_all(:,1)), mean(FullData_plugin_all(:,3)), mean(FullData_plugin_all(:,2)), ...
                        mean(HalfData_plugin_all(:,1)), mean(HalfData_plugin_all(:,3)), mean(HalfData_plugin_all(:,2)) ];
sempluginA1_monkey  = [ std(FullData_plugin_all(:,1))/sqrt(size(FullData_plugin_all,1)), ...
                        std(FullData_plugin_all(:,3))/sqrt(size(FullData_plugin_all,1)), ...
                        std(FullData_plugin_all(:,2))/sqrt(size(FullData_plugin_all,1)), ...
                        std(HalfData_plugin_all(:,1))/sqrt(size(HalfData_plugin_all,1)), ...
                        std(HalfData_plugin_all(:,3))/sqrt(size(HalfData_plugin_all,1)), ...
                        std(HalfData_plugin_all(:,2))/sqrt(size(HalfData_plugin_all,1)) ];

meanqeA1_monkey = [ mean(FullData_qe_all(:,1)), mean(FullData_qe_all(:,3)), mean(FullData_qe_all(:,2)), ...
                    mean(HalfData_qe_all(:,1)), mean(HalfData_qe_all(:,3)), mean(HalfData_qe_all(:,2)) ];
semqeA1_monkey  = [ std(FullData_qe_all(:,1))/sqrt(size(FullData_qe_all,1)), ...
                    std(FullData_qe_all(:,3))/sqrt(size(FullData_qe_all,1)), ...
                    std(FullData_qe_all(:,2))/sqrt(size(FullData_qe_all,1)), ...
                    std(HalfData_qe_all(:,1))/sqrt(size(HalfData_qe_all,1)), ...
                    std(HalfData_qe_all(:,3))/sqrt(size(HalfData_qe_all,1)), ...
                    std(HalfData_qe_all(:,2))/sqrt(size(HalfData_qe_all,1)) ];

meanshuffSubA1_monkey = [ mean(FullData_shuffSub_all(:,1)), mean(FullData_shuffSub_all(:,3)), mean(FullData_shuffSub_all(:,2)), ...
                          mean(HalfData_shuffSub_all(:,1)), mean(HalfData_shuffSub_all(:,3)), mean(HalfData_shuffSub_all(:,2)) ];
semshuffSubA1_monkey  = [ std(FullData_shuffSub_all(:,1))/sqrt(size(FullData_shuffSub_all,1)), ...
                          std(FullData_shuffSub_all(:,3))/sqrt(size(FullData_shuffSub_all,1)), ...
                          std(FullData_shuffSub_all(:,2))/sqrt(size(FullData_shuffSub_all,1)), ...
                          std(HalfData_shuffSub_all(:,1))/sqrt(size(HalfData_shuffSub_all,1)), ...
                          std(HalfData_shuffSub_all(:,3))/sqrt(size(HalfData_shuffSub_all,1)), ...
                          std(HalfData_shuffSub_all(:,2))/sqrt(size(HalfData_shuffSub_all,1)) ];

meanWeightedA1_monkey = [ mean(FullData_weighted_all_A1m(:,1)), mean(FullData_weighted_all_A1m(:,3)), mean(FullData_weighted_all_A1m(:,2)), ...
                          mean(HalfData_weighted_all_A1m(:,1)), mean(HalfData_weighted_all_A1m(:,3)), mean(HalfData_weighted_all_A1m(:,2)) ];
semWeightedA1_monkey  = [ std(FullData_weighted_all_A1m(:,1))/sqrt(size(FullData_weighted_all_A1m,1)), ...
                          std(FullData_weighted_all_A1m(:,3))/sqrt(size(FullData_weighted_all_A1m,1)), ...
                          std(FullData_weighted_all_A1m(:,2))/sqrt(size(FullData_weighted_all_A1m,1)), ...
                          std(HalfData_weighted_all_A1m(:,1))/sqrt(size(HalfData_weighted_all_A1m,1)), ...
                          std(HalfData_weighted_all_A1m(:,3))/sqrt(size(HalfData_weighted_all_A1m,1)), ...
                          std(HalfData_weighted_all_A1m(:,2))/sqrt(size(HalfData_weighted_all_A1m,1)) ];


%% ============================= Plot =============================
figure('Units', 'centimeters', 'Position', [1, 1, 12, 12]);
tiledlayout(3,4,"TileSpacing","compact");

% -------- Row: A1 --------
maxV = 0.3;
subplotIdx = 1;

plotPIDterms(meanpluginA1, sempluginA1, subplotIdx, 'plugin', 0, maxV, false);
hold on; % dadd_pvalue(P_joint_A1_plugin, [1 3 5], [2 4 6], 0.005, 0.92*maxV);
subplotIdx = subplotIdx + 1;
xticklabels([]);
yticks([0, 0.1, 0.2, 0.3]);

plotPIDterms(meanqeA1, semqeA1, subplotIdx, 'qe', 0, maxV, false);
hold on; %add_pvalue(P_joint_A1_qe, [1 3 5], [2 4 6], 0.005, 0.92*maxV);
subplotIdx = subplotIdx + 1;
xticklabels([]);
ylabel([]);
yticklabels([]);

plotPIDterms(meanshuffSubA1, semshuffSubA1, subplotIdx, 'shuffSub', 0, maxV, false);
hold on; %add_pvalue(P_joint_A1_shuffSub, [1 3 5], [2 4 6], 0.005, 0.92*maxV);
subplotIdx = subplotIdx + 1;
xticklabels([]);
ylabel([]);
yticklabels([]);

plotPIDterms(meanWeightedA1, semWeightedA1, subplotIdx, 'merged', 0, maxV, false);
hold on; %add_pvalue(P_joint_A1_weighted, [1 3 5], [2 4 6], 0.005, 0.92*maxV);
xticklabels([]);
ylabel([]);
yticklabels([]);

% -------- Row: CA1 --------
maxV = 0.15;
subplotIdx = subplotIdx + 1;

plotPIDterms(meanpluginCA1, sempluginCA1, subplotIdx, 'plugin', 0, maxV, false);
hold on; %add_pvalue(P_joint_CA1_plugin, [1 3 5], [2 4 6], 0.003, 0.92*maxV);
subplotIdx = subplotIdx + 1;
xticklabels([]);
title([]);
yticks([0, 0.05, 0.1, 0.15]);

plotPIDterms(meanqeCA1, semqeCA1, subplotIdx, 'qe', 0, maxV, false);
hold on; %add_pvalue(P_joint_CA1_qe, [1 3 5], [2 4 6], 0.003, 0.92*maxV);
subplotIdx = subplotIdx + 1;
xticklabels([]);
title([]);
ylabel([]);
yticklabels([]);

plotPIDterms(meanshuffSubCA1, semshuffSubCA1, subplotIdx, 'shuffSub', 0, maxV, false);
hold on; %add_pvalue(P_joint_CA1_shuffSub, [1 3 5], [2 4 6], 0.003, 0.92*maxV);
subplotIdx = subplotIdx + 1;
xticklabels([]);
title([]);
ylabel([]);
yticklabels([]);

plotPIDterms(meanWeightedCA1, semWeightedCA1, subplotIdx, 'merged', 0, maxV, false);
hold on; %add_pvalue(P_joint_CA1_weighted, [1 3 5], [2 4 6], 0.003, 0.92*maxV);
xticklabels([]);
title([]);
ylabel([]);
yticklabels([]);

% -------- Row: A1_monkey --------
maxV = 0.1; minV = 0;
subplotIdx = subplotIdx + 1;

plotPIDterms(meanpluginA1_monkey, sempluginA1_monkey, subplotIdx, 'plugin', minV, maxV, false);
hold on; %add_pvalue(P_joint_A1_m_plugin, [1 3 5], [2 4 6], 0.002, 0.92*maxV);
subplotIdx = subplotIdx + 1;
title([]);
yticks([0, 0.05, 0.1]);

plotPIDterms(meanqeA1_monkey, semqeA1_monkey, subplotIdx, 'qe', minV, maxV, false);
hold on; %add_pvalue(P_joint_A1_m_qe, [1 3 5], [2 4 6], 0.002, 0.92*maxV);
subplotIdx = subplotIdx + 1;
title([]);
ylabel([]);
yticklabels([]);

plotPIDterms(meanshuffSubA1_monkey, semshuffSubA1_monkey, subplotIdx, 'shuffSub', minV, maxV, false);
hold on; %add_pvalue(P_joint_A1_m_shuffSub, [1 3 5], [2 4 6], 0.002, 0.92*maxV);
subplotIdx = subplotIdx + 1;
title([]);
ylabel([]);
yticklabels([]);

plotPIDterms(meanWeightedA1_monkey, semWeightedA1_monkey, subplotIdx, 'merged', minV, maxV, false);
hold on; %add_pvalue(P_joint_A1_m_weighted, [1 3 5], [2 4 6], 0.002, 0.92*maxV);
title([]);
ylabel([]);
yticklabels([]);

filename = ['Figures/Figure_S13_' chosenAtom '.svg'];
saveas(gcf, filename);
filename = ['Figures/Figure_S13_' chosenAtom '.png'];
print(filename, '-dpng', '-r300');

%% ===== Transposed Markdown Table: Only Red & Syn, Δ and Normalized Δ =====
methods = {'plugin', 'qe', 'shuffSub', 'Weighted'};  % <-- replaced qeshuffSub
areas = {'A1', 'CA1', 'A1_monkey'};
components = {'Red', 'Syn'};
component_indices = [2, 3];

absTable = {};  normTable = {};

header = {'Component'};
for a = 1:length(areas)
    for m = 1:length(methods)
        header{end+1} = sprintf('%s-%s', areas{a}, methods{m});
    end
end

for c = 1:length(components)
    compRowAbs = {components{c}};
    compRowNorm = {components{c}};
    for a = 1:length(areas)
        for m = 1:length(methods)
            meanVar = eval(sprintf('mean%s%s', methods{m}, areas{a}));
            fullVal = meanVar(component_indices(c)*2 - 1);
            halfVal = meanVar(component_indices(c)*2);
            delta = halfVal - fullVal;
            normDelta = delta / fullVal;
            compRowAbs{end+1}  = sprintf('%.4f',  delta);
            compRowNorm{end+1} = sprintf('%.2f%%', 100*normDelta);
        end
    end
    absTable  = [absTable ; compRowAbs ];
    normTable = [normTable; compRowNorm];
end

fprintf('### Table RXXa: Absolute Δ(Half - Full) in bits\n');
fprintf('| %s |\n', strjoin(header, ' | '));
fprintf('|%s|\n', repmat(' --- |', 1, numel(header)));
for i = 1:size(absTable,1),  fprintf('| %s |\n', strjoin(absTable(i,:), ' | ')); end

fprintf('\n### Table RXXb: Δ Normalized to Full (%%)\n');
fprintf('| %s |\n', strjoin(header, ' | '));
fprintf('|%s|\n', repmat(' --- |', 1, numel(header)));
for i = 1:size(normTable,1), fprintf('| %s |\n', strjoin(normTable(i,:), ' | ')); end

%% ===== helper plotting & weighting =====
function ax = plotPIDterms(meanVals, semVals, subplotIdx, plotTitle, minLim, maxLim, isPercent)
    % meanVals order expected:
    % [F-Joint, F-Syn, F-Red, H-Joint, H-Syn, H-Red]
    phaseLabels = {'F-Joint','F-Syn','F-Red','H-Joint','H-Syn','H-Red'};
    baseColors = [0.9290, 0.6940, 0.1250;   % Joint
                  0,       0.4470, 0.7410; % Red (used for "Red")
                  0.4660,  0.6740, 0.1880];% Syn

    % Build 6×3 colors: first 3 = Joint,Syn,Red; next 3 repeat Joint,Syn,Red
    c6 = [baseColors([1 3 2],:); baseColors([1 3 2],:)];

    ax = nexttile(subplotIdx);
    b = bar(meanVals, 'FaceColor', 'flat', 'EdgeColor', 'k');
    b.CData = c6;
    hold on;
    errorbar(1:numel(meanVals), meanVals, semVals, 'k', 'LineStyle', 'none', 'LineWidth', 1.5);
    xticks(1:numel(meanVals)); xticklabels(phaseLabels);
    if isPercent, ylabel('% of Joint Info'); else, ylabel('Information [bits]'); end
    title(plotTitle, 'Interpreter','none');
    ylim([minLim, maxLim]); grid off; box off; hold off;
end

function add_pvalue(p_values, x1, x2, text_offset, line_heights)
    % broadcast scalars to the length of p_values
    n = numel(p_values);
    if isscalar(x1),           x1 = repmat(x1, 1, n);           end
    if isscalar(x2),           x2 = repmat(x2, 1, n);           end
    if isscalar(line_heights), line_heights = repmat(line_heights, 1, n); end

    for i = 1:n
        stars = get_stars(p_values(i));
        line([x1(i), x2(i)], [line_heights(i), line_heights(i)], 'Color', 'k', 'LineWidth', 1);
        text((x1(i)+x2(i))/2, line_heights(i)+text_offset, stars, ...
             'HorizontalAlignment','center','Color','k','FontSize',8);
    end
end


function stars = get_stars(p_value)
if p_value < 0.001, stars = '***';
elseif p_value < 0.01, stars = '**';
elseif p_value < 0.05, stars = '*';
else, stars = 'n.s';
end
end

function [weighted_mat, w_used] = apply_adaptive_weighting(qe_mat, shuff_mat, max_info, get_weight)
    if isempty(qe_mat)
        weighted_mat = qe_mat; w_used = []; return;
    end
    joint_qe    = qe_mat(:,1);
    joint_shuff = shuff_mat(:,1);
    joint_avg   = (joint_qe + joint_shuff)/2;

    info_frac = joint_avg ./ max_info;
    w_used = arrayfun(get_weight, info_frac);

    weighted_mat = w_used .* qe_mat + (1 - w_used) .* shuff_mat;
end
