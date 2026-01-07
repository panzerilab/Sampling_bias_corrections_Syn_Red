% clc, clear, close all;

set(0, 'DefaultTextFontName', 'Arial');
set(0, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultUicontrolFontName', 'Arial');
set(0, 'DefaultUitableFontName', 'Arial');
set(0, 'DefaultUipanelFontName', 'Arial');
set(0, 'DefaultLegendFontName', 'Arial');

customColors_info_PID = [
    0.9290, 0.6940, 0.1250;  % Color for 'Joint'
    0, 0.4470, 0.7410;       % Color for 'Red'
    0.4660, 0.6740, 0.1880;  % Color for 'Syn'
    0.5, 0.5, 0.5;           % Color for 'Unq'
];

% Joint, Syn, Red, Unq

chosenAtom = 'Syn';
%% Prepare A1 data
filename = ['../Results/Real_data_analysis/A1_PID_' chosenAtom '.mat'];
load(filename);

plugin_all      = [];
qe_all          = [];
qeshuffSub_all  = [];
shuffSub_all    = [];
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

plugin_all      = [];
qe_all          = [];
qeshuffSub_all  = [];
shuffSub_all    = [];
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



for k = 1:3
    if k == 1 
        i = 1;
    elseif k == 2
        i = 3; 
    elseif k == 3
        i = 2;
    end 
       
    [h, P_joint_A1_plugin(k), ci, stats] = ttest(FullData_plugin_all(:,i), HalfData_plugin_all(:,i));
    [h, P_joint_A1_qe(k), ci, stats] = ttest(FullData_qe_all(:,i), HalfData_qe_all(:,i));
    [h, P_joint_A1_shuffSub(k), ci, stats] = ttest(FullData_shuffSub_all(:,i), HalfData_shuffSub_all(:,i));
    [h, P_joint_A1_qeshuffSub(k), ci, stats] = ttest(FullData_qeshuffSub_all(:,i), HalfData_qeshuffSub_all(:,i));
end 

meanpluginA1 = [mean(FullData_plugin_all(:,1)), mean(HalfData_plugin_all(:,1)), ...
                mean(FullData_plugin_all(:,3)), mean(HalfData_plugin_all(:,3)),...
                mean(FullData_plugin_all(:,2)), mean(HalfData_plugin_all(:,2))];

sempluginA1 = [ std(FullData_plugin_all(:,1))/sqrt(size(FullData_plugin_all,1)), std(HalfData_plugin_all(:,1))/sqrt(size(HalfData_plugin_all,1)), ...
                std(FullData_plugin_all(:,3))/sqrt(size(FullData_plugin_all,1)), std(HalfData_plugin_all(:,3))/sqrt(size(HalfData_plugin_all,1)),...
                std(FullData_plugin_all(:,2))/sqrt(size(FullData_plugin_all,1)), std(HalfData_plugin_all(:,2))/sqrt(size(HalfData_plugin_all,1))];


meanqeA1 = [mean(FullData_qe_all(:,1)), mean(HalfData_qe_all(:,1)), ...
            mean(FullData_qe_all(:,3)), mean(HalfData_qe_all(:,3)),...
            mean(FullData_qe_all(:,2)), mean(HalfData_qe_all(:,2))];

semqeA1 = [std(FullData_qe_all(:,1))/sqrt(size(FullData_qe_all,1)), std(HalfData_qe_all(:,1))/sqrt(size(HalfData_qe_all,1)), ...
           std(FullData_qe_all(:,3))/sqrt(size(FullData_qe_all,1)), std(HalfData_qe_all(:,3))/sqrt(size(HalfData_qe_all,1)),...
           std(FullData_qe_all(:,2))/sqrt(size(FullData_qe_all,1)), std(HalfData_qe_all(:,2))/sqrt(size(HalfData_qe_all,1))];


meanshuffSubA1 = [mean(FullData_shuffSub_all(:,1)), mean(HalfData_shuffSub_all(:,1)), ...
                  mean(FullData_shuffSub_all(:,3)), mean(HalfData_shuffSub_all(:,3)),...
                  mean(FullData_shuffSub_all(:,2)), mean(HalfData_shuffSub_all(:,2))];

semshuffSubA1 = [std(FullData_shuffSub_all(:,1))/sqrt(size(FullData_shuffSub_all,1)), std(HalfData_shuffSub_all(:,1))/sqrt(size(HalfData_shuffSub_all,1)), ...
                 std(FullData_shuffSub_all(:,3))/sqrt(size(FullData_shuffSub_all,1)), std(HalfData_shuffSub_all(:,3))/sqrt(size(HalfData_shuffSub_all,1)),...
                 std(FullData_shuffSub_all(:,2))/sqrt(size(FullData_shuffSub_all,1)), std(HalfData_shuffSub_all(:,2))/sqrt(size(HalfData_shuffSub_all,1))];

meanqeshuffSubA1 = [mean(FullData_qeshuffSub_all(:,1)), mean(HalfData_qeshuffSub_all(:,1)), ...
                    mean(FullData_qeshuffSub_all(:,3)), mean(HalfData_qeshuffSub_all(:,3)),...
                    mean(FullData_qeshuffSub_all(:,2)), mean(HalfData_qeshuffSub_all(:,2))];

semqeshuffSubA1 = [std(FullData_qeshuffSub_all(:,1))/sqrt(size(FullData_qeshuffSub_all,1)), std(HalfData_qeshuffSub_all(:,1))/sqrt(size(HalfData_qeshuffSub_all,1)), ...
                   std(FullData_qeshuffSub_all(:,3))/sqrt(size(FullData_qeshuffSub_all,1)), std(HalfData_qeshuffSub_all(:,3))/sqrt(size(HalfData_qeshuffSub_all,1)),...
                   std(FullData_qeshuffSub_all(:,2))/sqrt(size(FullData_qeshuffSub_all,1)), std(HalfData_qeshuffSub_all(:,2))/sqrt(size(HalfData_qeshuffSub_all,1))];

%% Prepare CA1 data
filename = ['../Results/Real_data_analysis/CA1_PID_' chosenAtom '.mat'];
load(filename);

plugin_all      = [];
qe_all          = [];
qeshuffSub_all  = [];
shuffSub_all    = [];
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

plugin_all      = [];
qe_all          = [];
qeshuffSub_all  = [];
shuffSub_all    = [];
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

for k = 1:3
    if k == 1 
        i = 1;
    elseif k == 2
        i = 3; 
    elseif k == 3
        i = 2;
    end 
       
    [h, P_joint_CA1_plugin(k), ci, stats] = ttest(FullData_plugin_all(:,i), HalfData_plugin_all(:,i));
    [h, P_joint_CA1_qe(k), ci, stats] = ttest(FullData_qe_all(:,i), HalfData_qe_all(:,i));
    [h, P_joint_CA1_shuffSub(k), ci, stats] = ttest(FullData_shuffSub_all(:,i), HalfData_shuffSub_all(:,i));
    [h, P_joint_CA1_qeshuffSub(k), ci, stats] = ttest(FullData_qeshuffSub_all(:,i), HalfData_qeshuffSub_all(:,i));
end 

meanpluginCA1 = [mean(FullData_plugin_all(:,1)), mean(HalfData_plugin_all(:,1)), ...
                mean(FullData_plugin_all(:,3)), mean(HalfData_plugin_all(:,3)),...
                mean(FullData_plugin_all(:,2)), mean(HalfData_plugin_all(:,2))];

sempluginCA1 = [ std(FullData_plugin_all(:,1))/sqrt(size(FullData_plugin_all,1)), std(HalfData_plugin_all(:,1))/sqrt(size(HalfData_plugin_all,1)), ...
                std(FullData_plugin_all(:,3))/sqrt(size(FullData_plugin_all,1)), std(HalfData_plugin_all(:,3))/sqrt(size(HalfData_plugin_all,1)),...
                std(FullData_plugin_all(:,2))/sqrt(size(FullData_plugin_all,1)), std(HalfData_plugin_all(:,2))/sqrt(size(HalfData_plugin_all,1))];


meanqeCA1 = [mean(FullData_qe_all(:,1)), mean(HalfData_qe_all(:,1)), ...
            mean(FullData_qe_all(:,3)), mean(HalfData_qe_all(:,3)),...
            mean(FullData_qe_all(:,2)), mean(HalfData_qe_all(:,2))];

semqeCA1 = [std(FullData_qe_all(:,1))/sqrt(size(FullData_qe_all,1)), std(HalfData_qe_all(:,1))/sqrt(size(HalfData_qe_all,1)), ...
           std(FullData_qe_all(:,3))/sqrt(size(FullData_qe_all,1)), std(HalfData_qe_all(:,3))/sqrt(size(HalfData_qe_all,1)),...
           std(FullData_qe_all(:,2))/sqrt(size(FullData_qe_all,1)), std(HalfData_qe_all(:,2))/sqrt(size(HalfData_qe_all,1))];


meanshuffSubCA1 = [mean(FullData_shuffSub_all(:,1)), mean(HalfData_shuffSub_all(:,1)), ...
                  mean(FullData_shuffSub_all(:,3)), mean(HalfData_shuffSub_all(:,3)),...
                  mean(FullData_shuffSub_all(:,2)), mean(HalfData_shuffSub_all(:,2))];

semshuffSubCA1 = [std(FullData_shuffSub_all(:,1))/sqrt(size(FullData_shuffSub_all,1)), std(HalfData_shuffSub_all(:,1))/sqrt(size(HalfData_shuffSub_all,1)), ...
                 std(FullData_shuffSub_all(:,3))/sqrt(size(FullData_shuffSub_all,1)), std(HalfData_shuffSub_all(:,3))/sqrt(size(HalfData_shuffSub_all,1)),...
                 std(FullData_shuffSub_all(:,2))/sqrt(size(FullData_shuffSub_all,1)), std(HalfData_shuffSub_all(:,2))/sqrt(size(HalfData_shuffSub_all,1))];

meanqeshuffSubCA1 = [mean(FullData_qeshuffSub_all(:,1)), mean(HalfData_qeshuffSub_all(:,1)), ...
                    mean(FullData_qeshuffSub_all(:,3)), mean(HalfData_qeshuffSub_all(:,3)),...
                    mean(FullData_qeshuffSub_all(:,2)), mean(HalfData_qeshuffSub_all(:,2))];

semqeshuffSubCA1 = [std(FullData_qeshuffSub_all(:,1))/sqrt(size(FullData_qeshuffSub_all,1)), std(HalfData_qeshuffSub_all(:,1))/sqrt(size(HalfData_qeshuffSub_all,1)), ...
                   std(FullData_qeshuffSub_all(:,3))/sqrt(size(FullData_qeshuffSub_all,1)), std(HalfData_qeshuffSub_all(:,3))/sqrt(size(HalfData_qeshuffSub_all,1)),...
                   std(FullData_qeshuffSub_all(:,2))/sqrt(size(FullData_qeshuffSub_all,1)), std(HalfData_qeshuffSub_all(:,2))/sqrt(size(HalfData_qeshuffSub_all,1))];


%% Prepare A1_monkey data
filename = ['../Results/Real_data_analysis/KayserA1_PID_' chosenAtom '_half.mat'];
load(filename);
plugin_all      = [];
qe_all          = [];
qeshuffSub_all  = [];
shuffSub_all    = [];
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

filename = ['../Results/Real_data_analysis/KayserA1_PID_' chosenAtom '.mat'];
load(filename);
plugin_all      = [];
qe_all          = [];
qeshuffSub_all  = [];
shuffSub_all    = [];
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

for k = 1:3
    if k == 1 
        i = 1;
    elseif k == 2
        i = 3; 
    elseif k == 3
        i = 2;
    end 
       
    [h, P_joint_A1_m_plugin(k), ci, stats] = ttest(FullData_plugin_all(:,i), HalfData_plugin_all(:,i));
    [h, P_joint_A1_m_qe(k), ci, stats] = ttest(FullData_qe_all(:,i), HalfData_qe_all(:,i));
    [h, P_joint_A1_m_shuffSub(k), ci, stats] = ttest(FullData_shuffSub_all(:,i), HalfData_shuffSub_all(:,i));
    [h, P_joint_A1_m_qeshuffSub(k), ci, stats] = ttest(FullData_qeshuffSub_all(:,i), HalfData_qeshuffSub_all(:,i));
end 


meanpluginA1_monkey = [mean(FullData_plugin_all(:,1)), mean(HalfData_plugin_all(:,1)), ...
                mean(FullData_plugin_all(:,3)), mean(HalfData_plugin_all(:,3)),...
                mean(FullData_plugin_all(:,2)), mean(HalfData_plugin_all(:,2))];

sempluginA1_monkey = [ std(FullData_plugin_all(:,1))/sqrt(size(FullData_plugin_all,1)), std(HalfData_plugin_all(:,1))/sqrt(size(HalfData_plugin_all,1)), ...
                std(FullData_plugin_all(:,3))/sqrt(size(FullData_plugin_all,1)), std(HalfData_plugin_all(:,3))/sqrt(size(HalfData_plugin_all,1)),...
                std(FullData_plugin_all(:,2))/sqrt(size(FullData_plugin_all,1)), std(HalfData_plugin_all(:,2))/sqrt(size(HalfData_plugin_all,1))];


meanqeA1_monkey = [mean(FullData_qe_all(:,1)), mean(HalfData_qe_all(:,1)), ...
            mean(FullData_qe_all(:,3)), mean(HalfData_qe_all(:,3)),...
            mean(FullData_qe_all(:,2)), mean(HalfData_qe_all(:,2))];

semqeA1_monkey = [std(FullData_qe_all(:,1))/sqrt(size(FullData_qe_all,1)), std(HalfData_qe_all(:,1))/sqrt(size(HalfData_qe_all,1)), ...
           std(FullData_qe_all(:,3))/sqrt(size(FullData_qe_all,1)), std(HalfData_qe_all(:,3))/sqrt(size(HalfData_qe_all,1)),...
           std(FullData_qe_all(:,2))/sqrt(size(FullData_qe_all,1)), std(HalfData_qe_all(:,2))/sqrt(size(HalfData_qe_all,1))];


meanshuffSubA1_monkey = [mean(FullData_shuffSub_all(:,1)), mean(HalfData_shuffSub_all(:,1)), ...
                  mean(FullData_shuffSub_all(:,3)), mean(HalfData_shuffSub_all(:,3)),...
                  mean(FullData_shuffSub_all(:,2)), mean(HalfData_shuffSub_all(:,2))];

semshuffSubA1_monkey = [std(FullData_shuffSub_all(:,1))/sqrt(size(FullData_shuffSub_all,1)), std(HalfData_shuffSub_all(:,1))/sqrt(size(HalfData_shuffSub_all,1)), ...
                 std(FullData_shuffSub_all(:,3))/sqrt(size(FullData_shuffSub_all,1)), std(HalfData_shuffSub_all(:,3))/sqrt(size(HalfData_shuffSub_all,1)),...
                 std(FullData_shuffSub_all(:,2))/sqrt(size(FullData_shuffSub_all,1)), std(HalfData_shuffSub_all(:,2))/sqrt(size(HalfData_shuffSub_all,1))];

meanqeshuffSubA1_monkey = [mean(FullData_qeshuffSub_all(:,1)), mean(HalfData_qeshuffSub_all(:,1)), ...
                    mean(FullData_qeshuffSub_all(:,3)), mean(HalfData_qeshuffSub_all(:,3)),...
                    mean(FullData_qeshuffSub_all(:,2)), mean(HalfData_qeshuffSub_all(:,2))];

semqeshuffSubA1_monkey = [std(FullData_qeshuffSub_all(:,1))/sqrt(size(FullData_qeshuffSub_all,1)), std(HalfData_qeshuffSub_all(:,1))/sqrt(size(HalfData_qeshuffSub_all,1)), ...
                   std(FullData_qeshuffSub_all(:,3))/sqrt(size(FullData_qeshuffSub_all,1)), std(HalfData_qeshuffSub_all(:,3))/sqrt(size(HalfData_qeshuffSub_all,1)),...
                   std(FullData_qeshuffSub_all(:,2))/sqrt(size(FullData_qeshuffSub_all,1)), std(HalfData_qeshuffSub_all(:,2))/sqrt(size(HalfData_qeshuffSub_all,1))];

%% Plot
figure('Position', [100, 100, 600, 600]); 
tiledlayout(3,4);

maxV = 0.35;
subplotIdx = 1;
plotPIDterms(meanpluginA1, sempluginA1, subplotIdx, 'plugin', 0, maxV, false)
subplotIdx = subplotIdx + 1;
plotPIDterms(meanqeA1, semqeA1, subplotIdx, 'qe', 0, maxV, false)
subplotIdx = subplotIdx + 1;
plotPIDterms(meanshuffSubA1, semshuffSubA1, subplotIdx, 'shuffSub' ,0 ,maxV, false)
subplotIdx = subplotIdx + 1;
plotPIDterms(meanqeshuffSubA1, semqeshuffSubA1, subplotIdx, 'qe shuffSub', 0,maxV, false)


maxV = 0.15;
subplotIdx = subplotIdx + 1;
plotPIDterms(meanpluginCA1, sempluginCA1, subplotIdx, 'plugin', 0,maxV, false)
subplotIdx = subplotIdx + 1;
plotPIDterms(meanqeCA1, semqeCA1, subplotIdx, 'qe', 0,maxV, false)
subplotIdx = subplotIdx + 1;
plotPIDterms(meanshuffSubCA1, semshuffSubCA1, subplotIdx, 'shuffSub',0,maxV, false)
subplotIdx = subplotIdx + 1;
plotPIDterms(meanqeshuffSubCA1, semqeshuffSubCA1, subplotIdx, 'qe shuffSub',0,maxV, false)

maxV = 0.1;
minV = 0;
subplotIdx = subplotIdx + 1;
plotPIDterms(meanpluginA1_monkey, sempluginA1_monkey, subplotIdx, 'plugin', minV, maxV, false)
subplotIdx = subplotIdx + 1;
plotPIDterms(meanqeA1_monkey, semqeA1_monkey, subplotIdx, 'qe', minV,maxV, false)
subplotIdx = subplotIdx + 1;
plotPIDterms(meanshuffSubA1_monkey, semshuffSubA1_monkey, subplotIdx, 'shuffSub',minV,maxV, false)
subplotIdx = subplotIdx + 1;
plotPIDterms(meanqeshuffSubA1_monkey, semqeshuffSubA1_monkey, subplotIdx, 'qe shuffSub',minV,maxV, false)


%% === Save ===
outDir = fullfile(fileparts(pwd), 'Figures_mat');  

if ~exist(outDir, 'dir')
    mkdir(outDir);
end

exportgraphics(gcf, fullfile(outDir, 'Figure_S15.svg'), 'ContentType','vector');

%% Transposed Markdown Table for Real Data: Only Red & Syn, Δ and Normalized Δ

methods = {'plugin', 'qe', 'shuffSub', 'qeshuffSub'};
areas = {'A1', 'CA1', 'A1_monkey'};
components = {'Red', 'Syn'};
component_indices = [2, 3];  % Indices for Red and Syn in original mean arrays

absTable = {};  % Absolute difference
normTable = {}; % Normalized difference

% Build header
header = {'Component'};
for a = 1:length(areas)
    for m = 1:length(methods)
        header{end+1} = sprintf('%s-%s', areas{a}, methods{m});
    end
end

% Loop over each component (Red, Syn)
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
            
            compRowAbs{end+1} = sprintf('%.4f', delta);
            compRowNorm{end+1} = sprintf('%.2f%%', 100 * normDelta);
        end
    end
    
    absTable = [absTable; compRowAbs];
    normTable = [normTable; compRowNorm];
end

% Output as markdown
fprintf('### Table RXXa: Absolute Δ(Half - Full) in bits\n');
fprintf('| %s |\n', strjoin(header, ' | '));
fprintf('|%s|\n', repmat(' --- |', 1, numel(header)));
for i = 1:size(absTable,1)
    fprintf('| %s |\n', strjoin(absTable(i,:), ' | '));
end

fprintf('\n### Table RXXb: Δ Normalized to Full (%%)\n');
fprintf('| %s |\n', strjoin(header, ' | '));
fprintf('|%s|\n', repmat(' --- |', 1, numel(header)));
for i = 1:size(normTable,1)
    fprintf('| %s |\n', strjoin(normTable(i,:), ' | '));
end


%%


function plotPIDterms(meanJoint, semJoint, subplotIdx, plotTitle, minLim, maxLim, isPercent)
    barLabels = {'Joint','Red', 'Syn'};
    colors = [0.9290, 0.6940, 0.1250;  % Originalfarben
          0.9645, 0.8470, 0.5625;  % Hellerer Gelbton
          0, 0.4470, 0.7410;       % Original Blau
          0.5, 0.6720, 0.8705;     % Helleres Blau
          0.4660, 0.6740, 0.1880;  % Original Grün
          0.7330, 0.8370, 0.5940]; % Helleres Grün
    nexttile(subplotIdx);
    b = bar(meanJoint, 'FaceColor', 'flat', 'EdgeColor', 'k');
    for i = 1:length(meanJoint)
        b.CData(i, :) = colors(i, :);
    end
    hold on;
    errorbar(1:length(meanJoint), meanJoint, semJoint, 'k', 'LineStyle', 'none', 'LineWidth', 1.5);
    xticks(1.5:2:length(meanJoint)); 
    xticklabels(barLabels);   
    if isPercent
        ylabel('% of Joint Info');
    else
        ylabel('Information [bits]');
    end
    title(plotTitle, 'Interpreter', 'none');
    ylim([minLim, maxLim])
    grid on;
    box on;
    hold off;
end

function add_pvalue(p_values, x1, x2, text_heights, line_heights)
for i = 1:numel(p_values)
    stars = get_stars(p_values(i));
    line([x1(i), x2(i)], [line_heights(i),line_heights(i)], 'Color', 'black', 'LineWidth', 1);
    text((x1(i) + x2(i)) / 2, line_heights(i) + text_heights, stars, 'HorizontalAlignment', 'center', 'Color', 'black', 'FontSize', 12);
end
end

function stars = get_stars(p_value)
if p_value < 0.001
    stars = '***';
elseif p_value < 0.01
    stars = '**';
elseif p_value < 0.05
    stars = '*';
else
    stars = 'n.s';
end
end

