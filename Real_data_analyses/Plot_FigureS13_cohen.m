clc, clear, close all;

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

    cohenD_A1_plugin(k)     = computeCohensD_independent(HalfData_plugin_all(:,i), FullData_plugin_all(:,i));
    cohenD_A1_qe(k)         = computeCohensD_independent(HalfData_qe_all(:,i), FullData_qe_all(:,i));
    cohenD_A1_shuffSub(k)   = computeCohensD_independent(HalfData_shuffSub_all(:,i), FullData_shuffSub_all(:,i));
    cohenD_A1_qeshuffSub(k) = computeCohensD_independent(HalfData_qeshuffSub_all(:,i), FullData_qeshuffSub_all(:,i));

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
    cohenD_CA1_plugin(k)     = computeCohensD_independent(HalfData_plugin_all(:,i),FullData_plugin_all(:,i));
    cohenD_CA1_qe(k)         = computeCohensD_independent(HalfData_qe_all(:,i),FullData_qe_all(:,i));
    cohenD_CA1_shuffSub(k)   = computeCohensD_independent(HalfData_shuffSub_all(:,i),FullData_shuffSub_all(:,i));
    cohenD_CA1_qeshuffSub(k) = computeCohensD_independent(HalfData_qeshuffSub_all(:,i),FullData_qeshuffSub_all(:,i));

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
    cohenD_A1_monkey_plugin(k)     = computeCohensD_independent(HalfData_plugin_all(:,i),FullData_plugin_all(:,i));
    cohenD_A1_monkey_qe(k)         = computeCohensD_independent(HalfData_qe_all(:,i),FullData_qe_all(:,i));
    cohenD_A1_monkey_shuffSub(k)   = computeCohensD_independent(HalfData_shuffSub_all(:,i),FullData_shuffSub_all(:,i));
    cohenD_A1_monkey_qeshuffSub(k) = computeCohensD_independent(HalfData_qeshuffSub_all(:,i),FullData_qeshuffSub_all(:,i));
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

%% Transposed Markdown Tables: Half - Full, Normalized Δ, and Cohen's d
methods = {'plugin', 'qe', 'shuffSub', 'qeshuffSub'};
areas = {'A1', 'CA1', 'A1_monkey'};
components = {'Red', 'Syn'};
component_indices = [2, 3];  % Red = 2, Syn = 3

% Loop over each component separately (Red, Syn)
for c = 1:length(components)
    compName = components{c};
    i = component_indices(c);  % 2 for Red, 3 for Syn

    % Initialize tables
    absTable = {};
    normTable = {};
    cohenTable = {};

    % Create headers: methods
    methodHeader = [{'Dataset'}, methods];

    % Loop over areas (rows)
    for a = 1:length(areas)
        area = areas{a};

        % Initialize row entries
        rowAbs = {area};
        rowNorm = {area};
        rowCohen = {area};

        % Loop over methods (columns)
        for m = 1:length(methods)
            method = methods{m};

            % Mean and SEM
            meanVar = eval(sprintf('mean%s%s', method, area));
            semVar  = eval(sprintf('sem%s%s', method, area));
            fullVal = meanVar(2*i - 1);
            halfVal = meanVar(2*i);
            fullSEM = semVar(2*i - 1);
            halfSEM = semVar(2*i);

            % Δ and SEM
            delta = halfVal - fullVal;
            semDelta = sqrt(fullSEM^2 + halfSEM^2);
            rowAbs{end+1} = sprintf('%.4f ± %.4f', delta, semDelta);

            % Normalized to qeshuffSub Full
            refMean = eval(sprintf('meanqeshuffSub%s', area));
            refVal = refMean(2*i - 1);
            normDelta = 100 * (delta / refVal);
            normSEM   = 100 * (semDelta / refVal);
            rowNorm{end+1} = sprintf('%.2f%% ± %.2f%%', normDelta, normSEM);

            % Cohen's d
            cohenVar = eval(sprintf('cohenD_%s_%s', area, method));
            d = cohenVar(i);
            rowCohen{end+1} = sprintf('%.2f', d);
        end

        % Append to tables
        absTable = [absTable; rowAbs];
        normTable = [normTable; rowNorm];
        cohenTable = [cohenTable; rowCohen];
    end

    % ----------- Print Markdown Tables -------------
    fprintf('\n### Table: Absolute Δ(Half - Full) ± SEM [bits] - %s\n', compName);
    fprintf('| %s |\n', strjoin(methodHeader, ' | '));
    fprintf('|%s|\n', repmat(' --- |', 1, numel(methodHeader)));
    for i = 1:size(absTable,1)
        fprintf('| %s |\n', strjoin(absTable(i,:), ' | '));
    end

    fprintf('\n### Table: Normalized Δ(Half - Full) to qeshuffSub Full (%% ± SEM) - %s\n', compName);
    fprintf('| %s |\n', strjoin(methodHeader, ' | '));
    fprintf('|%s|\n', repmat(' --- |', 1, numel(methodHeader)));
    for i = 1:size(normTable,1)
        fprintf('| %s |\n', strjoin(normTable(i,:), ' | '));
    end

    fprintf('\n### Table: Cohen''s d for Δ(Half - Full) - %s\n', compName);
    fprintf('| %s |\n', strjoin(methodHeader, ' | '));
    fprintf('|%s|\n', repmat(' --- |', 1, numel(methodHeader)));
    for i = 1:size(cohenTable,1)
        fprintf('| %s |\n', strjoin(cohenTable(i,:), ' | '));
    end
end





%%


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

function d = computeCohensD(x1, x2)
    diff = x1 - x2;
    d = mean(diff) / std(diff);
end

function d = computeCohensD_independent(x1, x2)
    n1 = length(x1);
    n2 = length(x2);
    s1 = std(x1);
    s2 = std(x2);

    pooled_std = sqrt( ((n1 - 1)*s1^2 + (n2 - 1)*s2^2) / (n1 + n2 - 2) );
    d = (mean(x1) - mean(x2)) / pooled_std;
end
