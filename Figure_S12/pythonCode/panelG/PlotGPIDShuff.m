% Load the CSV file
load("FigR1E/discrete_shuffSub_results.mat");

dataCA1 = readtable('resultsCA1_shuffCorr.csv');
dataA1 = readtable('resultsA1_shuffCorr.csv');
dataPPC = readtable('resultsPPC_shuffCorr.csv');

% Extract mean and SEM for Ijoint
mean_IjointCA1 = dataCA1.mean_Ijoint;
sem_IjointCA1 = dataCA1.sem_Ijoint;
mean_IjointCA1_gpid = dataCA1.mean_Ijoint_gpid;
sem_IjointCA1_gpid  = dataCA1.sem_Ijoint_gpid;
mean_IjointCA1_discrete = mean(dataShuff.IJoint_CA1);
sem_IjointCA1_discrete = std(dataShuff.IJoint_CA1) / sqrt(numel(dataShuff.IJoint_CA1));

mean_IjointA1 = dataA1.mean_Ijoint;
sem_IjointA1  = dataA1.sem_Ijoint;
mean_IjointA1_gpid = dataA1.mean_Ijoint_gpid;
sem_IjointA1_gpid  = dataA1.sem_Ijoint_gpid;
mean_IjointA1_discrete = mean(dataShuff.IJoint_A1);
sem_IjointA1_discrete = std(dataShuff.IJoint_A1) / sqrt(numel(dataShuff.IJoint_A1));

mean_IjointPPC = dataPPC.mean_Ijoint;
sem_IjointPPC = dataPPC.sem_Ijoint;
mean_IjointPPC_gpid = dataPPC.mean_Ijoint_gpid;
sem_IjointPPC_gpid  = dataPPC.sem_Ijoint_gpid;
mean_IjointPPC_discrete = mean(dataShuff.IJoint_PPC);
sem_IjointPPC_discrete = std(dataShuff.IJoint_PPC) / sqrt(numel(dataShuff.IJoint_PPC));

idx = 1;

figure;
hold on;
bar(idx, mean_IjointCA1_discrete, 'FaceColor', 'b');
errorbar(idx, mean_IjointCA1_discrete, sem_IjointCA1_discrete, 'k', 'LineStyle', 'none', 'LineWidth', 2);idx = idx+1;
bar(idx, mean_IjointCA1, 'FaceColor', 'b');
errorbar(idx, mean_IjointCA1, sem_IjointCA1, 'k', 'LineStyle', 'none', 'LineWidth', 2);idx = idx+1;
bar(idx, mean_IjointCA1_gpid, 'FaceColor', 'b');
errorbar(idx, mean_IjointCA1_gpid, sem_IjointCA1_gpid, 'k', 'LineStyle', 'none', 'LineWidth', 2);idx = idx+1;

bar(idx, mean_IjointA1_discrete, 'FaceColor', 'b');
errorbar(idx, mean_IjointA1_discrete, sem_IjointA1_discrete, 'k', 'LineStyle', 'none', 'LineWidth', 2);idx = idx+1;
bar(idx, mean_IjointA1, 'FaceColor', 'b');
errorbar(idx, mean_IjointA1, sem_IjointA1, 'k', 'LineStyle', 'none', 'LineWidth', 2);idx = idx+1;
bar(idx, mean_IjointA1_gpid, 'FaceColor', 'b');
errorbar(idx, mean_IjointA1_gpid, sem_IjointA1_gpid, 'k', 'LineStyle', 'none', 'LineWidth', 2);idx = idx+1;

bar(idx, mean_IjointPPC_discrete, 'FaceColor', 'b');
errorbar(idx, mean_IjointPPC_discrete, sem_IjointPPC_discrete, 'k', 'LineStyle', 'none', 'LineWidth', 2);idx = idx+1;
bar(idx, mean_IjointPPC, 'FaceColor', 'b');
errorbar(idx, mean_IjointPPC, sem_IjointPPC, 'k', 'LineStyle', 'none', 'LineWidth', 2);idx = idx+1;
bar(idx, mean_IjointPPC_gpid, 'FaceColor', 'b');
errorbar(idx, mean_IjointPPC_gpid, sem_IjointPPC_gpid, 'k', 'LineStyle', 'none', 'LineWidth', 2);
 hold off;