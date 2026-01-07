

%%
clc; clear;

x = -5:1:100;
color_poiss =  [135, 112, 54]./255;
color_gauss = [224, 200, 98] ./ 255;
fontSize = 18; linewidth = 1.5;
axisLineWidth = 1.5; tickFontSize = 18;

x_short = -5:1:14;  
lambdas_panel1 = [0.5, 2, 5];
y_poiss_panel1 = arrayfun(@(l) poisspdf(x_short, l), lambdas_panel1, 'UniformOutput', false);
y_gauss_panel1 = arrayfun(@(l) normpdf(x_short, l, sqrt(l)), lambdas_panel1, 'UniformOutput', false);

lambdas_2 = [2, 7]+3;
y_2 = cellfun(@(l) poisspdf(x, l), num2cell(lambdas_2), 'UniformOutput', false);
z_poiss_2 = mean(cell2mat(y_2'), 1);
spikes = x.*z_poiss_2;
meani = sum(spikes);
vari = sum(((x-meani).^2).*z_poiss_2);
z_gauss_2 = normpdf(x, sum(spikes), sqrt(vari));

lambdas_4 = [2, 7, 12, 17]+3;
y_4 = cellfun(@(l) poisspdf(x, l), num2cell(lambdas_4), 'UniformOutput', false);
z_poiss_4 = mean(cell2mat(y_4'), 1);
spikes = x.*z_poiss_4;
meani = sum(spikes);
vari = sum(((x-meani).^2).*z_poiss_4);
z_gauss_4 = normpdf(x, sum(spikes), sqrt(vari));

lambdas_12 = linspace(2, 57, 12)+3;
y_12 = cellfun(@(l) poisspdf(x, l), num2cell(lambdas_12), 'UniformOutput', false);
z_poiss_12 = mean(cell2mat(y_12'), 1);
spikes = x.*z_poiss_12;
meani = sum(spikes);
vari = sum(((x-meani).^2).*z_poiss_12);
z_gauss_12 = normpdf(x, sum(spikes), sqrt(vari));


figure('Units', 'centimeters', 'Position', [1, 1, 40, 10]);

subplot(1,4,1); hold on
for i = 1:3
    plot(x_short, y_poiss_panel1{i}, 'Color', color_poiss, 'LineWidth', 2);
    plot(x_short, y_gauss_panel1{i}, 'Color', color_gauss, 'LineWidth', 2);
end
title('P(r|s)', 'FontSize', fontSize)
xlabel('Spike count'); ylabel('Probability')
ylim([0 0.7]); xlim([-5 20])
set(gca, 'FontSize', tickFontSize, 'LineWidth', axisLineWidth)

subplot(1,4,2); hold on
plot(x, z_poiss_2, 'Color', color_poiss, 'LineWidth', linewidth);
plot(x, z_gauss_2, 'Color', color_gauss, 'LineWidth', linewidth);
title('P(r)', 'FontSize', fontSize)
xlabel('Spike count'); ylabel('Probability')
ylim([0 0.2]); xlim([-5 70])
set(gca, 'FontSize', tickFontSize, 'LineWidth', axisLineWidth)
 
subplot(1,4,3); hold on
plot(x, z_poiss_4, 'Color', color_poiss, 'LineWidth', linewidth);
plot(x, z_gauss_4, 'Color', color_gauss, 'LineWidth', linewidth);
title('P(r)', 'FontSize', fontSize)
xlabel('Spike count'); ylabel('Probability')
ylim([0 0.2]); xlim([-5 70])
set(gca, 'FontSize', tickFontSize, 'LineWidth', axisLineWidth)

subplot(1,4,4); hold on
plot(x, z_poiss_12, 'Color', color_poiss, 'LineWidth', linewidth);
plot(x, z_gauss_12, 'Color', color_gauss, 'LineWidth', linewidth);
title('P(r)', 'FontSize', fontSize)
xlabel('Spike count'); ylabel('Probability')
ylim([0 0.2]); xlim([-5 70])
set(gca, 'FontSize', tickFontSize, 'LineWidth', axisLineWidth)


outDir = fullfile(fileparts(pwd), 'Figures_mat');   % ein Verzeichnis über pwd
set(gcf, 'Renderer', 'painters');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end
exportgraphics(gcf, fullfile(outDir, 'Figure_S12_AB.svg'), 'ContentType','vector');
