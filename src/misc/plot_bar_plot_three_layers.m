function plot_bar_plot_three_layers(layer1value, layer2value, layer3value)

layer1value(isnan(layer1value)) = [];
layer2value(isnan(layer2value)) = [];
layer3value(isnan(layer3value)) = [];
data = [mean(layer1value), mean(layer2value), mean(layer3value)];
err = [[std(layer1value)/sqrt(length(layer1value)), std(layer2value)/sqrt(length(layer2value)), std(layer3value)/sqrt(length(layer3value))]];
figure;  % Creates a new figure
b = bar(data,'BarWidth', 0.8);
b.FaceColor = 'flat';
set(gca, 'XTickLabel',{'L2/3','L4','L5'})
b.CData(1,:) = [0.8500, 0.3250, 0.0980];
b.CData(2,:) = [0, 0.4470, 0.7410];
b.CData(3,:) = [0.9290, 0.6940, 0.1250];
% b(1).FaceColor = [0.8500, 0.3250, 0.0980];  % Color for Variable 1
% b(2).FaceColor = [0, 0.4470, 0.7410];  % Color for Variable 2
% b(3).FaceColor = [0.9290, 0.6940, 0.1250];  % Color for Variable 3% Plots grouped bar chart
ylabel('radius (nm)');  % Y-axis label for the values of the variables
title('neck redius across layers');  % Chart title
numGroups = size(data, 1);
numBars = size(data, 2);
% Add error bars
hold on;
for i = 1:numBars
    % X positions for error bars
    x = b.XEndPoints(i);
    % Add error bars to each group
    errorbar(x, data(:, i), err(:, i), '.k', 'LineWidth', 1.5);
end


bar_plot_mean = cell(3,1);
bar_plot_mean{1} = data(:,1);
bar_plot_mean{2} = data(:,2);
bar_plot_mean{3} = data(:,3);
bar_plot_test = cell(3,2);
bar_plot_test{1,1} = layer1value;
bar_plot_test{2,1} = layer2value;
bar_plot_test{3,1} = layer3value;
hold on;
[position_cell, pvalueArray] = plot_significance_score_cell_cell(bar_plot_mean, bar_plot_test, b);
sigstar(position_cell,pvalueArray)



end