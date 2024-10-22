function plot_bar_plot_three_layers_two_group(layer1valueG1,layer1valueG2, layer2valueG1, layer2valueG2, layer3valueG1, layer3valueG2)
layer1valueG1(isnan(layer1valueG1)) = [];
layer1valueG2(isnan(layer1valueG2)) = [];
layer2valueG1(isnan(layer2valueG1)) = [];
layer2valueG2(isnan(layer2valueG2)) = [];
layer3valueG1(isnan(layer3valueG1)) = [];
layer3valueG2(isnan(layer3valueG2)) = [];

data = [[mean(layer1valueG1), mean(layer2valueG1), mean(layer3valueG1)];
    [mean(layer1valueG2), mean(layer2valueG2), mean(layer3valueG2)]];
err = [[std(layer1valueG1)/sqrt(length(layer1valueG1)), std(layer2valueG1)/sqrt(length(layer2valueG1)), std(layer3valueG1)/sqrt(length(layer3valueG1))];
    [std(layer1valueG2)/sqrt(length(layer1valueG2)), std(layer2valueG2)/sqrt(length(layer2valueG2)), std(layer3valueG2)/sqrt(length(layer3valueG2))]];

figure;  % Creates a new figure
b = bar(data,'BarWidth', 0.8);
b(1).FaceColor = [0.8500, 0.3250, 0.0980];  % Color for Variable 1
b(2).FaceColor = [0, 0.4470, 0.7410];  % Color for Variable 2
b(3).FaceColor = [0.9290, 0.6940, 0.1250];  % Color for Variable 3% Plots grouped bar chart
ylabel('radius (nm)');  % Y-axis label for the values of the variables
title('neck redius in thick and thin dendrites');  % Chart title
numGroups = size(data, 1);
numBars = size(data, 2);
% Add error bars
hold on;
for i = 1:numBars
    % X positions for error bars
    x = b(i).XEndPoints;
    % Add error bars to each group
    errorbar(x, data(:, i), err(:, i), '.k', 'LineWidth', 1.5);
end

groupNames = {'smaller dendrite radius', 'larger dendrite radius'}; % Custom names for each group of bars
xticks(1:length(groupNames)); % Set the x-axis ticks to match the number of groups
xticklabels(groupNames);



bar_plot_mean = cell(3,1);
bar_plot_mean{1} = data(:,1);
bar_plot_mean{2} = data(:,2);
bar_plot_mean{3} = data(:,3);
bar_plot_test = cell(3,2);
bar_plot_test{1,1} = layer1valueG1;
bar_plot_test{1,2} = layer1valueG2;
bar_plot_test{2,1} = layer2valueG1;
bar_plot_test{2,2} = layer2valueG2;
bar_plot_test{3,1} = layer3valueG1;
bar_plot_test{3,2} = layer3valueG2;
hold on;
[position_cell, pvalueArray] = plot_significance_score_cell_cell(bar_plot_mean, bar_plot_test, b);
sigstar(position_cell,pvalueArray)
legend('L2/3', 'L4', 'L5')



end


