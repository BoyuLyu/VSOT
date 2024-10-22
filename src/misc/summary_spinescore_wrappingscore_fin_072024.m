clear
% singleSynHeadVolume = [];
% singleSynHeadMeanRadius = [];
% singleSynNeckLength = [];
% singleSynNeckSection = [];
% singleSynNeckMeanRadius = [];
% singleDendriteDensity = [];
% doubleSynHeadVolume = [];
% doubleSynMeanHeadRadius = [];
% doubleSynNeckLength = [];
% doubleSynNeckSection = [];
% doubleSynNeckMeanRadius = [];
% doubleDendriteDensity = [];
% singleSynapticCleftSize = [];
% sinsperimeterRatio = [];
% sinsperimeterWeightedWrappingArea = [];
% sinspostSynapseTouchingArea = [];
% sinspostSynapseTouchingRatio = [];
% sinspreSynapseTouchingArea = [];
% sinspreSynapseTouchingRatio = [];
% singleSynHeadNeckTouchingArea = [];
% singleSynHeadNeckTouchingRatio = [];
% doubleSynapticCleftSize = [];
% dousperimeterRatio = [];
% dousperimeterWeightedWrappingArea = [];
% douspostSynapseTouchingArea = [];
% douspostSynapseTouchingRatio = [];
% douspreSynapseTouchingArea = [];
% douspreSynapseTouchingRatio = [];
% doubleSynHeadNeckTouchingArea = [];
% doubleSynHeadNeckTouchingRatio = [];
% doubleSynPosition = [];
% spineCoordinate = [];
% singleSynapse_astroID_neuronID_folder = [];
%[2,4,5,6,7,8,9,10,12,13,14,15,16,17,18,19]
%% go over all the folders to collect the quantification results
rootFolder = '/work/boyu/EM_astrocyte/astro_11_33_16_16_40/';
neuronListFolder = '/work/boyu/EM_astrocyte/astro_11_33_64_64_80/';
addpath('../resources/sigstar-master')
addpath('../resources/Violinplot-Matlab/')


nfList = [17:22, 31:33];
layer5_summary = summary_wrapping_score_subfunc(nfList, rootFolder, neuronListFolder);
layer5_dendrite_level1_summary = summary_dendrite_score_level1(nfList, rootFolder, neuronListFolder);
nfList = [11:16, 29:30];
layer4_summary = summary_wrapping_score_subfunc(nfList, rootFolder, neuronListFolder);
layer4_dendrite_level1_summary = summary_dendrite_score_level1(nfList, rootFolder, neuronListFolder);

nfList = [23:28];
layer2_summary = summary_wrapping_score_subfunc(nfList, rootFolder, neuronListFolder);
layer2_dendrite_level1_summary = summary_dendrite_score_level1(nfList, rootFolder, neuronListFolder);


%% assign all the values to each vector
spine_dendrite_type_label5 = layer5_dendrite_level1_summary.spine_dendrite_type_label; % [1: filopodia, 2: mushroom, 3: long thin, 4: thin, 5: stubby, dendrite_label] Ns x 2
dendrite_length5  = layer5_dendrite_level1_summary.dendrite_length; % Nd x 1
dendrite_radius5 = layer5_dendrite_level1_summary.dendrite_radius; % Nd x 1 
dendrite_each_spine_type_number5 = layer5_dendrite_level1_summary.dendrite_each_spine_type_number; % specifically add the 6th type as the branched spine, simply to check between different layers Nd x 6
neuron_total_type_number5 = layer5_dendrite_level1_summary.neuron_total_type_number; % Nn x 6

spine_dendrite_type_label4 = layer4_dendrite_level1_summary.spine_dendrite_type_label; % [1: filopodia, 2: mushroom, 3: long thin, 4: thin, 5: stubby, dendrite_label] Ns x 2
dendrite_length4  = layer4_dendrite_level1_summary.dendrite_length; % Nd x 1
dendrite_radius4 = layer4_dendrite_level1_summary.dendrite_radius; % Nd x 1 
dendrite_each_spine_type_number4 = layer4_dendrite_level1_summary.dendrite_each_spine_type_number; % specifically add the 6th type as the branched spine, simply to check between different layers Nd x 6
neuron_total_type_number4 = layer4_dendrite_level1_summary.neuron_total_type_number; % Nn x 6

spine_dendrite_type_label2 = layer2_dendrite_level1_summary.spine_dendrite_type_label; % [1: filopodia, 2: mushroom, 3: long thin, 4: thin, 5: stubby, dendrite_label] Ns x 2
dendrite_length2  = layer2_dendrite_level1_summary.dendrite_length; % Nd x 1
dendrite_radius2 = layer2_dendrite_level1_summary.dendrite_radius; % Nd x 1 
dendrite_each_spine_type_number2 = layer2_dendrite_level1_summary.dendrite_each_spine_type_number; % specifically add the 6th type as the branched spine, simply to check between different layers Nd x 6
neuron_total_type_number2 = layer2_dendrite_level1_summary.neuron_total_type_number; % Nn x 6


singleSynHeadVolume5 = layer5_summary.singleSynHeadVolume;
singleSynMeanHeadRadius5 = layer5_summary.singleSynMeanHeadRadius;
singleSynNeckLength5 = layer5_summary.singleSynNeckLength;
singleSynNeckSection5 = layer5_summary.singleSynNeckSection;
singleSynNeckMeanRadius5 = layer5_summary.singleSynNeckMeanRadius;
doubleSynHeadVolume5 = layer5_summary.doubleSynHeadVolume;
doubleSynMeanHeadRadius5 = layer5_summary.doubleSynMeanHeadRadius;
doubleSynNeckLength5 = layer5_summary.doubleSynNeckLength;
doubleSynNeckSection5 = layer5_summary.doubleSynNeckSection;
doubleSynNeckMeanRadius5 = layer5_summary.doubleSynNeckMeanRadius;
singleSynapticCleftSize5 = layer5_summary.singleSynapticCleftSize;
% singleDendriteDensity5 = layer5_summary.singleDendriteDensity;
sinsperimeterRatio5 = layer5_summary.sinsperimeterRatio;
sinsperimeterWeightedWrappingArea5 = layer5_summary.sinsperimeterWeightedWrappingArea;
sinspostSynapseTouchingArea5 = layer5_summary.sinspostSynapseTouchingArea;
sinspostSynapseTouchingRatio5 = layer5_summary.sinspostSynapseTouchingRatio;
sinspreSynapseTouchingArea5 = layer5_summary.sinspreSynapseTouchingArea;
sinspreSynapseTouchingRatio5 = layer5_summary.sinspreSynapseTouchingRatio;
singleSynHeadNeckTouchingArea5 = layer5_summary.singleSynHeadNeckTouchingArea;
singleSynHeadNeckTouchingRatio5 = layer5_summary.singleSynHeadNeckTouchingRatio;
% doubleDendriteDensity5 = layer5_summary.doubleDendriteDensity;
doubleSynapticCleftSize5 = layer5_summary.doubleSynapticCleftSize;
dousperimeterRatio5 = layer5_summary.dousperimeterRatio;
dousperimeterWeightedWrappingArea5 = layer5_summary.dousperimeterWeightedWrappingArea;
douspostSynapseTouchingArea5  = layer5_summary.douspostSynapseTouchingArea;
douspostSynapseTouchingRatio5 = layer5_summary.douspostSynapseTouchingRatio;
douspreSynapseTouchingArea5  = layer5_summary.douspreSynapseTouchingArea;
douspreSynapseTouchingRatio5 = layer5_summary.douspreSynapseTouchingRatio;
doubleSynHeadNeckTouchingArea5 = layer5_summary.doubleSynHeadNeckTouchingArea;
doubleSynHeadNeckTouchingRatio5 = layer5_summary.doubleSynHeadNeckTouchingRatio;
% doubleSynPosition5 = layer5_summary.doubleSynPosition;
% spineCoordinate5  =layer5_summary.spineCoordinate;





singleSynHeadVolume2 = layer2_summary.singleSynHeadVolume;
singleSynMeanHeadRadius2 = layer2_summary.singleSynMeanHeadRadius;
singleSynNeckLength2 = layer2_summary.singleSynNeckLength;
singleSynNeckSection2 = layer2_summary.singleSynNeckSection;
singleSynNeckMeanRadius2 = layer2_summary.singleSynNeckMeanRadius;
doubleSynHeadVolume2 = layer2_summary.doubleSynHeadVolume;
doubleSynMeanHeadRadius2 = layer2_summary.doubleSynMeanHeadRadius;
doubleSynNeckLength2 = layer2_summary.doubleSynNeckLength;
doubleSynNeckSection2 = layer2_summary.doubleSynNeckSection;
doubleSynNeckMeanRadius2 = layer2_summary.doubleSynNeckMeanRadius;
singleSynapticCleftSize2 = layer2_summary.singleSynapticCleftSize;
sinsperimeterRatio2 = layer2_summary.sinsperimeterRatio;
sinsperimeterWeightedWrappingArea2 = layer2_summary.sinsperimeterWeightedWrappingArea;
sinspostSynapseTouchingArea2 = layer2_summary.sinspostSynapseTouchingArea;
sinspostSynapseTouchingRatio2 = layer2_summary.sinspostSynapseTouchingRatio;
sinspreSynapseTouchingArea2 = layer2_summary.sinspreSynapseTouchingArea;
sinspreSynapseTouchingRatio2 = layer2_summary.sinspreSynapseTouchingRatio;
singleSynHeadNeckTouchingArea2 = layer2_summary.singleSynHeadNeckTouchingArea;
singleSynHeadNeckTouchingRatio2 = layer2_summary.singleSynHeadNeckTouchingRatio;
doubleSynapticCleftSize2 = layer2_summary.doubleSynapticCleftSize;
dousperimeterRatio2 = layer2_summary.dousperimeterRatio;
dousperimeterWeightedWrappingArea2 = layer2_summary.dousperimeterWeightedWrappingArea;
douspostSynapseTouchingArea2  = layer2_summary.douspostSynapseTouchingArea;
douspostSynapseTouchingRatio2 = layer2_summary.douspostSynapseTouchingRatio;
douspreSynapseTouchingArea2 = layer2_summary.douspreSynapseTouchingArea;
douspreSynapseTouchingRatio2 = layer2_summary.douspreSynapseTouchingRatio;
doubleSynHeadNeckTouchingArea2 = layer2_summary.doubleSynHeadNeckTouchingArea;
doubleSynHeadNeckTouchingRatio2 = layer2_summary.doubleSynHeadNeckTouchingRatio;
% doubleSynPosition2 = layer2_summary.doubleSynPosition;
% spineCoordinate2  =layer2_summary.spineCoordinate;
% singleDendriteDensity2 = layer2_summary.singleDendriteDensity;
% doubleDendriteDensity2 = layer2_summary.doubleDendriteDensity;

singleSynHeadVolume4 = layer4_summary.singleSynHeadVolume;
singleSynMeanHeadRadius4 = layer4_summary.singleSynMeanHeadRadius;
singleSynNeckLength4 = layer4_summary.singleSynNeckLength;
singleSynNeckSection4 = layer4_summary.singleSynNeckSection;
singleSynNeckMeanRadius4 = layer4_summary.singleSynNeckMeanRadius;
doubleSynHeadVolume4 = layer4_summary.doubleSynHeadVolume;
doubleSynMeanHeadRadius4 = layer4_summary.doubleSynMeanHeadRadius;
doubleSynNeckLength4 = layer4_summary.doubleSynNeckLength;
doubleSynNeckSection4 = layer4_summary.doubleSynNeckSection;
doubleSynNeckMeanRadius4 = layer4_summary.doubleSynNeckMeanRadius;
singleSynapticCleftSize4 = layer4_summary.singleSynapticCleftSize;
sinsperimeterRatio4 = layer4_summary.sinsperimeterRatio;
sinsperimeterWeightedWrappingArea4 = layer4_summary.sinsperimeterWeightedWrappingArea;
sinspostSynapseTouchingArea4 = layer4_summary.sinspostSynapseTouchingArea;
sinspostSynapseTouchingRatio4 = layer4_summary.sinspostSynapseTouchingRatio;
sinspreSynapseTouchingArea4 = layer4_summary.sinspreSynapseTouchingArea;
sinspreSynapseTouchingRatio4 = layer4_summary.sinspreSynapseTouchingRatio;
singleSynHeadNeckTouchingArea4 = layer4_summary.singleSynHeadNeckTouchingArea;
singleSynHeadNeckTouchingRatio4 = layer4_summary.singleSynHeadNeckTouchingRatio;
doubleSynapticCleftSize4 = layer4_summary.doubleSynapticCleftSize;
dousperimeterRatio4 = layer4_summary.dousperimeterRatio;
dousperimeterWeightedWrappingArea4 = layer4_summary.dousperimeterWeightedWrappingArea;
douspostSynapseTouchingArea4  = layer4_summary.douspostSynapseTouchingArea;
douspostSynapseTouchingRatio4 = layer4_summary.douspostSynapseTouchingRatio;
douspreSynapseTouchingArea4  = layer4_summary.douspreSynapseTouchingArea;
douspreSynapseTouchingRatio4 = layer4_summary.douspreSynapseTouchingRatio;
doubleSynHeadNeckTouchingArea4 = layer4_summary.doubleSynHeadNeckTouchingArea;
doubleSynHeadNeckTouchingRatio4 = layer4_summary.doubleSynHeadNeckTouchingRatio;

%% Neck length
singleSynNeckLength2(singleSynNeckLength2(:) == 0) = [];
singleSynNeckLength4(singleSynNeckLength4(:) == 0) = [];
singleSynNeckLength5(singleSynNeckLength5(:) == 0) = [];
singleSynNeckLength2(singleSynNeckLength2 > quantile(singleSynNeckLength2, 0.99))  = [];
singleSynNeckLength4(singleSynNeckLength4 > quantile(singleSynNeckLength4, 0.99))  = [];
singleSynNeckLength5(singleSynNeckLength5 > quantile(singleSynNeckLength5, 0.99))  = [];


x = [singleSynNeckLength2(:);singleSynNeckLength4(:);singleSynNeckLength5(:)];
% x = x*sqrt(16^2 + 40^2)/1000;
y = [zeros(length(singleSynNeckLength2(:)),1);ones(length(singleSynNeckLength4(:)),1);ones(length(singleSynNeckLength5(:)),1).*2];

figure; violinplot(x,y,'ShowMean',true,'ShowMedian',true); title(['dendrite spine neck length'])
set(gca,'xtick',[1,2,3],'xticklabel',{'L2/3', 'L4', 'L5'}); ylabel('nm')
% id2 = randperm(length(singleSynNeckLength2), 3000);
% id4 = randperm(length(singleSynNeckLength4), 3000);
% id5 = randperm(length(singleSynNeckLength5), 3000);
[h,p1,ci,stats] = ttest2(singleSynNeckLength2(:),singleSynNeckLength4(:));
[h,p2,ci,stats] = ttest2(singleSynNeckLength4(:),singleSynNeckLength5(:));
[h,p3,ci,stats] = ttest2(singleSynNeckLength2(:),singleSynNeckLength5(:));
sigstar({[1,2],[2,3], [1,3]},[p1,p2,p3])


% bar plot

singleSynNeckLength2(singleSynNeckLength2 == 0 |isinf(singleSynNeckLength2)|isnan(singleSynNeckLength2)) = [];
singleSynNeckLength2(singleSynNeckLength2 > quantile(singleSynNeckLength2, 0.99)) = [];
singleSynNeckLength4(singleSynNeckLength4 == 0 |isinf(singleSynNeckLength4)|isnan(singleSynNeckLength4)) = [];
singleSynNeckLength4(singleSynNeckLength4 > quantile(singleSynNeckLength4, 0.99)) = [];
singleSynNeckLength5(singleSynNeckLength5 == 0 |isinf(singleSynNeckLength5)|isnan(singleSynNeckLength5)) = [];
singleSynNeckLength5(singleSynNeckLength5 > quantile(singleSynNeckLength5, 0.99)) = [];

data = [mean(singleSynNeckLength2), mean(singleSynNeckLength4), mean(singleSynNeckLength5)];
err = [[std(singleSynNeckLength2)/sqrt(length(singleSynNeckLength2)), std(singleSynNeckLength4)/sqrt(length(singleSynNeckLength4))...
, std(singleSynNeckLength5)/sqrt(length(singleSynNeckLength5))]];
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
title('dendrite spine neck length');  % Chart title
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
bar_plot_test{1,1} = singleSynNeckLength2;
bar_plot_test{2,1} = singleSynNeckLength4;
bar_plot_test{3,1} = singleSynNeckLength5;
hold on;
[position_cell, pvalueArray] = plot_significance_score_cell_cell(bar_plot_mean, bar_plot_test, b);
sigstar(position_cell,pvalueArray)




% bar plot of neck length comparing between dendrite with large radius and small radius

singleSynNeckLength2_small = singleSynNeckLength2(spine_small_radius_id_2);
singleSynNeckLength2_small(singleSynNeckLength2_small == 0) = [];
singleSynNeckLength2_small(singleSynNeckLength2_small > quantile(singleSynNeckLength2_small, 0.99)) = [];
singleSynNeckLength2_small(isnan(singleSynNeckLength2_small)) = [];
singleSynNeckLength2_small(isinf(singleSynNeckLength2_small)) = [];

singleSynNeckLength4_small = singleSynNeckLength4(spine_small_radius_id_4);
singleSynNeckLength4_small(singleSynNeckLength4_small == 0) = [];
singleSynNeckLength4_small(singleSynNeckLength4_small > quantile(singleSynNeckLength4_small, 0.99)) = [];
singleSynNeckLength4_small(isnan(singleSynNeckLength4_small)) = [];
singleSynNeckLength4_small(isinf(singleSynNeckLength4_small)) = [];

singleSynNeckLength5_small = singleSynNeckLength5(spine_small_radius_id_5);
singleSynNeckLength5_small(singleSynNeckLength5_small == 0) = [];
singleSynNeckLength5_small(singleSynNeckLength5_small > quantile(singleSynNeckLength5_small, 0.99)) = [];
singleSynNeckLength5_small(isnan(singleSynNeckLength5_small)) = [];
singleSynNeckLength5_small(isinf(singleSynNeckLength5_small)) = [];

singleSynNeckLength2_large = singleSynNeckLength2(spine_large_radius_id_2);
singleSynNeckLength2_large(singleSynNeckLength2_large == 0) = [];
singleSynNeckLength2_large(singleSynNeckLength2_large > quantile(singleSynNeckLength2_large, 0.99)) = [];
singleSynNeckLength2_large(isnan(singleSynNeckLength2_large)) = [];
singleSynNeckLength2_large(isinf(singleSynNeckLength2_large)) = [];

singleSynNeckLength4_large = singleSynNeckLength4(spine_large_radius_id_4);
singleSynNeckLength4_large(singleSynNeckLength4_large == 0) = [];
singleSynNeckLength4_large(singleSynNeckLength4_large > quantile(singleSynNeckLength4_large, 0.99)) = [];
singleSynNeckLength4_large(isnan(singleSynNeckLength4_large)) = [];
singleSynNeckLength4_large(isinf(singleSynNeckLength4_large)) = [];

singleSynNeckLength5_large = singleSynNeckLength5(spine_large_radius_id_5);
singleSynNeckLength5_large(singleSynNeckLength5_large == 0) = [];
singleSynNeckLength5_large(singleSynNeckLength5_large > quantile(singleSynNeckLength5_large, 0.99)) = [];
singleSynNeckLength5_large(isnan(singleSynNeckLength5_large)) = [];
singleSynNeckLength5_large(isinf(singleSynNeckLength5_large)) = [];

data = [[mean(singleSynNeckLength2_small), mean(singleSynNeckLength4_small), mean(singleSynNeckLength5_small)];
    [mean(singleSynNeckLength2_large), mean(singleSynNeckLength4_large), mean(singleSynNeckLength5_large)]];
err = [[std(singleSynNeckLength2_small)/sqrt(length(singleSynNeckLength2_small)), std(singleSynNeckLength4_small)/sqrt(length(singleSynNeckLength4_small)), std(singleSynNeckLength5_small)/sqrt(length(singleSynNeckLength5_small))];
    [std(singleSynNeckLength2_large)/sqrt(length(singleSynNeckLength2_large)), std(singleSynNeckLength4_large)/sqrt(length(singleSynNeckLength4_large)), std(singleSynNeckLength5_large)/sqrt(length(singleSynNeckLength5_large))]];

figure;  % Creates a new figure
b = bar(data,'BarWidth', 0.8);
b(1).FaceColor = [0.8500, 0.3250, 0.0980];  % Color for Variable 1
b(2).FaceColor = [0, 0.4470, 0.7410];  % Color for Variable 2
b(3).FaceColor = [0.9290, 0.6940, 0.1250];  % Color for Variable 3% Plots grouped bar chart
ylabel('radius (nm)');  % Y-axis label for the values of the variables
title('neck length in thick and thin dendrites');  % Chart title
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
bar_plot_test{1,1} = singleSynNeckLength2_small;
bar_plot_test{1,2} = singleSynNeckLength2_large;
bar_plot_test{2,1} = singleSynNeckLength4_small;
bar_plot_test{2,2} = singleSynNeckLength4_large;
bar_plot_test{3,1} = singleSynNeckLength5_small;
bar_plot_test{3,2} = singleSynNeckLength5_large;
hold on;
[position_cell, pvalueArray] = plot_significance_score_cell_cell(bar_plot_mean, bar_plot_test, b);
sigstar(position_cell,pvalueArray)
legend('L2/3', 'L4', 'L5')







% NECK LENGTH COMPARE INSIDE EACH LAYERS   
singleSynNeckLength2_small = singleSynNeckLength2(spine_small_radius_id_2);
singleSynNeckLength2_small(singleSynNeckLength2_small == 0) = [];
singleSynNeckLength2_small(singleSynNeckLength2_small > quantile(singleSynNeckLength2_small, 0.99)) = [];
singleSynNeckLength2_small(isnan(singleSynNeckLength2_small)) = [];
singleSynNeckLength2_small(isinf(singleSynNeckLength2_small)) = [];

singleSynNeckLength4_small = singleSynNeckLength4(spine_small_radius_id_4);
singleSynNeckLength4_small(singleSynNeckLength4_small == 0) = [];
singleSynNeckLength4_small(singleSynNeckLength4_small > quantile(singleSynNeckLength4_small, 0.99)) = [];
singleSynNeckLength4_small(isnan(singleSynNeckLength4_small)) = [];
singleSynNeckLength4_small(isinf(singleSynNeckLength4_small)) = [];

singleSynNeckLength5_small = singleSynNeckLength5(spine_small_radius_id_5);
singleSynNeckLength5_small(singleSynNeckLength5_small == 0) = [];
singleSynNeckLength5_small(singleSynNeckLength5_small > quantile(singleSynNeckLength5_small, 0.99)) = [];
singleSynNeckLength5_small(isnan(singleSynNeckLength5_small)) = [];
singleSynNeckLength5_small(isinf(singleSynNeckLength5_small)) = [];

singleSynNeckLength2_large = singleSynNeckLength2(spine_large_radius_id_2);
singleSynNeckLength2_large(singleSynNeckLength2_large == 0) = [];
singleSynNeckLength2_large(singleSynNeckLength2_large > quantile(singleSynNeckLength2_large, 0.99)) = [];
singleSynNeckLength2_large(isnan(singleSynNeckLength2_large)) = [];
singleSynNeckLength2_large(isinf(singleSynNeckLength2_large)) = [];

singleSynNeckLength4_large = singleSynNeckLength4(spine_large_radius_id_4);
singleSynNeckLength4_large(singleSynNeckLength4_large == 0) = [];
singleSynNeckLength4_large(singleSynNeckLength4_large > quantile(singleSynNeckLength4_large, 0.99)) = [];
singleSynNeckLength4_large(isnan(singleSynNeckLength4_large)) = [];
singleSynNeckLength4_large(isinf(singleSynNeckLength4_large)) = [];

singleSynNeckLength5_large = singleSynNeckLength5(spine_large_radius_id_5);
singleSynNeckLength5_large(singleSynNeckLength5_large == 0) = [];
singleSynNeckLength5_large(singleSynNeckLength5_large > quantile(singleSynNeckLength5_large, 0.99)) = [];
singleSynNeckLength5_large(isnan(singleSynNeckLength5_large)) = [];
singleSynNeckLength5_large(isinf(singleSynNeckLength5_large)) = [];

plot_bar_plot_three_layers_two_group_withinlayers(singleSynNeckLength2_small,singleSynNeckLength2_large,...
    singleSynNeckLength4_small, singleSynNeckLength4_large, singleSynNeckLength5_small, singleSynNeckLength5_large)
ylabel('radius (nm)');  % Y-axis label for the values of the variables
title('layer-specific neck length in thick and thin dendrites');  % Chart title






%% head volume
singleSynHeadVolume2(singleSynHeadVolume2 == 0) = [];
singleSynHeadVolume4(singleSynHeadVolume4 == 0) = [];
singleSynHeadVolume5(singleSynHeadVolume5 == 0) = [];
singleSynHeadVolume2 = singleSynHeadVolume2.*16.*16.*40/10^9;
singleSynHeadVolume4 = singleSynHeadVolume4.*16.*16.*40/10^9;
singleSynHeadVolume5 = singleSynHeadVolume5.*16.*16.*40/10^9;
% x = [singleSynHeadVolume2(:);singleSynHeadVolume4(:);singleSynHeadVolume5(:)];
% x = x.*16.*16.*40/10^9;
% y = [zeros(length(singleSynHeadVolume2(:)),1);ones(length(singleSynHeadVolume4(:)),1);ones(length(singleSynHeadVolume5(:)),1).*2];
data = [mean(singleSynHeadVolume2), mean(singleSynHeadVolume4), mean(singleSynHeadVolume5)];
err = [[std(singleSynHeadVolume2)/sqrt(length(singleSynHeadVolume2)), std(singleSynHeadVolume4)/sqrt(length(singleSynHeadVolume4))...
, std(singleSynHeadVolume5)/sqrt(length(singleSynHeadVolume5))]];
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
ylabel('volume(\mu m^3)');
title(['volume of dendrite spine head'])
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
bar_plot_test{1,1} = singleSynHeadVolume2;
bar_plot_test{2,1} = singleSynHeadVolume4;
bar_plot_test{3,1} = singleSynHeadVolume5;
hold on;
[position_cell, pvalueArray] = plot_significance_score_cell_cell(bar_plot_mean, bar_plot_test, b);
sigstar(position_cell,pvalueArray)





figure; violinplot(x,y,'ShowMean',true,'ShowMedian',true); title(['volume of dendrite spine head'])
set(gca,'xtick',[1,2,3],'xticklabel',{'L2', 'L4', 'L5'})
ylabel('volume(\mu m^3)')
[h,p1,ci,stats] = ttest2(singleSynHeadVolume2(:),singleSynHeadVolume4(:));
[h,p2,ci,stats] = ttest2(singleSynHeadVolume4(:),singleSynHeadVolume5(:));
[h,p3,ci,stats] = ttest2(singleSynHeadVolume2(:),singleSynHeadVolume5(:));
sigstar({[1,2],[2,3], [1,3]},[p1,p2,p3])

% bar plot
std5 = std(singleSynHeadVolume5,1);
types5_mean = mean(singleSynHeadVolume5,1);

std4 = std(singleSynHeadVolume4, 1);
types4_mean = mean(singleSynHeadVolume4,1);

std2 = std(singleSynHeadVolume2, 1);
types2_mean = mean(singleSynHeadVolume2,1);

data = [types2_mean, types4_mean,types5_mean];
err = [std2, std4, std5];
figure;  % Creates a new figure
b = bar(data,'BarWidth', 0.8);
% b(1).FaceColor = [0.8500, 0.3250, 0.0980];  % Color for Variable 1
% b(2).FaceColor = [0, 0.4470, 0.7410];  % Color for Variable 2
% b(3).FaceColor = [0.9290, 0.6940, 0.1250];  % Color for Variable 3% Plots grouped bar chart
% xlabel({'Contact', 'contact at head', 'contact at neck'});  % X-axis label for categories/groups
ylabel('volume(\mu m^3)');  % Y-axis label for the values of the variables
title('Distribution of spine types in each neuron');  % Chart title


[f1, xi1] = ksdensity(singleSynHeadVolume2);  % KDE for dataset 1
[f2, xi2] = ksdensity(singleSynHeadVolume4);  % KDE for dataset 2
[f3, xi3] = ksdensity(singleSynHeadVolume5);  % KDE for dataset 3
figure;
plot(xi1, f1, 'b-', 'LineWidth', 2); hold on;
plot(xi2, f2, 'r-', 'LineWidth', 2); hold on;
plot(xi3, f3, 'g-', 'LineWidth', 2)
legend('L2/3', 'L4', 'L5')

title('distribution of the dendrite spine head volume')
ylabel('probability')
xlabel('volume (\mu m^3)')



% since we can see clearly the two modes in layer 5, which might be caused
% by the apparent thicker apical dendrites. Next, split the whole dendrites
% into thick dendrites and thin dendrites.
reset_values_summary_statistics;
[small_radius_id_2, large_radius_id_2] = gmm_fit_radius(dendrite_radius2);%the label of dendrites
[small_radius_id_4, large_radius_id_4] = gmm_fit_radius(dendrite_radius4);
[small_radius_id_5, large_radius_id_5] = gmm_fit_radius(dendrite_radius5);

aa = [dendrite_radius2;dendrite_radius4;dendrite_radius5];
[small_radius_id_aa, large_radius_id_aa] = gmm_fit_radius(aa);
split_point = max(aa(small_radius_id_aa));
small_radius_id_2 = find(dendrite_radius2 < split_point & dendrite_radius2 > 0);
small_radius_id_4 = find(dendrite_radius4 < split_point & dendrite_radius4 > 0);
small_radius_id_5 = find(dendrite_radius5 < split_point & dendrite_radius5 > 0);
large_radius_id_2 = find(dendrite_radius2 > split_point & dendrite_radius2 > 0);
large_radius_id_4 = find(dendrite_radius4 > split_point & dendrite_radius4 > 0);
large_radius_id_5 = find(dendrite_radius5 > split_point & dendrite_radius5 > 0);

spine_dendrite_idx2 = label2idx(spine_dendrite_type_label2(:,2));
spine_dendrite_idx2 = spine_dendrite_idx2(:);
spine_dendrite_idx4 = label2idx(spine_dendrite_type_label4(:,2));
spine_dendrite_idx4 = spine_dendrite_idx4(:);
spine_dendrite_idx5 = label2idx(spine_dendrite_type_label5(:,2));
spine_dendrite_idx5 = spine_dendrite_idx5(:);
spine_small_radius_id_2 = cell2mat(spine_dendrite_idx2(small_radius_id_2));
spine_small_radius_id_4 = cell2mat(spine_dendrite_idx4(small_radius_id_4));
spine_small_radius_id_5 = cell2mat(spine_dendrite_idx5(small_radius_id_5));

spine_large_radius_id_2 = cell2mat(spine_dendrite_idx2(large_radius_id_2));
spine_large_radius_id_4 = cell2mat(spine_dendrite_idx4(large_radius_id_4));
spine_large_radius_id_5 = cell2mat(spine_dendrite_idx5(large_radius_id_5));

% [small_radius_id_aa, large_radius_id_aa] = gmm_fit_radius(aa);
% split_point = max(aa(small_radius_id_aa));
% small_radius_id_2 = find(dendrite_radius2 < split_point & dendrite_radius2 > 0);
% small_radius_id_4 = find(dendrite_radius4 < split_point & dendrite_radius4 > 0);
% small_radius_id_5 = find(dendrite_radius5 < split_point & dendrite_radius5 > 0);
plot_violin_L23L4L5(singleSynHeadVolume2(spine_large_radius_id_2), singleSynHeadVolume4(spine_large_radius_id_4),...
    singleSynHeadVolume5(spine_large_radius_id_5), 'head_volume')
title('volume of dendrite spine head (larger dendritic radius)')
plot_violin_L23L4L5(singleSynHeadVolume2(spine_small_radius_id_2), singleSynHeadVolume4(spine_small_radius_id_4),...
    singleSynHeadVolume5(spine_small_radius_id_5), 'head_volume')
title('volume of dendrite spine head (smaller dendritic radius)')


singleSynHeadVolume2 = singleSynHeadVolume2.*16.*16.*40/10^9;
singleSynHeadVolume4 = singleSynHeadVolume4.*16.*16.*40/10^9;
singleSynHeadVolume5 = singleSynHeadVolume5.*16.*16.*40/10^9;
singleSynHeadVolume2_small = singleSynHeadVolume2(spine_small_radius_id_2);
singleSynHeadVolume2_small(singleSynHeadVolume2_small == 0) = [];
singleSynHeadVolume4_small = singleSynHeadVolume4(spine_small_radius_id_4);
singleSynHeadVolume4_small(singleSynHeadVolume4_small == 0) = [];
singleSynHeadVolume5_small = singleSynHeadVolume5(spine_small_radius_id_5);
singleSynHeadVolume5_small(singleSynHeadVolume5_small == 0) = [];
singleSynHeadVolume2_large = singleSynHeadVolume2(spine_large_radius_id_2);
singleSynHeadVolume2_large(singleSynHeadVolume2_large == 0) = [];
singleSynHeadVolume4_large = singleSynHeadVolume4(spine_large_radius_id_4);
singleSynHeadVolume4_large(singleSynHeadVolume4_large == 0) = [];
singleSynHeadVolume5_large = singleSynHeadVolume5(spine_large_radius_id_5);
singleSynHeadVolume5_large(singleSynHeadVolume5_large == 0) = [];

data = [[mean(singleSynHeadVolume2_small), mean(singleSynHeadVolume4_small), mean(singleSynHeadVolume5_small)];
    [mean(singleSynHeadVolume2_large), mean(singleSynHeadVolume4_large), mean(singleSynHeadVolume5_large)]];
err = [[std(singleSynHeadVolume2_small)/sqrt(length(singleSynHeadVolume2_small)), std(singleSynHeadVolume4_small)/sqrt(length(singleSynHeadVolume4_small)), std(singleSynHeadVolume5_small)/sqrt(length(singleSynHeadVolume5_small))];
    [std(singleSynHeadVolume2_large)/sqrt(length(singleSynHeadVolume2_large)), std(singleSynHeadVolume4_large)/sqrt(length(singleSynHeadVolume4_large)), std(singleSynHeadVolume5_large)/sqrt(length(singleSynHeadVolume5_large))]];

% data = [[mean(singleSynHeadVolume2(spine_large_radius_id_2)),mean(singleSynHeadVolume4(spine_large_radius_id_4)), mean(singleSynHeadVolume4(spine_large_radius_id_4)) ];
%     [mean(singleSynHeadVolume2(spine_small_radius_id_4)),mean(singleSynHeadVolume4(spine_small_radius_id_4)), mean(singleSynHeadVolume4(spine_small_radius_id_4)) ]];
% err = [[std(singleSynHeadVolume2(spine_large_radius_id_2)),std(singleSynHeadVolume4(spine_large_radius_id_4)), std(singleSynHeadVolume4(spine_large_radius_id_4)) ];
%     [std(singleSynHeadVolume2(spine_small_radius_id_4)),std(singleSynHeadVolume4(spine_small_radius_id_4)), std(singleSynHeadVolume4(spine_small_radius_id_4)) ]];
figure;  % Creates a new figure
b = bar(data,'BarWidth', 0.8);
b(1).FaceColor = [0.8500, 0.3250, 0.0980];  % Color for Variable 1
b(2).FaceColor = [0, 0.4470, 0.7410];  % Color for Variable 2
b(3).FaceColor = [0.9290, 0.6940, 0.1250];  % Color for Variable 3% Plots grouped bar chart
ylabel('volume (\mu m^3)');  % Y-axis label for the values of the variables
title('head volume in thick and thin dendrites');  % Chart title
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
bar_plot_test{1,1} = singleSynHeadVolume2_small;
bar_plot_test{1,2} = singleSynHeadVolume2_large;
bar_plot_test{2,1} = singleSynHeadVolume4_small;
bar_plot_test{2,2} = singleSynHeadVolume4_large;
bar_plot_test{3,1} = singleSynHeadVolume5_small;
bar_plot_test{3,2} = singleSynHeadVolume5_large;
hold on;
[position_cell, pvalueArray] = plot_significance_score_cell_cell(bar_plot_mean, bar_plot_test, b);
sigstar(position_cell,pvalueArray)
legend('L2/3', 'L4', 'L5')






reset_values_summary_statistics;
singleSynHeadVolume2 = singleSynHeadVolume2.*16.*16.*40/10^9;
singleSynHeadVolume4 = singleSynHeadVolume4.*16.*16.*40/10^9;
singleSynHeadVolume5 = singleSynHeadVolume5.*16.*16.*40/10^9;
singleSynHeadVolume2_small = singleSynHeadVolume2(spine_small_radius_id_2);
singleSynHeadVolume2_small(singleSynHeadVolume2_small == 0) = [];
singleSynHeadVolume4_small = singleSynHeadVolume4(spine_small_radius_id_4);
singleSynHeadVolume4_small(singleSynHeadVolume4_small == 0) = [];
singleSynHeadVolume5_small = singleSynHeadVolume5(spine_small_radius_id_5);
singleSynHeadVolume5_small(singleSynHeadVolume5_small == 0) = [];
singleSynHeadVolume2_large = singleSynHeadVolume2(spine_large_radius_id_2);
singleSynHeadVolume2_large(singleSynHeadVolume2_large == 0) = [];
singleSynHeadVolume4_large = singleSynHeadVolume4(spine_large_radius_id_4);
singleSynHeadVolume4_large(singleSynHeadVolume4_large == 0) = [];
singleSynHeadVolume5_large = singleSynHeadVolume5(spine_large_radius_id_5);
singleSynHeadVolume5_large(singleSynHeadVolume5_large == 0) = [];
plot_bar_plot_three_layers_two_group_withinlayers(singleSynHeadVolume2_small,singleSynHeadVolume2_large,...
    singleSynHeadVolume4_small, singleSynHeadVolume4_large, singleSynHeadVolume5_small, singleSynHeadVolume5_large)
ylabel('volume (\mu m^3)');  % Y-axis label for the values of the variables
title('layer-specific head volume comparison');  % Chart title





%% whether dendrite radius and dendrite spine head volume are correlated
id2 = find(dendrite_radius2 > 0 & singleSynHeadVolume2 > 0);
id4 = find(dendrite_radius4 > 0 & singleSynHeadVolume4 > 0);
id5 = find(dendrite_radius5 > 0 & singleSynHeadVolume5 > 0);

plot_regression_L2L4L5(dendrite_radius2(id2), dendrite_radius4(id4), dendrite_radius5(id5), ...
    singleSynHeadVolume2(id2), singleSynHeadVolume4(id4), singleSynHeadVolume5(id5))






figure; histogram(sorted_dendrite_radius2(small_group)); hold on; histogram(sorted_dendrite_radius2(large_group))




groupNames = {'L2/3', 'L4', 'L5'}; % Custom names for each group of bars
xticks(1:length(groupNames)); % Set the x-axis ticks to match the number of groups
xticklabels(groupNames);

hold on;
numGroups = size(data, 1);
numBars = size(data, 2);

% Add error bars
for i = 1:numBars
    % X positions for error bars
    x = b(i).XEndPoints;
    % Add error bars to each group
    errorbar(x, data(:, i), err(:, i), '.k', 'LineWidth', 1.5);
end
hold on;


types5_reorg = [types5(:,1),types5(:,3), types5(:,4), types5(:,5), types5(:,2), types5(:,6)];
types4_reorg = [types4(:,1),types4(:,3), types4(:,4), types4(:,5), types4(:,2), types4(:,6)];
types2_reorg = [types2(:,1),types2(:,3), types2(:,4), types2(:,5), types2(:,2), types2(:,6)];
bar_plot_mean = cell(3,1);
bar_plot_mean{1} = types2_mean;
bar_plot_mean{2} = types4_mean;
bar_plot_mean{3} = types5_mean;
bar_plot_test = cell(3,1);
bar_plot_test{1} = types2_reorg;
bar_plot_test{2} = types4_reorg;
bar_plot_test{3} = types5_reorg;
hold on;
[position_cell, pvalueArray] = plot_significance_score(bar_plot_mean, bar_plot_test, b);
sigstar(position_cell,pvalueArray)



%% head radius
singleSynMeanHeadRadius2 = singleSynMeanHeadRadius2(singleSynMeanHeadRadius2~=0);
singleSynMeanHeadRadius4 = singleSynMeanHeadRadius4(singleSynMeanHeadRadius4~=0);
singleSynMeanHeadRadius5 = singleSynMeanHeadRadius5(singleSynMeanHeadRadius5~=0);
singleSynMeanHeadRadius2(singleSynMeanHeadRadius2> quantile(singleSynMeanHeadRadius2, 0.99)) = [];
singleSynMeanHeadRadius4(singleSynMeanHeadRadius4> quantile(singleSynMeanHeadRadius4, 0.99)) = [];
singleSynMeanHeadRadius5(singleSynMeanHeadRadius5> quantile(singleSynMeanHeadRadius5, 0.99)) = [];
x = [singleSynMeanHeadRadius2(:);singleSynMeanHeadRadius4(:);singleSynMeanHeadRadius5(:)];

y = [zeros(length(singleSynMeanHeadRadius2(:)),1);ones(length(singleSynMeanHeadRadius4(:)),1);ones(length(singleSynMeanHeadRadius5(:)),1).*2];

figure; violinplot(x,y,'ShowMean',true,'ShowMedian',true); title(['dendrite spine head radius'])
set(gca,'xtick',[1,2,3],'xticklabel',{'L2/3', 'L4', 'L5'}); ylabel('radius(nm)')
% figure; violinplot(x,y,'ShowMean',true,'ShowMedian',false); title(['dendrite spine neck section'])
% set(gca,'xtick',[1,2,3],'xticklabel',{'L2/3', 'L4', 'L5'})
% id2 = randperm(length(singleSynMeanHeadRadius2), 3000);
% id4 = randperm(length(singleSynMeanHeadRadius4), 3000);
% id5 = randperm(length(singleSynMeanHeadRadius5), 3000);
[h,p1,ci,stats] = ttest2(singleSynMeanHeadRadius2,singleSynMeanHeadRadius4);
[h,p2,ci,stats] = ttest2(singleSynMeanHeadRadius4,singleSynMeanHeadRadius5);
[h,p3,ci,stats] = ttest2(singleSynMeanHeadRadius2,singleSynMeanHeadRadius5);
sigstar({[1,2],[2,3], [1,3]},[p1,p2,p3])
% histogram
figure; f1 = histogram(singleSynMeanHeadRadius2,100, "Normalization","probability"); 
hold on; histogram(singleSynMeanHeadRadius4,100,"Normalization","probability", "BinEdges",f1.BinEdges); 
hold on; histogram(singleSynMeanHeadRadius5,100, "Normalization","probability", "BinEdges",f1.BinEdges); 

legend('L2/3', 'L4', 'L5')
title('histogram of the radius for the max section in head')
xlabel('nm')







data = [mean(singleSynMeanHeadRadius2), mean(singleSynMeanHeadRadius4), mean(singleSynMeanHeadRadius5)];
err = [[std(singleSynMeanHeadRadius2)/sqrt(length(singleSynMeanHeadRadius2)), std(singleSynMeanHeadRadius4)/sqrt(length(singleSynMeanHeadRadius4)),...
 std(singleSynMeanHeadRadius5)/sqrt(length(singleSynMeanHeadRadius5))]];
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
title('dendrite spine head radius');  % Chart title
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
bar_plot_test{1,1} = singleSynMeanHeadRadius2;
bar_plot_test{2,1} = singleSynMeanHeadRadius4;
bar_plot_test{3,1} = singleSynMeanHeadRadius5;
hold on;
[position_cell, pvalueArray] = plot_significance_score_cell_cell(bar_plot_mean, bar_plot_test, b);
sigstar(position_cell,pvalueArray)
% split the head 

singleSynMeanHeadRadius2_small = singleSynMeanHeadRadius2(spine_small_radius_id_2);
singleSynMeanHeadRadius2_small(singleSynMeanHeadRadius2_small == 0) = [];
singleSynMeanHeadRadius2_small(singleSynMeanHeadRadius2_small > quantile(singleSynMeanHeadRadius2_small, 0.99)) = [];
singleSynMeanHeadRadius4_small = singleSynMeanHeadRadius4(spine_small_radius_id_4);
singleSynMeanHeadRadius4_small(singleSynMeanHeadRadius4_small == 0) = [];
singleSynMeanHeadRadius4_small(singleSynMeanHeadRadius4_small > quantile(singleSynMeanHeadRadius4_small, 0.99)) = [];
singleSynMeanHeadRadius5_small = singleSynMeanHeadRadius5(spine_small_radius_id_5);
singleSynMeanHeadRadius5_small(singleSynMeanHeadRadius5_small == 0) = [];
singleSynMeanHeadRadius5_small(singleSynMeanHeadRadius5_small > quantile(singleSynMeanHeadRadius5_small, 0.99)) = [];
singleSynMeanHeadRadius2_large = singleSynMeanHeadRadius2(spine_large_radius_id_2);
singleSynMeanHeadRadius2_large(singleSynMeanHeadRadius2_large == 0) = [];
singleSynMeanHeadRadius2_large(singleSynMeanHeadRadius2_large > quantile(singleSynMeanHeadRadius2_large, 0.99)) = [];
singleSynMeanHeadRadius4_large = singleSynMeanHeadRadius4(spine_large_radius_id_4);
singleSynMeanHeadRadius4_large(singleSynMeanHeadRadius4_large == 0) = [];
singleSynMeanHeadRadius4_large(singleSynMeanHeadRadius4_large > quantile(singleSynMeanHeadRadius4_large, 0.99)) = [];
singleSynMeanHeadRadius5_large = singleSynMeanHeadRadius5(spine_large_radius_id_5);
singleSynMeanHeadRadius5_large(singleSynMeanHeadRadius5_large == 0) = [];
singleSynMeanHeadRadius5_large(singleSynMeanHeadRadius5_large > quantile(singleSynMeanHeadRadius5_large, 0.99)) = [];

data = [[mean(singleSynMeanHeadRadius2_small), mean(singleSynMeanHeadRadius4_small), mean(singleSynMeanHeadRadius5_small)];
    [mean(singleSynMeanHeadRadius2_large), mean(singleSynMeanHeadRadius4_large), mean(singleSynMeanHeadRadius5_large)]];
err = [[std(singleSynMeanHeadRadius2_small)/sqrt(length(singleSynMeanHeadRadius2_small)), std(singleSynMeanHeadRadius4_small)/sqrt(length(singleSynMeanHeadRadius4_small)), std(singleSynMeanHeadRadius5_small)/sqrt(length(singleSynMeanHeadRadius5_small))];
    [std(singleSynMeanHeadRadius2_large)/sqrt(length(singleSynMeanHeadRadius2_large)), std(singleSynMeanHeadRadius4_large)/sqrt(length(singleSynMeanHeadRadius4_large)), std(singleSynMeanHeadRadius5_large)/sqrt(length(singleSynMeanHeadRadius5_large))]];


% data = [[mean(singleSynHeadVolume2_small), mean(singleSynHeadVolume4_small), mean(singleSynHeadVolume5_small)];
%     [mean(singleSynHeadVolume2_large), mean(singleSynHeadVolume4_large), mean(singleSynHeadVolume5_large)]];
% err = [[std(singleSynHeadVolume2_small)/sqrt(length(singleSynHeadVolume2_small)), std(singleSynHeadVolume4_small)/sqrt(length(singleSynHeadVolume4_small)), std(singleSynHeadVolume5_small)/sqrt(length(singleSynHeadVolume5_small))];
%     [std(singleSynHeadVolume2_large)/sqrt(length(singleSynHeadVolume2_large)), std(singleSynHeadVolume4_large)/sqrt(length(singleSynHeadVolume4_large)), std(singleSynHeadVolume5_large)/sqrt(length(singleSynHeadVolume5_large))]];

% data = [[mean(singleSynHeadVolume2(spine_large_radius_id_2)),mean(singleSynHeadVolume4(spine_large_radius_id_4)), mean(singleSynHeadVolume4(spine_large_radius_id_4)) ];
%     [mean(singleSynHeadVolume2(spine_small_radius_id_4)),mean(singleSynHeadVolume4(spine_small_radius_id_4)), mean(singleSynHeadVolume4(spine_small_radius_id_4)) ]];
% err = [[std(singleSynHeadVolume2(spine_large_radius_id_2)),std(singleSynHeadVolume4(spine_large_radius_id_4)), std(singleSynHeadVolume4(spine_large_radius_id_4)) ];
%     [std(singleSynHeadVolume2(spine_small_radius_id_4)),std(singleSynHeadVolume4(spine_small_radius_id_4)), std(singleSynHeadVolume4(spine_small_radius_id_4)) ]];
figure;  % Creates a new figure
b = bar(data,'BarWidth', 0.8);
b(1).FaceColor = [0.8500, 0.3250, 0.0980];  % Color for Variable 1
b(2).FaceColor = [0, 0.4470, 0.7410];  % Color for Variable 2
b(3).FaceColor = [0.9290, 0.6940, 0.1250];  % Color for Variable 3% Plots grouped bar chart
ylabel('radius (nm)');  % Y-axis label for the values of the variables
title('head redius in thick and thin dendrites');  % Chart title
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
bar_plot_test{1,1} = singleSynMeanHeadRadius2_small;
bar_plot_test{1,2} = singleSynMeanHeadRadius2_large;
bar_plot_test{2,1} = singleSynMeanHeadRadius4_small;
bar_plot_test{2,2} = singleSynMeanHeadRadius4_large;
bar_plot_test{3,1} = singleSynMeanHeadRadius5_small;
bar_plot_test{3,2} = singleSynMeanHeadRadius5_large;
hold on;
[position_cell, pvalueArray] = plot_significance_score_cell_cell(bar_plot_mean, bar_plot_test, b);
sigstar(position_cell,pvalueArray)
legend('L2/3', 'L4', 'L5')





% check within the same layer
singleSynMeanHeadRadius2_small = singleSynMeanHeadRadius2(spine_small_radius_id_2);
singleSynMeanHeadRadius2_small(singleSynMeanHeadRadius2_small == 0) = [];
singleSynMeanHeadRadius2_small(singleSynMeanHeadRadius2_small > quantile(singleSynMeanHeadRadius2_small, 0.99)) = [];
singleSynMeanHeadRadius4_small = singleSynMeanHeadRadius4(spine_small_radius_id_4);
singleSynMeanHeadRadius4_small(singleSynMeanHeadRadius4_small == 0) = [];
singleSynMeanHeadRadius4_small(singleSynMeanHeadRadius4_small > quantile(singleSynMeanHeadRadius4_small, 0.99)) = [];
singleSynMeanHeadRadius5_small = singleSynMeanHeadRadius5(spine_small_radius_id_5);
singleSynMeanHeadRadius5_small(singleSynMeanHeadRadius5_small == 0) = [];
singleSynMeanHeadRadius5_small(singleSynMeanHeadRadius5_small > quantile(singleSynMeanHeadRadius5_small, 0.99)) = [];
singleSynMeanHeadRadius2_large = singleSynMeanHeadRadius2(spine_large_radius_id_2);
singleSynMeanHeadRadius2_large(singleSynMeanHeadRadius2_large == 0) = [];
singleSynMeanHeadRadius2_large(singleSynMeanHeadRadius2_large > quantile(singleSynMeanHeadRadius2_large, 0.99)) = [];
singleSynMeanHeadRadius4_large = singleSynMeanHeadRadius4(spine_large_radius_id_4);
singleSynMeanHeadRadius4_large(singleSynMeanHeadRadius4_large == 0) = [];
singleSynMeanHeadRadius4_large(singleSynMeanHeadRadius4_large > quantile(singleSynMeanHeadRadius4_large, 0.99)) = [];
singleSynMeanHeadRadius5_large = singleSynMeanHeadRadius5(spine_large_radius_id_5);
singleSynMeanHeadRadius5_large(singleSynMeanHeadRadius5_large == 0) = [];
singleSynMeanHeadRadius5_large(singleSynMeanHeadRadius5_large > quantile(singleSynMeanHeadRadius5_large, 0.99)) = [];


plot_bar_plot_three_layers_two_group_withinlayers(singleSynMeanHeadRadius2_small,singleSynMeanHeadRadius2_large,...
    singleSynMeanHeadRadius4_small, singleSynMeanHeadRadius2_large, singleSynMeanHeadRadius5_small, singleSynMeanHeadRadius5_large)
ylabel('radius (nm)');  % Y-axis label for the values of the variablesneck
title('layer-specific head redius comparison (thick/ thin dendrite)');  % Chart title







%% neck thickness
neck_radius2 = [singleSynNeckMeanRadius2];
neck_radius4 = [singleSynNeckMeanRadius4];
neck_radius5 = [singleSynNeckMeanRadius5];

neck_radius2(neck_radius2 == 0 |isinf(neck_radius2)|isnan(neck_radius2)) = [];
neck_radius2(neck_radius2 > quantile(neck_radius2, 0.99)) = [];
neck_radius4(neck_radius4 == 0 |isinf(neck_radius4)|isnan(neck_radius4)) = [];
neck_radius4(neck_radius4 > quantile(neck_radius4, 0.99)) = [];
neck_radius5(neck_radius5 == 0 |isinf(neck_radius5)|isnan(neck_radius5)) = [];
neck_radius5(neck_radius5 > quantile(neck_radius5, 0.99)) = [];
% singleSynMeanHeadRadius5(singleSynMeanHeadRadius5>600) = [];
x = [neck_radius2(:);neck_radius4(:);neck_radius5(:)];

y = [zeros(length(neck_radius2(:)),1);ones(length(neck_radius4(:)),1);ones(length(neck_radius5(:)),1).*2];

figure; violinplot(x,y,'ShowMean',true,'ShowMedian',true); title(['dendrite spine neck radius'])
set(gca,'xtick',[1,2,3],'xticklabel',{'L2/3', 'L4', 'L5'}); ylabel('radius(nm)')
% figure; violinplot(x,y,'ShowMean',true,'ShowMedian',false); title(['dendrite spine neck section'])
% set(gca,'xtick',[1,2,3],'xticklabel',{'L2/3', 'L4', 'L5'})
% id2 = randperm(length(singleSynMeanHeadRadius2), 3000);
% id4 = randperm(length(singleSynMeanHeadRadius4), 3000);
% id5 = randperm(length(singleSynMeanHeadRadius5), 3000);
[h,p1,ci,stats] = ttest2(neck_radius2,neck_radius4);
[h,p2,ci,stats] = ttest2(neck_radius4,neck_radius5);
[h,p3,ci,stats] = ttest2(neck_radius2,neck_radius5);
sigstar({[1,2],[2,3], [1,3]},[p1,p2,p3])

% directly change to bar plot
neck_radius2 = [singleSynNeckMeanRadius2];
neck_radius4 = [singleSynNeckMeanRadius4];
neck_radius5 = [singleSynNeckMeanRadius5];

neck_radius2(neck_radius2 == 0 |isinf(neck_radius2)|isnan(neck_radius2)) = [];
neck_radius2(neck_radius2 > quantile(neck_radius2, 0.99)) = [];
neck_radius4(neck_radius4 == 0 |isinf(neck_radius4)|isnan(neck_radius4)) = [];
neck_radius4(neck_radius4 > quantile(neck_radius4, 0.99)) = [];
neck_radius5(neck_radius5 == 0 |isinf(neck_radius5)|isnan(neck_radius5)) = [];
neck_radius5(neck_radius5 > quantile(neck_radius5, 0.99)) = [];

data = [mean(neck_radius2), mean(neck_radius4), mean(neck_radius5)];
err = [[std(neck_radius2)/sqrt(length(neck_radius2)), std(neck_radius4)/sqrt(length(neck_radius4)), std(neck_radius5)/sqrt(length(neck_radius5))]];
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
bar_plot_test{1,1} = neck_radius2;
bar_plot_test{2,1} = neck_radius4;
bar_plot_test{3,1} = neck_radius5;
hold on;
[position_cell, pvalueArray] = plot_significance_score_cell_cell(bar_plot_mean, bar_plot_test, b);
sigstar(position_cell,pvalueArray)

ylim([0,180])




% bar plot of neck radius comparing between dendrite with large radius and small radius

singleSynNeckMeanRadius2_small = singleSynNeckMeanRadius2(spine_small_radius_id_2);
singleSynNeckMeanRadius2_small(singleSynNeckMeanRadius2_small == 0) = [];
singleSynNeckMeanRadius2_small(singleSynNeckMeanRadius2_small > quantile(singleSynNeckMeanRadius2_small, 0.99)) = [];
singleSynNeckMeanRadius2_small(isnan(singleSynNeckMeanRadius2_small)) = [];
singleSynNeckMeanRadius2_small(isinf(singleSynNeckMeanRadius2_small)) = [];

singleSynNeckMeanRadius4_small = singleSynNeckMeanRadius4(spine_small_radius_id_4);
singleSynNeckMeanRadius4_small(singleSynNeckMeanRadius4_small == 0) = [];
singleSynNeckMeanRadius4_small(singleSynNeckMeanRadius4_small > quantile(singleSynNeckMeanRadius4_small, 0.99)) = [];
singleSynNeckMeanRadius4_small(isnan(singleSynNeckMeanRadius4_small)) = [];
singleSynNeckMeanRadius4_small(isinf(singleSynNeckMeanRadius4_small)) = [];

singleSynNeckMeanRadius5_small = singleSynNeckMeanRadius5(spine_small_radius_id_5);
singleSynNeckMeanRadius5_small(singleSynNeckMeanRadius5_small == 0) = [];
singleSynNeckMeanRadius5_small(singleSynNeckMeanRadius5_small > quantile(singleSynNeckMeanRadius5_small, 0.99)) = [];
singleSynNeckMeanRadius5_small(isnan(singleSynNeckMeanRadius5_small)) = [];
singleSynNeckMeanRadius5_small(isinf(singleSynNeckMeanRadius5_small)) = [];

singleSynNeckMeanRadius2_large = singleSynNeckMeanRadius2(spine_large_radius_id_2);
singleSynNeckMeanRadius2_large(singleSynNeckMeanRadius2_large == 0) = [];
singleSynNeckMeanRadius2_large(singleSynNeckMeanRadius2_large > quantile(singleSynNeckMeanRadius2_large, 0.99)) = [];
singleSynNeckMeanRadius2_large(isnan(singleSynNeckMeanRadius2_large)) = [];
singleSynNeckMeanRadius2_large(isinf(singleSynNeckMeanRadius2_large)) = [];

singleSynNeckMeanRadius4_large = singleSynNeckMeanRadius4(spine_large_radius_id_4);
singleSynNeckMeanRadius4_large(singleSynNeckMeanRadius4_large == 0) = [];
singleSynNeckMeanRadius4_large(singleSynNeckMeanRadius4_large > quantile(singleSynNeckMeanRadius4_large, 0.99)) = [];
singleSynNeckMeanRadius4_large(isnan(singleSynNeckMeanRadius4_large)) = [];
singleSynNeckMeanRadius4_large(isinf(singleSynNeckMeanRadius4_large)) = [];

singleSynNeckMeanRadius5_large = singleSynNeckMeanRadius4(spine_large_radius_id_5);
singleSynNeckMeanRadius5_large(singleSynNeckMeanRadius5_large == 0) = [];
singleSynNeckMeanRadius5_large(singleSynNeckMeanRadius5_large > quantile(singleSynNeckMeanRadius5_large, 0.99)) = [];
singleSynNeckMeanRadius5_large(isnan(singleSynNeckMeanRadius5_large)) = [];
singleSynNeckMeanRadius5_large(isinf(singleSynNeckMeanRadius5_large)) = [];

data = [[mean(singleSynNeckMeanRadius2_small), mean(singleSynNeckMeanRadius4_small), mean(singleSynNeckMeanRadius5_small)];
    [mean(singleSynNeckMeanRadius2_large), mean(singleSynNeckMeanRadius4_large), mean(singleSynNeckMeanRadius5_large)]];
err = [[std(singleSynNeckMeanRadius2_small)/sqrt(length(singleSynNeckMeanRadius2_small)), std(singleSynNeckMeanRadius4_small)/sqrt(length(singleSynNeckMeanRadius4_small)), std(singleSynNeckMeanRadius5_small)/sqrt(length(singleSynNeckMeanRadius5_small))];
    [std(singleSynNeckMeanRadius2_large)/sqrt(length(singleSynNeckMeanRadius2_large)), std(singleSynNeckMeanRadius4_large)/sqrt(length(singleSynNeckMeanRadius4_large)), std(singleSynNeckMeanRadius5_large)/sqrt(length(singleSynNeckMeanRadius5_large))]];

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
bar_plot_test{1,1} = singleSynNeckMeanRadius2_small;
bar_plot_test{1,2} = singleSynNeckMeanRadius2_large;
bar_plot_test{2,1} = singleSynNeckMeanRadius4_small;
bar_plot_test{2,2} = singleSynNeckMeanRadius4_large;
bar_plot_test{3,1} = singleSynNeckMeanRadius5_small;
bar_plot_test{3,2} = singleSynNeckMeanRadius5_large;
hold on;
[position_cell, pvalueArray] = plot_significance_score_cell_cell(bar_plot_mean, bar_plot_test, b);
sigstar(position_cell,pvalueArray)
legend('L2/3', 'L4', 'L5')




% compare within each layer
singleSynNeckMeanRadius2_small = singleSynNeckMeanRadius2(spine_small_radius_id_2);
singleSynNeckMeanRadius2_small(singleSynNeckMeanRadius2_small == 0) = [];
singleSynNeckMeanRadius2_small(singleSynNeckMeanRadius2_small > quantile(singleSynNeckMeanRadius2_small, 0.99)) = [];
singleSynNeckMeanRadius2_small(isnan(singleSynNeckMeanRadius2_small)) = [];
singleSynNeckMeanRadius2_small(isinf(singleSynNeckMeanRadius2_small)) = [];

singleSynNeckMeanRadius4_small = singleSynNeckMeanRadius4(spine_small_radius_id_4);
singleSynNeckMeanRadius4_small(singleSynNeckMeanRadius4_small == 0) = [];
singleSynNeckMeanRadius4_small(singleSynNeckMeanRadius4_small > quantile(singleSynNeckMeanRadius4_small, 0.99)) = [];
singleSynNeckMeanRadius4_small(isnan(singleSynNeckMeanRadius4_small)) = [];
singleSynNeckMeanRadius4_small(isinf(singleSynNeckMeanRadius4_small)) = [];

singleSynNeckMeanRadius5_small = singleSynNeckMeanRadius5(spine_small_radius_id_5);
singleSynNeckMeanRadius5_small(singleSynNeckMeanRadius5_small == 0) = [];
singleSynNeckMeanRadius5_small(singleSynNeckMeanRadius5_small > quantile(singleSynNeckMeanRadius5_small, 0.99)) = [];
singleSynNeckMeanRadius5_small(isnan(singleSynNeckMeanRadius5_small)) = [];
singleSynNeckMeanRadius5_small(isinf(singleSynNeckMeanRadius5_small)) = [];

singleSynNeckMeanRadius2_large = singleSynNeckMeanRadius2(spine_large_radius_id_2);
singleSynNeckMeanRadius2_large(singleSynNeckMeanRadius2_large == 0) = [];
singleSynNeckMeanRadius2_large(singleSynNeckMeanRadius2_large > quantile(singleSynNeckMeanRadius2_large, 0.99)) = [];
singleSynNeckMeanRadius2_large(isnan(singleSynNeckMeanRadius2_large)) = [];
singleSynNeckMeanRadius2_large(isinf(singleSynNeckMeanRadius2_large)) = [];

singleSynNeckMeanRadius4_large = singleSynNeckMeanRadius4(spine_large_radius_id_4);
singleSynNeckMeanRadius4_large(singleSynNeckMeanRadius4_large == 0) = [];
singleSynNeckMeanRadius4_large(singleSynNeckMeanRadius4_large > quantile(singleSynNeckMeanRadius4_large, 0.99)) = [];
singleSynNeckMeanRadius4_large(isnan(singleSynNeckMeanRadius4_large)) = [];
singleSynNeckMeanRadius4_large(isinf(singleSynNeckMeanRadius4_large)) = [];

singleSynNeckMeanRadius5_large = singleSynNeckMeanRadius4(spine_large_radius_id_5);
singleSynNeckMeanRadius5_large(singleSynNeckMeanRadius5_large == 0) = [];
singleSynNeckMeanRadius5_large(singleSynNeckMeanRadius5_large > quantile(singleSynNeckMeanRadius5_large, 0.99)) = [];
singleSynNeckMeanRadius5_large(isnan(singleSynNeckMeanRadius5_large)) = [];
singleSynNeckMeanRadius5_large(isinf(singleSynNeckMeanRadius5_large)) = [];

plot_bar_plot_three_layers_two_group_withinlayers(singleSynNeckMeanRadius2_small, singleSynNeckMeanRadius2_large, ...
    singleSynNeckMeanRadius4_small, singleSynNeckMeanRadius4_large, singleSynNeckMeanRadius5_small, singleSynNeckMeanRadius5_large)

ylabel('radius (nm)');  % Y-axis label for the values of the variables
title('layer-specific neck radius comparison');  % Chart title





%% neck_section
neck_section2 = [singleSynNeckSection2];
neck_section4 = [singleSynNeckSection4];
neck_section5 = [singleSynNeckSection5];
neck_section2(neck_section2 == 0 |isinf(neck_section2)|isnan(neck_section2)) = [];
neck_section2(neck_section2 > quantile(neck_section2, 0.99)) = [];
neck_section4(neck_section4 == 0 |isinf(neck_section4)|isnan(neck_section4)) = [];
neck_section4(neck_section4 > quantile(neck_section4, 0.99)) = [];
neck_section5(neck_section5 == 0 |isinf(neck_section5)|isnan(neck_section5)) = [];
neck_section5(neck_section5 > quantile(neck_section5, 0.99)) = [];
% singleSynMeanHeadRadius5(singleSynMeanHeadRadius5>600) = [];
x = [neck_section2(:);neck_section4(:);neck_section5(:)];

y = [zeros(length(neck_section2(:)),1);ones(length(neck_section4(:)),1);ones(length(neck_section5(:)),1).*2];

figure; violinplot(x,y,'ShowMean',true,'ShowMedian',true); title(['dendrite spine neck section'])
set(gca,'xtick',[1,2,3],'xticklabel',{'L2/3', 'L4', 'L5'}); ylabel('Area(nm^2)')
% figure; violinplot(x,y,'ShowMean',true,'ShowMedian',false); title(['dendrite spine neck section'])
% set(gca,'xtick',[1,2,3],'xticklabel',{'L2/3', 'L4', 'L5'})
% id2 = randperm(length(singleSynMeanHeadRadius2), 3000);
% id4 = randperm(length(singleSynMeanHeadRadius4), 3000);
% id5 = randperm(length(singleSynMeanHeadRadius5), 3000);
[h,p1,ci,stats] = ttest2(neck_section2,neck_section4);
[h,p2,ci,stats] = ttest2(neck_section4,neck_section5);
[h,p3,ci,stats] = ttest2(neck_section2,neck_section5);
sigstar({[1,2],[2,3], [1,3]},[p1,p2,p3])



%% correlation between spine volume and synaptic cleft size
reset_values_summary_statistics;
id2 = find(singleSynHeadVolume2 > 0 & (singleSynapticCleftSize2 > 0) );
id4 = find(singleSynHeadVolume4 > 0 & (singleSynapticCleftSize4 > 0));
id5 = find(singleSynHeadVolume5 > 0 & (singleSynapticCleftSize5 > 0));
singleSynHeadVolume2 = singleSynHeadVolume2(id2).*16.*16.*40/10^9;
singleSynHeadVolume4 = singleSynHeadVolume4(id4).*16.*16.*40/10^9;
singleSynHeadVolume5 = singleSynHeadVolume5(id5).*16.*16.*40/10^9;
singleSynapticCleftSize2 = singleSynapticCleftSize2(id2).*16.*16.*40/10^9;
singleSynapticCleftSize4 = singleSynapticCleftSize4(id4).*16.*16.*40/10^9;
singleSynapticCleftSize5 = singleSynapticCleftSize5(id5).*16.*16.*40/10^9;
% singleSynapticCleftSize2(id2);singleSynapticCleftSize4(id4);singleSynapticCleftSize5(id5)
x1 = singleSynHeadVolume2;
x2 = singleSynHeadVolume4;
x3 = singleSynHeadVolume5;
y1 = singleSynapticCleftSize2;
y2 = singleSynapticCleftSize4;
y3 = singleSynapticCleftSize5;

figure; hold on; % Hold on to plot all on the same figure
plot(x1, y1, 'r.');
plot(x2, y2, 'g.');
plot(x3, y3, 'b.');

% Fitting regression lines
[p1,S1] = polyfit(x1, y1, 1);
[p2, S2] = polyfit(x2, y2, 1);
[p3, S3] = polyfit(x3, y3, 1);

% Calculating the regression line values
y1_fit = polyval(p1, x1);
y2_fit = polyval(p2, x2);
y3_fit = polyval(p3, x3);

% Plotting the regression lines
plot(x1, y1_fit, 'r-', 'LineWidth', 1.5);
plot(x2, y2_fit, 'g-', 'LineWidth', 1.5);
plot(x3, y3_fit, 'b-', 'LineWidth', 1.5);

% Adding legend and labels
legend('L2/3', 'L4', 'L5', 'Location', 'best');
xlabel('head volume (\mu m^3)');
ylabel('synaptic cleft volume (\mu m^3)');
title('dendrite spine head volume vs synaptic cleft volume');
grid on; % Optional: Add grid lines for better readability
hold off;



x1_range = linspace(min(x1), max(x1), 100);
y1_range = linspace(min(y1), max(y1), 100);
[X1, Y1] = meshgrid(x1_range, y1_range);    
density1 = ksdensity([x1(:), y1(:)], [X1(:), Y1(:)]);
density1 = reshape(density1, size(X1)); 

x2_range = linspace(min(x2), max(x2), 100);
y2_range = linspace(min(y2), max(y2), 100);
[X2, Y2] = meshgrid(x2_range, y2_range);
density2 = ksdensity([x2(:), y2(:)], [X2(:), Y2(:)]);
density2 = reshape(density2, size(X2));

x3_range = linspace(min(x3), max(x3), 100);
y3_range = linspace(min(y3), max(y3), 100);
[X3, Y3] = meshgrid(x3_range, y3_range);
density3 = ksdensity([x3(:), y3(:)], [X3(:), Y3(:)]);
density3 = reshape(density3, size(X3));


figure;
contour(X1, Y1, density1, 'LineWidth', 2, 'LineColor', 'r');
hold on;
contour(X2, Y2, density2, 'LineWidth', 2,'LineColor', 'g');
hold on;
contour(X3, Y3, density3, 'LineWidth', 2, 'LineColor', 'b');

xlim([0,0.6])
ylim([0, 0.02])
legend({'L2/3', 'L4', 'L5'}, 'Location', 'best')
xlabel('head volume (\mu m^3)');
ylabel('synaptic cleft volume (\mu m^3)');
title('dendrite spine head volume vs synaptic cleft volume');


%% head radius vs neck section radius
id2 = find(singleSynMeanHeadRadius2~=0);
id4 = find(singleSynMeanHeadRadius4~=0);
id5 = find(singleSynMeanHeadRadius5~=0);
singleSynMeanHeadRadius2 = singleSynMeanHeadRadius2(id2);
singleSynMeanHeadRadius4 = singleSynMeanHeadRadius4(id4);
singleSynMeanHeadRadius5 = singleSynMeanHeadRadius5(id5);
singleSynNeckSection2 = singleSynNeckSection2(id2);
singleSynNeckSection4 = singleSynNeckSection4(id4);
singleSynNeckSection5 = singleSynNeckSection5(id5);
figure; scatter(singleSynMeanHeadRadius2, sqrt(singleSynNeckSection2/pi))

hold on; scatter(singleSynMeanHeadRadius4, sqrt(singleSynNeckSection4/pi))
hold on; scatter(singleSynMeanHeadRadius5, sqrt(singleSynNeckSection5/pi))

ratio2 = singleSynMeanHeadRadius2(:)./sqrt(singleSynNeckSection2/pi);
ratio2(isinf(ratio2)) = [];
ratio2(ratio2 > quantile(ratio2 , 0.99)) = [];
ratio4 = singleSynMeanHeadRadius4(:)./sqrt(singleSynNeckSection4/pi);
ratio4(isinf(ratio4)) = [];
ratio4(ratio4 > quantile(ratio4 , 0.99)) = [];
ratio5 = singleSynMeanHeadRadius5(:)./sqrt(singleSynNeckSection5/pi);
ratio5(isinf(ratio5)) = [];
ratio5(ratio5 > quantile(ratio5 , 0.99)) = [];


x = [ratio2;ratio4;ratio5];
y = [zeros(length(ratio2),1);ones(length(ratio4),1);ones(length(ratio5),1).*2];
figure; violinplot(x,y,'ShowMean',true,'ShowMedian',true); title(['ratio between the radius of head and neck'])
set(gca,'xtick',[1,2,3],'xticklabel',{'L2/3', 'L4', 'L5'});

[h,p1,ci,stats] = ttest2(ratio2(:),ratio4(:));
[h,p2,ci,stats] = ttest2(ratio4(:),ratio5(:));
[h,p3,ci,stats] = ttest2(ratio2(:),ratio5(:));
sigstar({[1,2],[2,3], [1,3]},[p1,p2,p3])


% plot bar plot
id2 = find(singleSynMeanHeadRadius2~=0);
id4 = find(singleSynMeanHeadRadius4~=0);
id5 = find(singleSynMeanHeadRadius5~=0);
singleSynMeanHeadRadius2 = singleSynMeanHeadRadius2(id2);
singleSynMeanHeadRadius4 = singleSynMeanHeadRadius4(id4);
singleSynMeanHeadRadius5 = singleSynMeanHeadRadius5(id5);
singleSynNeckSection2 = singleSynNeckSection2(id2);
singleSynNeckSection4 = singleSynNeckSection4(id4);
singleSynNeckSection5 = singleSynNeckSection5(id5);
ratio2 = singleSynMeanHeadRadius2(:)./sqrt(singleSynNeckSection2/pi);
ratio2(isinf(ratio2)) = [];
ratio2(isnan(ratio2)) = [];
ratio2(ratio2 > quantile(ratio2 , 0.99)) = [];
ratio4 = singleSynMeanHeadRadius4(:)./sqrt(singleSynNeckSection4/pi);
ratio4(isinf(ratio4)) = [];
ratio4(isnan(ratio4)) = [];
ratio4(ratio4 > quantile(ratio4 , 0.99)) = [];
ratio5 = singleSynMeanHeadRadius5(:)./sqrt(singleSynNeckSection5/pi);
ratio5(isinf(ratio5)) = [];
ratio5(isnan(ratio5)) = [];
ratio5(ratio5 > quantile(ratio5 , 0.99)) = [];

data = [mean(ratio2), mean(ratio4), mean(ratio5)];
err = [[std(ratio2)/sqrt(length(ratio2)), std(ratio4)/sqrt(length(ratio4)), std(ratio5)/sqrt(length(ratio5))]];
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
% ylabel('radius (nm)');  % Y-axis label for the values of the variables
title('ratio between the radius of head and neck');  % Chart title
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
bar_plot_test{1,1} = ratio2;
bar_plot_test{2,1} = ratio4;
bar_plot_test{3,1} = ratio5;
hold on;
[position_cell, pvalueArray] = plot_significance_score_cell_cell(bar_plot_mean, bar_plot_test, b);
sigstar(position_cell,pvalueArray)


% bar plot of neck radius comparing between dendrite with large radius and small radius

ratio2 = singleSynMeanHeadRadius2(:)./sqrt(singleSynNeckSection2/pi);

ratio4 = singleSynMeanHeadRadius4(:)./sqrt(singleSynNeckSection4/pi);

ratio5 = singleSynMeanHeadRadius5(:)./sqrt(singleSynNeckSection5/pi);




ratio2_small = ratio2(spine_small_radius_id_2);
ratio2_small(ratio2_small == 0) = [];
ratio2_small(ratio2_small > quantile(ratio2_small, 0.99)) = [];
ratio2_small(isnan(ratio2_small)) = [];
ratio2_small(isinf(ratio2_small)) = [];

ratio4_small = ratio4(spine_small_radius_id_4);
ratio4_small(ratio4_small == 0) = [];
ratio4_small(ratio4_small > quantile(ratio4_small, 0.99)) = [];
ratio4_small(isnan(ratio4_small)) = [];
ratio4_small(isinf(ratio4_small)) = [];

ratio5_small = ratio5(spine_small_radius_id_5);
ratio5_small(ratio5_small == 0) = [];
ratio5_small(ratio5_small > quantile(ratio5_small, 0.99)) = [];
ratio5_small(isnan(ratio5_small)) = [];
ratio5_small(isinf(ratio5_small)) = [];

ratio2_large = ratio2(spine_large_radius_id_2);
ratio2_large(ratio2_large == 0) = [];
ratio2_large(ratio2_large > quantile(ratio2_large, 0.99)) = [];
ratio2_large(isnan(ratio2_large)) = [];
ratio2_large(isinf(ratio2_large)) = [];

ratio4_large = ratio4(spine_large_radius_id_4);
ratio4_large(ratio4_large == 0) = [];
ratio4_large(ratio4_large > quantile(ratio4_large, 0.99)) = [];
ratio4_large(isnan(ratio4_large)) = [];
ratio4_large(isinf(ratio4_large)) = [];

ratio5_large = ratio5(spine_large_radius_id_5);
ratio5_large(ratio5_large == 0) = [];
ratio5_large(ratio5_large > quantile(ratio5_large, 0.99)) = [];
ratio5_large(isnan(ratio5_large)) = [];
ratio5_large(isinf(ratio5_large)) = [];

data = [[mean(ratio2_small), mean(ratio4_small), mean(ratio5_small)];
    [mean(ratio2_large), mean(ratio4_large), mean(ratio5_large)]];
err = [[std(ratio2_small)/sqrt(length(ratio2_small)), std(ratio4_small)/sqrt(length(ratio4_small)), std(ratio5_small)/sqrt(length(ratio5_small))];
    [std(ratio2_large)/sqrt(length(ratio2_large)), std(ratio4_large)/sqrt(length(ratio4_large)), std(ratio5_large)/sqrt(length(ratio5_large))]];

figure;  % Creates a new figure
b = bar(data,'BarWidth', 0.8);
b(1).FaceColor = [0.8500, 0.3250, 0.0980];  % Color for Variable 1
b(2).FaceColor = [0, 0.4470, 0.7410];  % Color for Variable 2
b(3).FaceColor = [0.9290, 0.6940, 0.1250];  % Color for Variable 3% Plots grouped bar chart
% ylabel('radius (nm)');  % Y-axis label for the values of the variables
title('ratio between the radius of head and neck');  % Chart title
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
bar_plot_test{1,1} = ratio2_small;
bar_plot_test{1,2} = ratio2_large;
bar_plot_test{2,1} = ratio4_small;
bar_plot_test{2,2} = ratio4_large;
bar_plot_test{3,1} = ratio5_small;
bar_plot_test{3,2} = ratio5_large;
hold on;
[position_cell, pvalueArray] = plot_significance_score_cell_cell(bar_plot_mean, bar_plot_test, b);
sigstar(position_cell,pvalueArray)
legend('L2/3', 'L4', 'L5')

% layer-specific comparison between the ratio of the head and neck radius
    reset_values_summary_statistics
    
    ratio2 = singleSynMeanHeadRadius2(:)./sqrt(singleSynNeckSection2/pi);
    
    ratio4 = singleSynMeanHeadRadius4(:)./sqrt(singleSynNeckSection4/pi);
    
    ratio5 = singleSynMeanHeadRadius5(:)./sqrt(singleSynNeckSection5/pi);
    
    
    
    
    ratio2_small = ratio2(spine_small_radius_id_2);
    ratio2_small(ratio2_small == 0) = [];
    ratio2_small(ratio2_small > quantile(ratio2_small, 0.99)) = [];
    ratio2_small(isnan(ratio2_small)) = [];
    ratio2_small(isinf(ratio2_small)) = [];
    
    ratio4_small = ratio4(spine_small_radius_id_4);
    ratio4_small(ratio4_small == 0) = [];
    ratio4_small(ratio4_small > quantile(ratio4_small, 0.99)) = [];
    ratio4_small(isnan(ratio4_small)) = [];
    ratio4_small(isinf(ratio4_small)) = [];
    
    ratio5_small = ratio5(spine_small_radius_id_5);
    ratio5_small(ratio5_small == 0) = [];
    ratio5_small(ratio5_small > quantile(ratio5_small, 0.99)) = [];
    ratio5_small(isnan(ratio5_small)) = [];
    ratio5_small(isinf(ratio5_small)) = [];
    
    ratio2_large = ratio2(spine_large_radius_id_2);
    ratio2_large(ratio2_large == 0) = [];
    ratio2_large(ratio2_large > quantile(ratio2_large, 0.99)) = [];
    ratio2_large(isnan(ratio2_large)) = [];
    ratio2_large(isinf(ratio2_large)) = [];
    
    ratio4_large = ratio4(spine_large_radius_id_4);
    ratio4_large(ratio4_large == 0) = [];
    ratio4_large(ratio4_large > quantile(ratio4_large, 0.99)) = [];
    ratio4_large(isnan(ratio4_large)) = [];
    ratio4_large(isinf(ratio4_large)) = [];
    
    ratio5_large = ratio5(spine_large_radius_id_5);
    ratio5_large(ratio5_large == 0) = [];
    ratio5_large(ratio5_large > quantile(ratio5_large, 0.99)) = [];
    ratio5_large(isnan(ratio5_large)) = [];
    ratio5_large(isinf(ratio5_large)) = [];
    
    
    plot_bar_plot_three_layers_two_group_withinlayers(ratio2_small, ratio2_large, ...
        ratio4_small, ratio4_large, ratio5_small, ratio5_large)
    title('ratio between the radius of head and neck');










%% astrocyte head neck contact ratio
% plot the bar plot of the contact, conact in head, contact in neck
reset_values_summary_statistics;
contact = [sum(singleSynHeadNeckTouchingRatio2(:) > 0)/length(singleSynHeadNeckTouchingRatio2(:)),sum(singleSynHeadNeckTouchingRatio4(:) > 0)/length(singleSynHeadNeckTouchingRatio4(:)) ...
    sum(singleSynHeadNeckTouchingRatio5(:) > 0)/length(singleSynHeadNeckTouchingRatio5(:))];
headcontact = [sum(singleSynHeadNeckTouchingRatio2(:,1) >0)/length(singleSynHeadNeckTouchingRatio2(:,1)),sum(singleSynHeadNeckTouchingRatio4(:,1) > 0)/length(singleSynHeadNeckTouchingRatio4(:,1)) ...
    sum(singleSynHeadNeckTouchingRatio5(:,1) > 0)/length(singleSynHeadNeckTouchingRatio5(:,1))];
neckcontact = [sum(singleSynHeadNeckTouchingRatio2(:,2) >0)/length(singleSynHeadNeckTouchingRatio2(:,2)),sum(singleSynHeadNeckTouchingRatio4(:,2) > 0)/length(singleSynHeadNeckTouchingRatio4(:,2)) ...
    sum(singleSynHeadNeckTouchingRatio5(:,2) > 0)/length(singleSynHeadNeckTouchingRatio5(:,2))];
data = [contact;headcontact;neckcontact];

figure;  % Creates a new figure
b = bar(data, 'BarWidth', 0.8);
b(1).FaceColor = [0.8500, 0.3250, 0.0980];  % Color for Variable 1
b(2).FaceColor = [0, 0.4470, 0.7410];  % Color for Variable 2
b(3).FaceColor = [0.9290, 0.6940, 0.1250];  % Color for Variable 3% Plots grouped bar chart
% xlabel({'Contact', 'contact at head', 'contact at neck'});  % X-axis label for categories/groups
ylabel('ratio');  % Y-axis label for the values of the variables
title('Contact ratio at head and neck');  % Chart title
% legend({'L2/3', 'L4', 'L5'}, 'Location', 'best');  % Legend

groupNames = {'Contact', 'contact at head', 'contact at neck'}; % Custom names for each group of bars
xticks(1:length(groupNames)); % Set the x-axis ticks to match the number of groups
xticklabels(groupNames);


% extract the ratio per dendrite

reset_values_summary_statistics;
spine_dendrite_idx2 = label2idx(spine_dendrite_type_label2(:,2));
spine_dendrite_idx2(cellfun(@length, spine_dendrite_idx2) <= 3) = [];
spine_dendrite_idx2 = spine_dendrite_idx2(:);
contact2 = cellfun(@(c) sum(singleSynHeadNeckTouchingRatio2(c) > 0)/length(singleSynHeadNeckTouchingRatio2(c)), spine_dendrite_idx2);
headcontact2 = cellfun(@(c) sum(singleSynHeadNeckTouchingRatio2(c,1) > 0)/length(singleSynHeadNeckTouchingRatio2(c,1)), spine_dendrite_idx2);
neckcontact2 = cellfun(@(c) sum(singleSynHeadNeckTouchingRatio2(c,2) > 0)/length(singleSynHeadNeckTouchingRatio2(c,2)), spine_dendrite_idx2);

spine_dendrite_idx4 = label2idx(spine_dendrite_type_label4(:,2));
spine_dendrite_idx4(cellfun(@length, spine_dendrite_idx4) <= 3) = [];
spine_dendrite_idx4 = spine_dendrite_idx4(:);
contact4 = cellfun(@(c) sum(singleSynHeadNeckTouchingRatio4(c) > 0)/length(singleSynHeadNeckTouchingRatio4(c)), spine_dendrite_idx4);
headcontact4 = cellfun(@(c) sum(singleSynHeadNeckTouchingRatio4(c,1) > 0)/length(singleSynHeadNeckTouchingRatio4(c,1)), spine_dendrite_idx4);
neckcontact4 = cellfun(@(c) sum(singleSynHeadNeckTouchingRatio4(c,2) > 0)/length(singleSynHeadNeckTouchingRatio4(c,2)), spine_dendrite_idx4);

spine_dendrite_idx5 = label2idx(spine_dendrite_type_label5(:,2));
spine_dendrite_idx5(cellfun(@length, spine_dendrite_idx5) <= 3) = [];
spine_dendrite_idx5 = spine_dendrite_idx5(:);
contact5 = cellfun(@(c) sum(singleSynHeadNeckTouchingRatio5(c) > 0)/length(singleSynHeadNeckTouchingRatio5(c)), spine_dendrite_idx5);
headcontact5 = cellfun(@(c) sum(singleSynHeadNeckTouchingRatio5(c,1) > 0)/length(singleSynHeadNeckTouchingRatio5(c,1)), spine_dendrite_idx5);
neckcontact5 = cellfun(@(c) sum(singleSynHeadNeckTouchingRatio5(c,2) > 0)/length(singleSynHeadNeckTouchingRatio5(c,2)), spine_dendrite_idx5);

data = [[mean(contact2), mean(contact4), mean(contact5)];
    [mean(headcontact2), mean(headcontact4), mean(headcontact5)];
    [mean(neckcontact2), mean(neckcontact4), mean(neckcontact5)]]; % Data for the bar plot
err = [[std(contact2)/sqrt(length(contact2)), std(contact4)/sqrt(length(contact4)), std(contact5)/sqrt(length(contact5))];
    [std(headcontact2)/sqrt(length(headcontact2)), std(headcontact4)/sqrt(length(headcontact4)), std(headcontact5)/sqrt(length(headcontact5))];
    [std(neckcontact2)/sqrt(length(neckcontact2)), std(neckcontact4)/sqrt(length(neckcontact4)), std(neckcontact5)/sqrt(length(neckcontact5))]]; % Error bars

figure;  % Creates a new figure
b = bar(data, 'BarWidth', 0.8);
b(1).FaceColor = [0.8500, 0.3250, 0.0980];  % Color for Variable 1
b(2).FaceColor = [0, 0.4470, 0.7410];  % Color for Variable 2
b(3).FaceColor = [0.9290, 0.6940, 0.1250];  % Color for Variable 3% Plots grouped bar chart
% xlabel({'Contact', 'contact at head', 'contact at neck'});  % X-axis label for categories/groups
ylabel('ratio');  % Y-axis label for the values of the variables
title('Contact ratio at head and neck');  % Chart title
% legend({'L2/3', 'L4', 'L5'}, 'Location', 'best');  % Legend

groupNames = {'Contact', 'contact at head', 'contact at neck'}; % Custom names for each group of bars
xticks(1:length(groupNames)); % Set the x-axis ticks to match the number of groups
xticklabels(groupNames);
    
hold on;
numGroups = size(data, 1);
numBars = size(data, 2);

% Add error bars
for i = 1:numBars
    % X positions for error bars
    x = b(i).XEndPoints;
    % Add error bars to each group
    errorbar(x, data(:, i), err(:, i), '.k', 'LineWidth', 1.5);
end

types2_mean = [mean(contact2), mean(headcontact2), mean(neckcontact2)];
types4_mean = [mean(contact4), mean(headcontact4), mean(neckcontact4)];
types5_mean = [mean(contact5), mean(headcontact5), mean(neckcontact5)];
types2_reorg = [contact2, headcontact2, neckcontact2];
types4_reorg = [contact4, headcontact4, neckcontact4];
types5_reorg = [contact5, headcontact5, neckcontact5];
bar_plot_mean = cell(3,1);
bar_plot_mean{1} = types2_mean;
bar_plot_mean{2} = types4_mean;
bar_plot_mean{3} = types5_mean;
bar_plot_test = cell(3,1);
bar_plot_test{1} = types2_reorg;
bar_plot_test{2} = types4_reorg;
bar_plot_test{3} = types5_reorg;
hold on;
[position_cell, pvalueArray] = plot_significance_score(bar_plot_mean, bar_plot_test, b);
sigstar(position_cell,pvalueArray)

hold off;
legend({'L2/3', 'L4', 'L5'}, 'Location', 'best');  % Legend

plot_bar_plot_three_layers_two_group_withinlayers(headcontact2, neckcontact2, ...
     headcontact4, neckcontact4, ...
    headcontact5, neckcontact5);
legend({'contact at head', 'contact at neck'}, 'Location', 'best');  % Legend
title('Contact ratio at head and neck (within-layer comparison)');  % Chart title
        ylabel('ratio')


%% weigted wrapping area??
reset_values_summary_statistics;
sinsperimeterWeightedWrappingArea2(sinsperimeterWeightedWrappingArea2 == 0) = [];
sinsperimeterWeightedWrappingArea4(sinsperimeterWeightedWrappingArea4 == 0) = [];
sinsperimeterWeightedWrappingArea5(sinsperimeterWeightedWrappingArea5 == 0) = [];
% x = [sinsperimeterWeightedWrappingArea2(:);sinsperimeterWeightedWrappingArea4(:);sinsperimeterWeightedWrappingArea5(:)];
% y = [zeros(length(sinsperimeterWeightedWrappingArea2(:)),1);ones(length(sinsperimeterWeightedWrappingArea4(:)),1);ones(length(sinsperimeterWeightedWrappingArea5(:)),1).*2];
% figure; violinplot(x,y,'ShowMean',true,'ShowMedian',true); title(['perimeter ratio'])
% set(gca,'xtick',[1,2,3],'xticklabel',{'L2/3', 'L4', 'L5'})
% [h,p1,ci,stats] = ttest2(sinsperimeterWeightedWrappingArea2,sinsperimeterWeightedWrappingArea4);
% [h,p2,ci,stats] = ttest2(sinsperimeterWeightedWrappingArea4,sinsperimeterWeightedWrappingArea5);
% [h,p3,ci,stats] = ttest2(sinsperimeterWeightedWrappingArea2,sinsperimeterWeightedWrappingArea5);

plot_bar_plot_three_layers(sinsperimeterWeightedWrappingArea2, sinsperimeterWeightedWrappingArea4, sinsperimeterWeightedWrappingArea5)

title('weighted wrapping area (5 \mum range)')
ylabel('number of voxels')
        
        
%% perimeter touching ratio
reset_values_summary_statistics;
sinsperimeterRatio2(sinsperimeterRatio2 == 0) = [];
sinsperimeterRatio4(sinsperimeterRatio4 == 0) = [];
sinsperimeterRatio5(sinsperimeterRatio5 == 0) = [];
x = [sinsperimeterRatio2(:);sinsperimeterRatio4(:);sinsperimeterRatio5(:)];
y = [zeros(length(sinsperimeterRatio2(:)),1);ones(length(sinsperimeterRatio4(:)),1);ones(length(sinsperimeterRatio5(:)),1).*2];
figure; violinplot(x,y,'ShowMean',true,'ShowMedian',true); title(['perimeter ratio'])
set(gca,'xtick',[1,2,3],'xticklabel',{'L2/3', 'L4', 'L5'})
[h,p1,ci,stats] = ttest2(sinsperimeterRatio2,sinsperimeterRatio4);
[h,p2,ci,stats] = ttest2(sinsperimeterRatio4,sinsperimeterRatio5);
[h,p3,ci,stats] = ttest2(sinsperimeterRatio2,sinsperimeterRatio5);

plot_bar_plot_three_layers(sinsperimeterRatio2, sinsperimeterRatio4, sinsperimeterRatio5)


reset_values_summary_statistics;
sinsperimeterRatio2_small = sinsperimeterRatio2(spine_small_radius_id_2);
sinsperimeterRatio2_small(sinsperimeterRatio2_small == 0) = [];
sinsperimeterRatio4_small = sinsperimeterRatio4(spine_small_radius_id_4);
sinsperimeterRatio4_small(sinsperimeterRatio4_small == 0) = [];
sinsperimeterRatio5_small = sinsperimeterRatio5(spine_small_radius_id_5);
sinsperimeterRatio5_small(sinsperimeterRatio5_small == 0) = [];
sinsperimeterRatio2_large = sinsperimeterRatio2(spine_large_radius_id_2);
sinsperimeterRatio2_large(sinsperimeterRatio2_large == 0) = [];
sinsperimeterRatio4_large = sinsperimeterRatio4(spine_large_radius_id_4);
sinsperimeterRatio4_large(sinsperimeterRatio4_large == 0) = [];
sinsperimeterRatio5_large = sinsperimeterRatio5(spine_large_radius_id_5);
sinsperimeterRatio5_large(sinsperimeterRatio5_large == 0) = [];

plot_bar_plot_three_layers_two_group(sinsperimeterRatio2_small, sinsperimeterRatio2_large, sinsperimeterRatio4_small, ...
    sinsperimeterRatio4_large, sinsperimeterRatio5_small, sinsperimeterRatio5_large)
title('perimeter ratio in thick and thin dendrites')
ylabel('ratio')

% layer-specific comparison of the perimeter ratio
reset_values_summary_statistics;
sinsperimeterRatio2_small = sinsperimeterRatio2(spine_small_radius_id_2);
sinsperimeterRatio2_small(sinsperimeterRatio2_small == 0) = [];
sinsperimeterRatio4_small = sinsperimeterRatio4(spine_small_radius_id_4);
sinsperimeterRatio4_small(sinsperimeterRatio4_small == 0) = [];
sinsperimeterRatio5_small = sinsperimeterRatio5(spine_small_radius_id_5);
sinsperimeterRatio5_small(sinsperimeterRatio5_small == 0) = [];
sinsperimeterRatio2_large = sinsperimeterRatio2(spine_large_radius_id_2);
sinsperimeterRatio2_large(sinsperimeterRatio2_large == 0) = [];
sinsperimeterRatio4_large = sinsperimeterRatio4(spine_large_radius_id_4);
sinsperimeterRatio4_large(sinsperimeterRatio4_large == 0) = [];
sinsperimeterRatio5_large = sinsperimeterRatio5(spine_large_radius_id_5);
sinsperimeterRatio5_large(sinsperimeterRatio5_large == 0) = [];

plot_bar_plot_three_layers_two_group_withinlayers(sinsperimeterRatio2_small, sinsperimeterRatio2_large, sinsperimeterRatio4_small, ...
    sinsperimeterRatio4_large, sinsperimeterRatio5_small, sinsperimeterRatio5_large)
title('perimeter ratio (within-layer comparison)')
ylabel('ratio')

%% Pre synapse contact ratio
sinspreSynapseTouchingRatio2(sinspreSynapseTouchingRatio2 == 0) = [];
sinspreSynapseTouchingRatio4(sinspreSynapseTouchingRatio4 == 0) = [];
sinspreSynapseTouchingRatio5(sinspreSynapseTouchingRatio5 == 0) = [];

x = [sinspreSynapseTouchingRatio2(:);sinspreSynapseTouchingRatio4(:);sinspreSynapseTouchingRatio5(:)];
y = [zeros(length(sinspreSynapseTouchingRatio2(:)),1);ones(length(sinspreSynapseTouchingRatio4(:)),1);ones(length(sinspreSynapseTouchingRatio5(:)),1).*2];

figure; violinplot(x,y,'ShowMean',true,'ShowMedian',true); title(['Pre synapse contact ratio'])
set(gca,'xtick',[1,2,3],'xticklabel',{'L2/3', 'L4', 'L5'})
[h,p1,ci,stats] = ttest2(sinspreSynapseTouchingRatio2(:),sinspreSynapseTouchingRatio4(:));
[h,p2,ci,stats] = ttest2(sinspreSynapseTouchingRatio4(:),sinspreSynapseTouchingRatio5(:));
[h,p3,ci,stats] = ttest2(sinspreSynapseTouchingRatio2(:),sinspreSynapseTouchingRatio5(:));
sigstar({[1,2],[2,3], [1,3]},[p1,p2,p3])



%barplot

plot_bar_plot_three_layers(sinspreSynapseTouchingRatio2, sinspreSynapseTouchingRatio4, sinspreSynapseTouchingRatio5)
ylabel('ratio');title('pre synapse contact ratio')

reset_values_summary_statistics

sinspreSynapseTouchingRatio2_small = sinspreSynapseTouchingRatio2(spine_small_radius_id_2);
sinspreSynapseTouchingRatio2_small(sinspreSynapseTouchingRatio2_small == 0) = [];
sinspreSynapseTouchingRatio4_small = sinspreSynapseTouchingRatio4(spine_small_radius_id_4);
sinspreSynapseTouchingRatio4_small(sinspreSynapseTouchingRatio4_small == 0) = [];
sinspreSynapseTouchingRatio5_small = sinspreSynapseTouchingRatio5(spine_small_radius_id_5);
sinspreSynapseTouchingRatio5_small(sinspreSynapseTouchingRatio5_small == 0) = [];

sinspreSynapseTouchingRatio2_large = sinspreSynapseTouchingRatio2(spine_large_radius_id_2);
sinspreSynapseTouchingRatio2_large(sinspreSynapseTouchingRatio2_large == 0) = [];
sinspreSynapseTouchingRatio4_large = sinspreSynapseTouchingRatio4(spine_large_radius_id_4);
sinspreSynapseTouchingRatio4_large(sinspreSynapseTouchingRatio4_large == 0) = [];
sinspreSynapseTouchingRatio5_large = sinspreSynapseTouchingRatio5(spine_large_radius_id_5);
sinspreSynapseTouchingRatio5_large(sinspreSynapseTouchingRatio5_large == 0) = [];

plot_bar_plot_three_layers_two_group(sinspreSynapseTouchingRatio2_small, sinspreSynapseTouchingRatio2_large, sinspreSynapseTouchingRatio4_small, ...
    sinspreSynapseTouchingRatio4_large, sinspreSynapseTouchingRatio5_small, sinspreSynapseTouchingRatio5_large)
ylabel('ratio');title('pre synapse contact ratio in thin and thick dendrites')


plot_bar_plot_three_layers_two_group_withinlayers(sinspreSynapseTouchingRatio2_small, sinspreSynapseTouchingRatio2_large, sinspreSynapseTouchingRatio4_small, ...
    sinspreSynapseTouchingRatio4_large, sinspreSynapseTouchingRatio5_small, sinspreSynapseTouchingRatio5_large)
ylabel('ratio');title('pre synapse contact ratio (within-layer comparison)')



%% Post Synapse touching ratio
sinspostSynapseTouchingRatio2(sinspostSynapseTouchingRatio2 == 0) = [];
sinspostSynapseTouchingRatio4(sinspostSynapseTouchingRatio4 == 0)= [];
sinspostSynapseTouchingRatio5(sinspostSynapseTouchingRatio5 == 0) = [];
x = [sinspostSynapseTouchingRatio2(:);sinspostSynapseTouchingRatio4(:);sinspostSynapseTouchingRatio5(:)];
y = [zeros(length(sinspostSynapseTouchingRatio2(:)),1);ones(length(sinspostSynapseTouchingRatio4(:)),1);ones(length(sinspostSynapseTouchingRatio5(:)),1).*2];

figure; violinplot(x,y,'ShowMean',true,'ShowMedian',true); title(['post-synapse contact ratio'])
set(gca,'xtick',[1,2,3],'xticklabel',{'L2/3', 'L4', 'L5'})
[h,p1,ci,stats] = ttest2(sinspostSynapseTouchingRatio2,sinspostSynapseTouchingRatio4);
[h,p2,ci,stats] = ttest2(sinspostSynapseTouchingRatio4,sinspostSynapseTouchingRatio5);
[h,p3,ci,stats] = ttest2(sinspostSynapseTouchingRatio2,sinspostSynapseTouchingRatio5);
sigstar({[1,2],[2,3], [1,3]},[p1,p2,p3])

%barplot

plot_bar_plot_three_layers(sinspostSynapseTouchingRatio2, sinspostSynapseTouchingRatio4, sinspostSynapseTouchingRatio5)
ylabel('ratio');title('post synapse contact ratio')
ylim([0,0.32])

reset_values_summary_statistics

sinspostSynapseTouchingRatio2_small = sinspostSynapseTouchingRatio2(spine_small_radius_id_2);
sinspostSynapseTouchingRatio2_small(sinspostSynapseTouchingRatio2_small == 0) = [];
sinspostSynapseTouchingRatio4_small = sinspostSynapseTouchingRatio4(spine_small_radius_id_4);
sinspostSynapseTouchingRatio4_small(sinspostSynapseTouchingRatio4_small == 0) = [];
sinspostSynapseTouchingRatio5_small = sinspostSynapseTouchingRatio5(spine_small_radius_id_5);
sinspostSynapseTouchingRatio5_small(sinspostSynapseTouchingRatio5_small == 0) = [];

sinspostSynapseTouchingRatio2_large = sinspostSynapseTouchingRatio2(spine_large_radius_id_2);
sinspostSynapseTouchingRatio2_large(sinspostSynapseTouchingRatio2_large == 0) = [];
sinspostSynapseTouchingRatio4_large = sinspostSynapseTouchingRatio4(spine_large_radius_id_4);
sinspostSynapseTouchingRatio4_large(sinspostSynapseTouchingRatio4_large == 0) = [];
sinspostSynapseTouchingRatio5_large = sinspostSynapseTouchingRatio5(spine_large_radius_id_5);
sinspostSynapseTouchingRatio5_large(sinspostSynapseTouchingRatio5_large == 0) = [];

plot_bar_plot_three_layers_two_group(sinspostSynapseTouchingRatio2_small, sinspostSynapseTouchingRatio2_large, sinspostSynapseTouchingRatio4_small, ...
    sinspostSynapseTouchingRatio4_large, sinspostSynapseTouchingRatio5_small, sinspostSynapseTouchingRatio5_large)
ylabel('ratio');title('post synapse contact ratio in thin and thick dendrites')

plot_bar_plot_three_layers_two_group_withinlayers(sinspostSynapseTouchingRatio2_small, sinspostSynapseTouchingRatio2_large, sinspostSynapseTouchingRatio4_small, ...
    sinspostSynapseTouchingRatio4_large, sinspostSynapseTouchingRatio5_small, sinspostSynapseTouchingRatio5_large)
title('post synapse contact ratio (within-layer comparison)')
ylabel('ratio');


%% ratio of the contact area between head and neck
singleSynHeadNeckTouchingArea2 = singleSynHeadNeckTouchingArea2((singleSynHeadNeckTouchingArea2(:,1) > 10) & (singleSynHeadNeckTouchingArea2(:,2) > 10), :);
singleSynHeadNeckTouchingArea4 = singleSynHeadNeckTouchingArea4((singleSynHeadNeckTouchingArea4(:,1) > 10) & (singleSynHeadNeckTouchingArea4(:,2) > 10), :);
singleSynHeadNeckTouchingArea5 = singleSynHeadNeckTouchingArea5((singleSynHeadNeckTouchingArea5(:,1) > 10) & (singleSynHeadNeckTouchingArea5(:,2) > 10), :);

ratio2 = singleSynHeadNeckTouchingArea2(:,1)./singleSynHeadNeckTouchingArea2(:,2);
ratio4 = singleSynHeadNeckTouchingArea4(:,1)./singleSynHeadNeckTouchingArea4(:,2);
ratio5 = singleSynHeadNeckTouchingArea5(:,1)./singleSynHeadNeckTouchingArea5(:,2);
ratio2(ratio2 > quantile(ratio2, 0.99)) = [];
ratio4(ratio4 > quantile(ratio4, 0.99)) = [];
ratio5(ratio5 > quantile(ratio5, 0.99)) = [];
x = [ratio2(:);ratio4(:);ratio5(:)];
x = log(x);
y = [zeros(length(ratio2(:)),1);ones(length(ratio4(:)),1);ones(length(ratio5(:)),1).*2];

figure; violinplot(x,y,'ShowMean',true,'ShowMedian',true); title(['ratio of the contact area between head and neck'])
set(gca,'xtick',[1,2,3],'xticklabel',{'L2/3', 'L4', 'L5'})
% id2 = randperm(length(sinspostSynapseTouchingRatio2), 3000);
% id4 = randperm(length(sinspostSynapseTouchingRatio4), 3000);
% id5 = randperm(length(sinspostSynapseTouchingRatio5), 3000);
[h,p1,ci,stats] = ttest2(ratio2,ratio4);
[h,p2,ci,stats] = ttest2(ratio4,ratio5);
[h,p3,ci,stats] = ttest2(ratio2,ratio5);
sigstar({[1,2],[2,3], [1,3]},[p1,p2,p3])

plot_bar_plot_three_layers(ratio2, ratio4, ratio5)
ylabel('ratio');title('ratio of the contact area between head and neck')



reset_values_summary_statistics

singleSynHeadNeckTouchingArea2_small = singleSynHeadNeckTouchingArea2(spine_small_radius_id_2,:);
singleSynHeadNeckTouchingArea2_small = singleSynHeadNeckTouchingArea2_small((singleSynHeadNeckTouchingArea2_small(:,1) > 10) & (singleSynHeadNeckTouchingArea2_small(:,2) > 10), :);
ratio2_small = singleSynHeadNeckTouchingArea2_small(:,1)./singleSynHeadNeckTouchingArea2_small(:,2);
ratio2_small(ratio2_small > quantile(ratio2_small, 0.99)) = [];
singleSynHeadNeckTouchingArea4_small = singleSynHeadNeckTouchingArea4(spine_small_radius_id_4,:);
singleSynHeadNeckTouchingArea4_small = singleSynHeadNeckTouchingArea4_small((singleSynHeadNeckTouchingArea4_small(:,1) > 10) & (singleSynHeadNeckTouchingArea4_small(:,2) > 10), :);
ratio4_small = singleSynHeadNeckTouchingArea4_small(:,1)./singleSynHeadNeckTouchingArea4_small(:,2);
ratio4_small(ratio4_small > quantile(ratio4_small, 0.99)) = [];
singleSynHeadNeckTouchingArea5_small = singleSynHeadNeckTouchingArea5(spine_small_radius_id_5,:);
singleSynHeadNeckTouchingArea5_small = singleSynHeadNeckTouchingArea5_small((singleSynHeadNeckTouchingArea5_small(:,1) > 10) & (singleSynHeadNeckTouchingArea5_small(:,2) > 10), :);
ratio5_small = singleSynHeadNeckTouchingArea5_small(:,1)./singleSynHeadNeckTouchingArea5_small(:,2);
ratio5_small(ratio5_small > quantile(ratio5_small, 0.99)) = [];

singleSynHeadNeckTouchingArea2_large = singleSynHeadNeckTouchingArea2(spine_large_radius_id_2,:);
singleSynHeadNeckTouchingArea2_large = singleSynHeadNeckTouchingArea2_large((singleSynHeadNeckTouchingArea2_large(:,1) > 10) & (singleSynHeadNeckTouchingArea2_large(:,2) > 10), :);
ratio2_large = singleSynHeadNeckTouchingArea2_large(:,1)./singleSynHeadNeckTouchingArea2_large(:,2);
ratio2_large(ratio2_large > quantile(ratio2_large, 0.99)) = [];
singleSynHeadNeckTouchingArea4_large = singleSynHeadNeckTouchingArea4(spine_large_radius_id_4,:);
singleSynHeadNeckTouchingArea4_large = singleSynHeadNeckTouchingArea4_large((singleSynHeadNeckTouchingArea4_large(:,1) > 10) & (singleSynHeadNeckTouchingArea4_large(:,2) > 10), :);
ratio4_large = singleSynHeadNeckTouchingArea4_large(:,1)./singleSynHeadNeckTouchingArea4_large(:,2);
ratio4_large(ratio4_large > quantile(ratio4_large, 0.99)) = [];
singleSynHeadNeckTouchingArea5_large = singleSynHeadNeckTouchingArea5(spine_large_radius_id_5,:);
singleSynHeadNeckTouchingArea5_large = singleSynHeadNeckTouchingArea5_large((singleSynHeadNeckTouchingArea5_large(:,1) > 10) & (singleSynHeadNeckTouchingArea5_large(:,2) > 10), :);
ratio5_large = singleSynHeadNeckTouchingArea5_large(:,1)./singleSynHeadNeckTouchingArea5_large(:,2);
ratio5_large(ratio5_large > quantile(ratio5_large, 0.99)) = [];

plot_bar_plot_three_layers_two_group(ratio2_small, ratio2_large, ratio4_small, ...
    ratio4_large, ratio5_small, ratio5_large)

ylabel('ratio');title('ratio of the contact area between head and neck in thick and thin dendrites')

% layer-specific comparison of the ratio that contact at head/ contact at neck
plot_bar_plot_three_layers_two_group_withinlayers(ratio2_small, ratio2_large, ratio4_small, ...
    ratio4_large, ratio5_small, ratio5_large)
title('ratio of the contact area between head and neck (within-layer comparison)')
ylabel('ratio')



figure; hs = histogram(log2(ratio2), 50,'Normalization','probability');
hold on; histogram(log2(ratio4), 50,'Normalization','probability', 'BinLimits',hs.BinLimits);
hold on; histogram(log2(ratio5), 50,'Normalization','probability','BinLimits',hs.BinLimits);
xlabel('log-ratio'); ylabel('frequency')
legend('L2/3', 'L4', 'L5')
title('Ratio of the contact area between head and neck')
%% correlation between neck length and dendrite radius

reset_values_summary_statistics;
spine_dendrite_idx2 = label2idx(spine_dendrite_type_label2(:,2));
dendrite_radius2(cellfun(@length, spine_dendrite_idx2) <= 3) = [];
spine_dendrite_idx2(cellfun(@length, spine_dendrite_idx2) <= 3) = [];
mean_neck_length_2 = cellfun(@(c) mean(singleSynNeckLength2(c), "omitnan"), spine_dendrite_idx2);
id2 = find(mean_neck_length_2 > 0 & mean_neck_length_2 < quantile(mean_neck_length_2, 0.99));
mean_neck_length_2 = mean_neck_length_2(id2);
dendrite_radius2 = dendrite_radius2(id2);

spine_dendrite_idx4 = label2idx(spine_dendrite_type_label4(:,2));
dendrite_radius4(cellfun(@length, spine_dendrite_idx4) <= 3) = [];
spine_dendrite_idx4(cellfun(@length, spine_dendrite_idx4) <= 3) = [];
mean_neck_length_4 = cellfun(@(c) mean(singleSynNeckLength4(c), "omitnan"), spine_dendrite_idx4);
id4 = find(mean_neck_length_4 > 0 & mean_neck_length_4 < quantile(mean_neck_length_4, 0.99));
mean_neck_length_4 = mean_neck_length_4(id4);
dendrite_radius4 = dendrite_radius4(id4);

spine_dendrite_idx5 = label2idx(spine_dendrite_type_label5(:,2));
dendrite_radius5(cellfun(@length, spine_dendrite_idx5) <= 3) = [];
spine_dendrite_idx5(cellfun(@length, spine_dendrite_idx5) <= 3) = [];
mean_neck_length_5 = cellfun(@(c) mean(singleSynNeckLength5(c), "omitnan"), spine_dendrite_idx5);
id5 = find(mean_neck_length_5 > 0 & mean_neck_length_5 < quantile(mean_neck_length_5, 0.99));
mean_neck_length_5 = mean_neck_length_5(id5);
dendrite_radius5 = dendrite_radius5(id5);

plot_kde_contour_3(dendrite_radius2, mean_neck_length_2, dendrite_radius4, mean_neck_length_4, dendrite_radius5, mean_neck_length_5)
title('joint distribution of dendrite radius and neck length')

%% correlation between head volume and contact ratio
reset_values_summary_statistics;
id2 = find((singleSynHeadVolume2 ~=0) & (singleSynapticCleftSize2 > 0));
id4 = find((singleSynHeadVolume4 ~= 0) & (singleSynapticCleftSize4 > 0));
id5 = find((singleSynHeadVolume5 ~= 0) & (singleSynapticCleftSize5 > 0));
singleSynHeadVolume2 = singleSynHeadVolume2(id2).*16.*16.*40/10^9;
singleSynHeadVolume4 = singleSynHeadVolume4(id4).*16.*16.*40/10^9;
singleSynHeadVolume5 = singleSynHeadVolume5(id5).*16.*16.*40/10^9;
singleSynapticCleftSize2 = singleSynapticCleftSize2(id2).*16.*16.*40/10^9;
singleSynapticCleftSize4 = singleSynapticCleftSize4(id4).*16.*16.*40/10^9;
singleSynapticCleftSize5 = singleSynapticCleftSize5(id5).*16.*16.*40/10^9;
sinsperimeterRatio2 = sinsperimeterRatio2(id2);
sinsperimeterRatio4 = sinsperimeterRatio4(id4);
sinsperimeterRatio5 = sinsperimeterRatio5(id5);
singleSynHeadNeckTouchingArea2 = singleSynHeadNeckTouchingArea2(id2,:).*16.*16.*40/10^9;
singleSynHeadNeckTouchingArea4 = singleSynHeadNeckTouchingArea4(id4,:).*16.*16.*40/10^9;
singleSynHeadNeckTouchingArea5 = singleSynHeadNeckTouchingArea5(id5,:).*16.*16.*40/10^9;

x1 = (singleSynapticCleftSize2.^(1/3));
x2 = (singleSynapticCleftSize4.^(1/3));
x3 = (singleSynapticCleftSize5.^(1/3));
y1 = (singleSynHeadNeckTouchingArea2(:,1) + singleSynHeadNeckTouchingArea2(:,2)).^(1/3);
y2 = (singleSynHeadNeckTouchingArea4(:,1) + singleSynHeadNeckTouchingArea4(:,2)).^(1/3);
y3 = (singleSynHeadNeckTouchingArea5(:,1) + singleSynHeadNeckTouchingArea5(:,2)).^(1/3);
x1(isnan(y1)) = [];
y1(isnan(y1)) = [];
id1 = find(y1> quantile(y1, 0.01) & y1 < quantile(y1, 0.99));
x2(isnan(y2)) = [];
y2(isnan(y2)) = [];
id2 = find(y2> quantile(y2, 0.01) & y2 < quantile(y2, 0.99));
x3(isnan(y3)) = [];
y3(isnan(y3)) = [];
id3 = find(y3 > quantile(y3, 0.01) & y3 < quantile(y3, 0.99));
plot_kde_contour_3(x1(id1), y1(id1), x2(id2), y2(id2), x3(id3), y3(id3))


plot_regression_L2L4L5(x1(id1), x2(id2), x3(id3), y1(id1), y2(id2), y3(id3))
%% Ratio of each type across layers grouped by each dendrite
% 1: filopodia, 2: mushroom, 3: long thin, 4: thin, 5: stubby, 6:branched

id5 = find(dendrite_length5 >0 & sum(dendrite_each_spine_type_number5,2) > 3);
dendrite_length5 = dendrite_length5(id5);
dendrite_each_spine_type_number5 = dendrite_each_spine_type_number5(id5,:);
types5 = dendrite_each_spine_type_number5./dendrite_length5.*1000;
std5 = std(types5, 1)./sqrt(size(types5,1));
types5_mean = mean(types5,1);


id4 = find(dendrite_length4 > 0 & sum(dendrite_each_spine_type_number4,2) > 3);
dendrite_length4 = dendrite_length4(id4);
dendrite_each_spine_type_number4 = dendrite_each_spine_type_number4(id4,:);
types4 = dendrite_each_spine_type_number4./dendrite_length4.*1000;
% types4 = types4./sum(types4, 2);
std4 = std(types4, 1)./sqrt(size(types4,1));
types4_mean = mean(types4,1);

id2 = find(dendrite_length2 >0 & sum(dendrite_each_spine_type_number2,2) > 3);
dendrite_length2 = dendrite_length2(id2);
dendrite_each_spine_type_number2 = dendrite_each_spine_type_number2(id2,:);
types2 = dendrite_each_spine_type_number2./dendrite_length2.*1000;
% types2 = types2./sum(types2, 2);
std2 = std(types2, 1)./sqrt(size(types2,1));
types2_mean = mean(types2,1);

filo = [types2_mean(1), types4_mean(1),types5_mean(1)];
longthin = [types2_mean(3), types4_mean(3),types5_mean(3)];
thin = [types2_mean(4), types4_mean(4),types5_mean(4)];
stubby = [types2_mean(5), types4_mean(5),types5_mean(5)];
mush = [types2_mean(2), types4_mean(2),types5_mean(2)];
branched = [types2_mean(6), types4_mean(6),types5_mean(6)];
data = [filo;longthin;thin;stubby;mush;branched];

filo_err = [std2(1), std4(1),std5(1)];
longthin_err = [std2(3), std4(3),std5(3)];
thin_err = [std2(4), std4(4),std5(4)];
stubby_err = [std2(5), std4(5),std5(5)];
mush_err = [std2(2), std4(2),std5(2)];
branched_err = [std2(6), std4(6),std5(6)];
err = [filo_err;longthin_err;thin_err;stubby_err;mush_err;branched_err];
figure;  % Creates a new figure
b = bar(data,'BarWidth', 0.8);
b(1).FaceColor = [0.8500, 0.3250, 0.0980];  % Color for Variable 1
b(2).FaceColor = [0, 0.4470, 0.7410];  % Color for Variable 2
b(3).FaceColor = [0.9290, 0.6940, 0.1250];  % Color for Variable 3% Plots grouped bar chart
% xlabel({'Contact', 'contact at head', 'contact at neck'});  % X-axis label for categories/groups
ylabel('Density \mum');  % Y-axis label for the values of the variables
title('Density of each spine type per dendrite');  % Chart title


groupNames = {'Filopodia', 'Long Thin', 'Thin', 'Stubby', 'Mushroom', 'Branched'}; % Custom names for each group of bars
xticks(1:length(groupNames)); % Set the x-axis ticks to match the number of groups
xticklabels(groupNames);

hold on;
numGroups = size(data, 1);
numBars = size(data, 2);

% Add error bars
for i = 1:numBars
    % X positions for error bars
    x = b(i).XEndPoints;
    % Add error bars to each group
    errorbar(x, data(:, i), err(:, i), '.k', 'LineWidth', 1.5);
end
types5_reorg = [types5(:,1),types5(:,3), types5(:,4), types5(:,5), types5(:,2), types5(:,6)];
types4_reorg = [types4(:,1),types4(:,3), types4(:,4), types4(:,5), types4(:,2), types4(:,6)];
types2_reorg = [types2(:,1),types2(:,3), types2(:,4), types2(:,5), types2(:,2), types2(:,6)];

bar_plot_mean = cell(3,1);
bar_plot_mean{1} = types2_mean;
bar_plot_mean{2} = types4_mean;
bar_plot_mean{3} = types5_mean;
bar_plot_test = cell(3,1);
bar_plot_test{1} = types2_reorg;
bar_plot_test{2} = types4_reorg;
bar_plot_test{3} = types5_reorg;
hold on;
[position_cell, pvalueArray] = plot_significance_score(bar_plot_mean, bar_plot_test, b);
sigstar(position_cell,pvalueArray)
hold off;
legend({'L2/3', 'L4', 'L5'}, 'Location', 'best');  % Legend

%%  Ratio of each type across layers grouped by neuron
types5 = neuron_total_type_number5./sum(neuron_total_type_number5, 2);

std5 = std(types5, 1)./sqrt(size(types5,1));
types5_mean = mean(types5,1);


types4 = neuron_total_type_number4./sum(neuron_total_type_number4, 2);
std4 = std(types4, 1)./sqrt(size(types4,1));
types4_mean = mean(types4,1);

types2 = neuron_total_type_number2./sum(neuron_total_type_number2, 2);
std2 = std(types2, 1)./sqrt(size(types2,1));
types2_mean = mean(types2,1);

filo = [types2_mean(1), types4_mean(1),types5_mean(1)];
longthin = [types2_mean(3), types4_mean(3),types5_mean(3)];
thin = [types2_mean(4), types4_mean(4),types5_mean(4)];
stubby = [types2_mean(5), types4_mean(5),types5_mean(5)];
mush = [types2_mean(2), types4_mean(2),types5_mean(2)];
branched = [types2_mean(6), types4_mean(6),types5_mean(6)];
data = [filo;longthin;thin;stubby;mush;branched];

filo_err = [std2(1), std4(1),std5(1)];
longthin_err = [std2(3), std4(3),std5(3)];
thin_err = [std2(4), std4(4),std5(4)];
stubby_err = [std2(5), std4(5),std5(5)];
mush_err = [std2(2), std4(2),std5(2)];
branched_err = [std2(6), std4(6),std5(6)];
err = [filo_err;longthin_err;thin_err;stubby_err;mush_err;branched_err];
figure;  % Creates a new figure
b = bar(data,'BarWidth', 0.8);
b(1).FaceColor = [0.8500, 0.3250, 0.0980];  % Color for Variable 1
b(2).FaceColor = [0, 0.4470, 0.7410];  % Color for Variable 2
b(3).FaceColor = [0.9290, 0.6940, 0.1250];  % Color for Variable 3% Plots grouped bar chart
ylabel('Probability of each types of dendrite spine');  % Y-axis label for the values of the variables
title('Distribution of spine types in each neuron');  % Chart title


groupNames = {'Filopodia', 'Long Thin', 'Thin', 'Stubby', 'Mushroom', 'Branched'}; % Custom names for each group of bars
xticks(1:length(groupNames)); % Set the x-axis ticks to match the number of groups
xticklabels(groupNames);

hold on;
numGroups = size(data, 1);
numBars = size(data, 2);

% Add error bars
for i = 1:numBars
    % X positions for error bars
    x = b(i).XEndPoints;
    % Add error bars to each group
    errorbar(x, data(:, i), err(:, i), '.k', 'LineWidth', 1.5);
end
hold on;

types5_reorg = [types5(:,1),types5(:,3), types5(:,4), types5(:,5), types5(:,2), types5(:,6)];
types4_reorg = [types4(:,1),types4(:,3), types4(:,4), types4(:,5), types4(:,2), types4(:,6)];
types2_reorg = [types2(:,1),types2(:,3), types2(:,4), types2(:,5), types2(:,2), types2(:,6)];
bar_plot_mean = cell(3,1);
bar_plot_mean{1} = types2_mean;
bar_plot_mean{2} = types4_mean;
bar_plot_mean{3} = types5_mean;
bar_plot_test = cell(3,1);
bar_plot_test{1} = types2_reorg;
bar_plot_test{2} = types4_reorg;
bar_plot_test{3} = types5_reorg;
hold on;
[position_cell, pvalueArray] = plot_significance_score(bar_plot_mean, bar_plot_test, b);
sigstar(position_cell,pvalueArray)
% x1 = b(1).XEndPoints(6); % X position of bar 1 in group 1
% x2 = b(2).XEndPoints(1); % X position of bar 2 in group 1
% x3 = b(3).XEndPoints(1); % X position of bar 2 in group 1
% yMax = max(data(1, 1:3), [], 'all') + max(err(:, 1:2), [], 'all'); % Maximum y-value plus some margin
% % Plot a line between the bars
% plot([x1, x2], [yMax, yMax], '-k', 'LineWidth', 1.5);
% 
% % Add an asterisk above the line
% text(mean([x1, x2]), yMax + 0.005, '*', 'HorizontalAlignment', 'center', 'FontSize', 14);


hold off;
legend({'L2/3', 'L4', 'L5'}, 'Location', 'best');  % Legend



%% check the correlation between the spine types and dendrite length
reset_values_summary_statistics
id5 = find(dendrite_length5 > 0 & dendrite_length5 < quantile(dendrite_length5, 0.99) & sum(dendrite_each_spine_type_number5,2) > 3);
dendrite_length5 = dendrite_length5(id5);
dendrite_each_spine_type_number5 = dendrite_each_spine_type_number5(id5,:);
density_5 = sum(dendrite_each_spine_type_number5, 2)./ dendrite_length5.*1000;
types5 = dendrite_each_spine_type_number5./dendrite_length5.*1000;
dendrite_radius5 = dendrite_radius5(id5);

id4 = find(dendrite_length4 > 0 & dendrite_length4 < quantile(dendrite_length4, 0.99) & sum(dendrite_each_spine_type_number4,2) > 3);
dendrite_length4 = dendrite_length4(id4);
dendrite_each_spine_type_number4 = dendrite_each_spine_type_number4(id4,:);
density_4 = sum(dendrite_each_spine_type_number4, 2)./ dendrite_length4.*1000;
types4 = dendrite_each_spine_type_number4./dendrite_length4.*1000;
dendrite_radius4 = dendrite_radius4(id4);

id2 = find(dendrite_length2 > 0 & dendrite_length2 < quantile(dendrite_length2, 0.99) & sum(dendrite_each_spine_type_number2,2) > 3);
dendrite_length2 = dendrite_length2(id2);
dendrite_each_spine_type_number2 = dendrite_each_spine_type_number2(id2,:);
density_2 = sum(dendrite_each_spine_type_number2, 2)./ dendrite_length2.*1000;
types2 = dendrite_each_spine_type_number2./dendrite_length2.*1000;
dendrite_radius2 = dendrite_radius2(id2);

% density distribution 

data = [mean(density_2), mean(density_4), mean(density_5)];
err = [[std(density_2)/sqrt(length(density_2)), std(density_4)/sqrt(length(density_4)), std(density_5)/sqrt(length(density_5))]];
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
% ylabel('radius (nm)');  % Y-axis label for the values of the variables
title('spine density across different layers');  % Chart title
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
bar_plot_test{1,1} = density_2;
bar_plot_test{2,1} = density_4;
bar_plot_test{3,1} = density_5;
hold on;
[position_cell, pvalueArray] = plot_significance_score_cell_cell(bar_plot_mean, bar_plot_test, b);
sigstar(position_cell,pvalueArray)


% dendrite radius across different layers

data = [mean(dendrite_radius2), mean(dendrite_radius4), mean(dendrite_radius5)];
err = [[std(dendrite_radius2)/sqrt(length(dendrite_radius2)), std(dendrite_radius4)/sqrt(length(dendrite_radius4)), std(dendrite_radius5)/sqrt(length(dendrite_radius5))]];
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
% ylabel('radius (nm)');  % Y-axis label for the values of the variables
title('dendrite radius across different layers');  % Chart title
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
bar_plot_test{1,1} = dendrite_radius2;
bar_plot_test{2,1} = dendrite_radius4;
bar_plot_test{3,1} = dendrite_radius5;
hold on;
[position_cell, pvalueArray] = plot_significance_score_cell_cell(bar_plot_mean, bar_plot_test, b);
sigstar(position_cell,pvalueArray)





% check by constraining the density to be within the same range
lowerThres = quantile([density_2(:);density_4(:); density_5(:)], 0.5);
upperThres = quantile([density_2(:);density_4(:); density_5(:)], 0.99);
types5 = types5((density_5 < upperThres) & (density_5 > lowerThres) ,:);
types2 = types2((density_2 < upperThres) & (density_2 > lowerThres),:);
types4 = types4((density_4 < upperThres) & (density_4 > lowerThres),:);


types5_reorg = [types5(:,1),types5(:,3), types5(:,4), types5(:,2)];
types4_reorg = [types4(:,1),types4(:,3), types4(:,4), types4(:,2)];
types2_reorg = [types2(:,1),types2(:,3), types2(:,4), types2(:,2)]; % stubby and branched are removed
types5_reorg = types5_reorg./max(types5_reorg, [], 2).*[1/4,2/4,3/4,4/4];
types5_sum = sum(types5_reorg, 2);
types4_reorg = types4_reorg./max(types4_reorg, [], 2).*[1/4,2/4,3/4,4/4];
types4_sum = sum(types4_reorg, 2);
types2_reorg = types2_reorg./max(types2_reorg, [], 2).*[1/4,2/4,3/4,4/4];
types2_sum = sum(types2_reorg, 2);



% x = [types2_sum(:);types4_sum(:);types5_sum(:)];
% y = [zeros(length(types2_sum(:)),1);ones(length(types4_sum(:)),1);ones(length(types5_sum(:)),1).*2];
% 
% figure; violinplot(x,y,'ShowMean',true,'ShowMedian',true); title(['spine type score different layers'])
% set(gca,'xtick',[1,2,3],'xticklabel',{'L2/3', 'L4', 'L5'})
% [h,p1,ci,stats] = ttest2(types2_sum,types4_sum);
% [h,p2,ci,stats] = ttest2(types4_sum,types5_sum);
% [h,p3,ci,stats] = ttest2(types2_sum,types5_sum);
% sigstar({[1,2],[2,3], [1,3]},[p1,p2,p3])
% ylim([0.2,1])

%bar plot
data = [mean(types2_sum), mean(types4_sum), mean(types5_sum)];
err = [[std(types2_sum)/sqrt(length(types2_sum)), std(types4_sum)/sqrt(length(types4_sum)), std(types5_sum)/sqrt(length(types5_sum))]];
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
% ylabel('radius (nm)');  % Y-axis label for the values of the variables
title('spine type score different layers');  % Chart title
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
bar_plot_test{1,1} = types2_sum;
bar_plot_test{2,1} = types4_sum;
bar_plot_test{3,1} = types5_sum;
hold on;
[position_cell, pvalueArray] = plot_significance_score_cell_cell(bar_plot_mean, bar_plot_test, b);
sigstar(position_cell,pvalueArray)


% check the correlation between density and spine types



id22 = find(density_2 < quantile(density_2, 0.99));
id42 = find(density_4 < quantile(density_4, 0.99));
id52 = find(density_5 < quantile(density_5, 0.99));
x = [density_2(id22);density_4(id42);density_5(id52)];
y = [zeros(length(density_2(id22)),1);ones(length(density_4(id42)),1);ones(length(density_5(id52)),1).*2];

figure; violinplot(x,y,'ShowMean',true,'ShowMedian',true); title(['spine density across different layers'])
set(gca,'xtick',[1,2,3],'xticklabel',{'L2/3', 'L4', 'L5'})
[h,p1,ci,stats] = ttest2(density_2(id22),density_4(id42));
[h,p2,ci,stats] = ttest2(density_4(id42),density_5(id52));
[h,p3,ci,stats] = ttest2(density_2(id22),density_5(id52));
sigstar({[1,2],[2,3], [1,3]},[p1,p2,p3])

ylabel('spine density per \mu m');

plot_regression_L2L4L5(types2_sum(id22), types4_sum(id42), types5_sum(id52), density_2(id22), density_4(id42), density_5(id52))
% check the correlation between dendrite radius and the score type

id22 = find(density_2 < quantile(density_2, 0.99));
id42 = find(density_4 < quantile(density_4, 0.99));
id52 = find(density_5 < quantile(density_5, 0.99));
plot_regression_L2L4L5(types2_sum(id22), types4_sum(id42), types5_sum(id52), ...
    dendrite_radius2(id22), dendrite_radius4(id42), dendrite_radius5(id52))

plot_kde_contour_3(types2_sum(id22),dendrite_radius2(id22),types4_sum(id42),dendrite_radius4(id42),types5_sum(id52),dendrite_radius5(id52));



% check the distribution of the joint probability of dendrite radius and
% the spine density of each corresponding dendrite. The dendrites with very
% large density could be the reason of the small length in its fragments
% which are removed.
id22 = find(density_2 < quantile(density_2, 0.6));
id42 = find(density_4 < quantile(density_4, 0.6));
id52 = find(density_5 < quantile(density_5, 0.6));
plot_kde_contour_3(density_2(id22),dendrite_radius2(id22), density_4(id42),dendrite_radius4(id42),density_5(id52), dendrite_radius5(id52))
% Based on the distribution of density, length, radius, remove the outliers
figure; scatter3(dendrite_radius2,dendrite_length2,density_2)
hold on; scatter3(dendrite_radius4,dendrite_length4,density_4)
hold on; scatter3(dendrite_radius5,dendrite_length5,density_5)
xlabel('radius'); ylabel('length');zlabel('density')
legend('L2/3', 'L4','L5')




% !!check the distribution of the joint probability of dendrite radius and select the points that are close to the center (remove the top 10% outliers) 
reset_values_summary_statistics
id5 = find(dendrite_length5 > 0 & dendrite_length5 < quantile(dendrite_length5, 0.99) & sum(dendrite_each_spine_type_number5,2) > 3);
dendrite_length5 = dendrite_length5(id5);
dendrite_each_spine_type_number5 = dendrite_each_spine_type_number5(id5,:);
density_5 = sum(dendrite_each_spine_type_number5, 2)./ dendrite_length5.*1000;
types5 = dendrite_each_spine_type_number5./dendrite_length5.*1000;
dendrite_radius5 = dendrite_radius5(id5);

id4 = find(dendrite_length4 > 0 & dendrite_length4 < quantile(dendrite_length4, 0.99) & sum(dendrite_each_spine_type_number4,2) > 3);
dendrite_length4 = dendrite_length4(id4);
dendrite_each_spine_type_number4 = dendrite_each_spine_type_number4(id4,:);
density_4 = sum(dendrite_each_spine_type_number4, 2)./ dendrite_length4.*1000;
types4 = dendrite_each_spine_type_number4./dendrite_length4.*1000;
dendrite_radius4 = dendrite_radius4(id4);

id2 = find(dendrite_length2 > 0 & dendrite_length2 < quantile(dendrite_length2, 0.99) & sum(dendrite_each_spine_type_number2,2) > 3);
dendrite_length2 = dendrite_length2(id2);
dendrite_each_spine_type_number2 = dendrite_each_spine_type_number2(id2,:);
density_2 = sum(dendrite_each_spine_type_number2, 2)./ dendrite_length2.*1000;
types2 = dendrite_each_spine_type_number2./dendrite_length2.*1000;
dendrite_radius2 = dendrite_radius2(id2);

% figure; scatter3([dendrite_radius2(:); dendrite_radius4(:); dendrite_radius5(:)],[dendrite_length2(:); dendrite_length4(:); dendrite_length5(:)],[density_2(:);density_4(:);density_5(:)])

types5_reorg = [types5(:,1),types5(:,3), types5(:,4), types5(:,2)];
types4_reorg = [types4(:,1),types4(:,3), types4(:,4), types4(:,2)];
types2_reorg = [types2(:,1),types2(:,3), types2(:,4), types2(:,2)]; % stubby and branched are removed
% types5_reorg = types5_reorg./max(types5_reorg, [], 2).*[1/4,2/4,3/4,4/4]; % normalized so that the maximum is 1
types5_reorg = types5_reorg./sum(types5_reorg, 2).*[1/4,2/4,3/4,4/4]; % normalized so that the maximum is 1
types5_sum = sum(types5_reorg, 2);
% types4_reorg = types4_reorg./max(types4_reorg, [], 2).*[1/4,2/4,3/4,4/4];
types4_reorg = types4_reorg./sum(types4_reorg, 2).*[1/4,2/4,3/4,4/4]; % normalized so that the maximum is 1
types4_sum = sum(types4_reorg, 2);
% types2_reorg = types2_reorg./max(types2_reorg, [], 2).*[1/4,2/4,3/4,4/4];
types2_reorg = types2_reorg./sum(types2_reorg, 2).*[1/4,2/4,3/4,4/4];
types2_sum = sum(types2_reorg, 2);

total_radius = [dendrite_radius2(:); dendrite_radius4(:); dendrite_radius5(:)];
max_radius = max(total_radius);
total_radius = total_radius./max(total_radius);
total_length = [dendrite_length2(:); dendrite_length4(:); dendrite_length5(:)];
max_length = max(total_length);
total_length = total_length./max(total_length);
total_density = [density_2(:);density_4(:);density_5(:)];
max_density = max(total_density);
total_density = total_density./max(total_density);
center = [mean(total_radius), mean(total_length),mean(total_density)];
dist_all = vecnorm([total_radius - center(1), total_length - center(2), total_density - center(3)], 2, 2);
threshold_ratio = 0.9;
selected_radius = total_radius(dist_all < quantile(dist_all , threshold_ratio)).*max_radius;
selected_length = total_length(dist_all < quantile(dist_all , threshold_ratio)).*max_length;
selected_density = total_density(dist_all < quantile(dist_all , threshold_ratio)).*max_density;
% hold on; scatter3(selected_radius,selected_length,selected_density)
id2 = find(dendrite_radius2 >= min(selected_radius) & dendrite_length2 >= min(selected_length) & density_2 >= min(selected_density)...
    & dendrite_radius2 <= max(selected_radius) & dendrite_length2 <= max(selected_length) & density_2 <= max(selected_density));
id4 = find(dendrite_radius4 >= min(selected_radius) & dendrite_length4 >= min(selected_length) & density_4 >= min(selected_density)...
    & dendrite_radius4 <= max(selected_radius) & dendrite_length4 <= max(selected_length) & density_4 <= max(selected_density));
id5 = find(dendrite_radius5 >= min(selected_radius) & dendrite_length5 >= min(selected_length) & density_5 >= min(selected_density)...
    & dendrite_radius5 <= max(selected_radius) & dendrite_length5 <= max(selected_length) & density_5 <= max(selected_density));
plot_kde_contour_3(types2_sum(id2),dendrite_radius2(id2),types4_sum(id4),dendrite_radius4(id4),types5_sum(id5),dendrite_radius5(id5));
legend('L2/3', 'L4', 'L5')
xlabel('spine type score'); ylabel('dendrite radius(nm)');title('Joint distribution of spine score and radius')

plot_kde_contour_3(density_2(id2),dendrite_radius2(id2),density_4(id4),dendrite_radius4(id4),density_5(id5),dendrite_radius5(id5));
legend('L2/3', 'L4', 'L5')
xlabel('spine density'); ylabel('dendrite radius(nm)');title('Joint distribution of spine density and radius')
xlim([0, 20])







%% spine type vs dendrite radius


% !!check the distribution of the joint probability of dendrite radius and select the points that are close to the center (remove the top 10% outliers) 
reset_values_summary_statistics
id5 = find(dendrite_length5 > 0 & dendrite_length5 < quantile(dendrite_length5, 0.99) & sum(dendrite_each_spine_type_number5,2) > 3);
dendrite_length5 = dendrite_length5(id5);
dendrite_each_spine_type_number5 = dendrite_each_spine_type_number5(id5,:);
density_5 = sum(dendrite_each_spine_type_number5, 2)./ dendrite_length5.*1000;
types5 = dendrite_each_spine_type_number5./dendrite_length5.*1000;
dendrite_radius5 = dendrite_radius5(id5);

id4 = find(dendrite_length4 > 0 & dendrite_length4 < quantile(dendrite_length4, 0.99) & sum(dendrite_each_spine_type_number4,2) > 3);
dendrite_length4 = dendrite_length4(id4);
dendrite_each_spine_type_number4 = dendrite_each_spine_type_number4(id4,:);
density_4 = sum(dendrite_each_spine_type_number4, 2)./ dendrite_length4.*1000;
types4 = dendrite_each_spine_type_number4./dendrite_length4.*1000;
dendrite_radius4 = dendrite_radius4(id4);

id2 = find(dendrite_length2 > 0 & dendrite_length2 < quantile(dendrite_length2, 0.99) & sum(dendrite_each_spine_type_number2,2) > 3);
dendrite_length2 = dendrite_length2(id2);
dendrite_each_spine_type_number2 = dendrite_each_spine_type_number2(id2,:);
density_2 = sum(dendrite_each_spine_type_number2, 2)./ dendrite_length2.*1000;
types2 = dendrite_each_spine_type_number2./dendrite_length2.*1000;
dendrite_radius2 = dendrite_radius2(id2);

% figure; scatter3([dendrite_radius2(:); dendrite_radius4(:); dendrite_radius5(:)],[dendrite_length2(:); dendrite_length4(:); dendrite_length5(:)],[density_2(:);density_4(:);density_5(:)])

types5_reorg = [types5(:,1),types5(:,3), types5(:,4), types5(:,2)];
types4_reorg = [types4(:,1),types4(:,3), types4(:,4), types4(:,2)];
types2_reorg = [types2(:,1),types2(:,3), types2(:,4), types2(:,2)]; % stubby and branched are removed
% types5_reorg = types5_reorg./max(types5_reorg, [], 2).*[1/4,2/4,3/4,4/4]; % normalized so that the maximum is 1
types5_reorg = types5_reorg./sum(types5_reorg, 2).*[1/4,2/4,3/4,4/4]; % normalized so that the maximum is 1
types5_sum = sum(types5_reorg, 2);
% types4_reorg = types4_reorg./max(types4_reorg, [], 2).*[1/4,2/4,3/4,4/4];
types4_reorg = types4_reorg./sum(types4_reorg, 2).*[1/4,2/4,3/4,4/4]; % normalized so that the maximum is 1
types4_sum = sum(types4_reorg, 2);
% types2_reorg = types2_reorg./max(types2_reorg, [], 2).*[1/4,2/4,3/4,4/4];
types2_reorg = types2_reorg./sum(types2_reorg, 2).*[1/4,2/4,3/4,4/4];
types2_sum = sum(types2_reorg, 2);

total_radius = [dendrite_radius2(:); dendrite_radius4(:); dendrite_radius5(:)];
max_radius = max(total_radius);
total_radius = total_radius./max(total_radius);
total_length = [dendrite_length2(:); dendrite_length4(:); dendrite_length5(:)];
max_length = max(total_length);
total_length = total_length./max(total_length);
total_density = [density_2(:);density_4(:);density_5(:)];
max_density = max(total_density);
total_density = total_density./max(total_density);
center = [mean(total_radius), mean(total_length),mean(total_density)];
dist_all = vecnorm([total_radius - center(1), total_length - center(2), total_density - center(3)], 2, 2);
threshold_ratio = 0.9;
selected_radius = total_radius(dist_all < quantile(dist_all , threshold_ratio)).*max_radius;
selected_length = total_length(dist_all < quantile(dist_all , threshold_ratio)).*max_length;
selected_density = total_density(dist_all < quantile(dist_all , threshold_ratio)).*max_density;
% hold on; scatter3(selected_radius,selected_length,selected_density)
id22 = find(dendrite_radius2 >= min(selected_radius) & dendrite_length2 >= min(selected_length) & density_2 >= min(selected_density)...
    & dendrite_radius2 <= max(selected_radius) & dendrite_length2 <= max(selected_length) & density_2 <= max(selected_density));
id42 = find(dendrite_radius4 >= min(selected_radius) & dendrite_length4 >= min(selected_length) & density_4 >= min(selected_density)...
    & dendrite_radius4 <= max(selected_radius) & dendrite_length4 <= max(selected_length) & density_4 <= max(selected_density));
id52 = find(dendrite_radius5 >= min(selected_radius) & dendrite_length5 >= min(selected_length) & density_5 >= min(selected_density)...
    & dendrite_radius5 <= max(selected_radius) & dendrite_length5 <= max(selected_length) & density_5 <= max(selected_density));
plot_kde_contour_3(types2_sum(id22),dendrite_radius2(id22),types4_sum(id42),dendrite_radius4(id42),types5_sum(id52),dendrite_radius5(id52));
legend('L2/3', 'L4', 'L5')
xlabel('spine type score'); ylabel('dendrite radius(nm)');title('Joint distribution of spine score and radius')

plot_kde_contour_3(density_2(id22),dendrite_radius2(id22),density_4(id42),dendrite_radius4(id42),density_5(id52),dendrite_radius5(id52));
legend('L2/3', 'L4', 'L5')
xlabel('spine density'); ylabel('dendrite radius(nm)');title('Joint distribution of spine density and radius')
xlim([0, 20])



% old


reset_values_summary_statistics
id5 = find(dendrite_radius5 >0 & sum(dendrite_each_spine_type_number5,2) > 3);
dendrite_radius5 = dendrite_radius5(id5);
dendrite_each_spine_type_number5 = dendrite_each_spine_type_number5(id5,:);

id4 = find(dendrite_radius4 > 0 & sum(dendrite_each_spine_type_number4,2) > 3);
dendrite_radius4 = dendrite_radius4(id4);
dendrite_each_spine_type_number4 = dendrite_each_spine_type_number4(id4,:);


id2 = find(dendrite_radius2 > 0 & sum(dendrite_each_spine_type_number2,2) > 3);
dendrite_radius2 = dendrite_radius2(id2);
dendrite_each_spine_type_number2 = dendrite_each_spine_type_number2(id2,:);

types5 = dendrite_each_spine_type_number5;
types4 = dendrite_each_spine_type_number4;
types2 = dendrite_each_spine_type_number2;

types5_reorg = [types5(:,1),types5(:,3), types5(:,4), types5(:,5), types5(:,2)];
types4_reorg = [types4(:,1),types4(:,3), types4(:,4), types4(:,5), types4(:,2)];
types2_reorg = [types2(:,1),types2(:,3), types2(:,4), types2(:,5), types2(:,2)];
types5_reorg = types5_reorg./sum(types5_reorg,2).*[1/5,2/5,3/5,4/5,1];
types5_sum = sum(types5_reorg, 2);
types4_reorg = types4_reorg./sum(types4_reorg,2).*[1/5,2/5,3/5,4/5,1];
types4_sum = sum(types4_reorg, 2);
types2_reorg = types2_reorg./sum(types2_reorg,2).*[1/5,2/5,3/5,4/5,1];
types2_sum = sum(types2_reorg, 2);
% 
% id22 = find(density_2 < quantile(density_2, 0.95));
% id42 = find(density_4 < quantile(density_4, 0.95));
% id52 = find(density_5 < quantile(density_5, 0.95));
plot_regression_L2L4L5(types2_sum, types4_sum, types5_sum, dendrite_radius2, dendrite_radius4, dendrite_radius5)
plot_regression_L2L4L5(dendrite_radius2, dendrite_radius4, dendrite_radius5,density_2, density_4, density_5)


x = [dendrite_radius2(:);dendrite_radius4(:);dendrite_radius5(:)];
y = [zeros(length(dendrite_radius2(:)),1);ones(length(dendrite_radius4(:)),1);ones(length(dendrite_radius5(:)),1).*2];

figure; violinplot(x,y,'ShowMean',true,'ShowMedian',true); title(['dendrite radius'])
set(gca,'xtick',[1,2,3],'xticklabel',{'L2/3', 'L4', 'L5'}); ylabel('radius (nm)')
[h,p1,ci,stats] = ttest2(dendrite_radius2,dendrite_radius4);
[h,p2,ci,stats] = ttest2(dendrite_radius4,dendrite_radius5);
[h,p3,ci,stats] = ttest2(dendrite_radius2,dendrite_radius5);
sigstar({[1,2],[2,3], [1,3]},[p1,p2,p3])


%% whether the score correlates with the astrocyte wrapping ratio
% 1: filopodia, 2: mushroom, 3: long thin, 4: thin, 5: stubby
reset_values_summary_statistics;            
spine_dendrite_idx2 = label2idx(spine_dendrite_type_label2(:,2));
spine_dendrite_idx2(cellfun(@length, spine_dendrite_idx2) <= 3) = [];
% spine_score_2 = cellfun(@(c) (sum(spine_dendrite_type_label2(c,1)== 1)*(1/5) + sum(spine_dendrite_type_label2(c,1)== 3)*(2/5) ...
% + sum(spine_dendrite_type_label2(c,1)== 4)*(3/5) + sum(spine_dendrite_type_label2(c,1)== 5)*(4/5) + sum(spine_dendrite_type_label2(c,1)== 2))/(length(c)), spine_dendrite_idx2);
spine_score_2 = cellfun(@(c) (sum(spine_dendrite_type_label2(c,1)== 1)*(1/4) + sum(spine_dendrite_type_label2(c,1)== 3)*(2/4) ...
 + sum(spine_dendrite_type_label2(c,1)== 4)*(3/4) + sum(spine_dendrite_type_label2(c,1)== 2))/(length(c)), spine_dendrite_idx2);
perimeter_ratio_2 = cellfun(@(c) mean(sinsperimeterRatio2(c), "omitnan"), spine_dendrite_idx2);
preSynapse_ratio_2 = cellfun(@(c) mean(sinspreSynapseTouchingRatio2(c), "omitnan"), spine_dendrite_idx2);
postSynapse_ratio_2 = cellfun(@(c) mean(sinspostSynapseTouchingRatio2(c), "omitnan"), spine_dendrite_idx2);
% headneck_ratio_2 = cellfun(@(c) mean(singleSynHeadNeckTouchingRatio2(c,1)/singleSynHeadNeckTouchingRatio2(c,2),"omitnan"), spine_dendrite_idx2);

spine_dendrite_idx4 = label2idx(spine_dendrite_type_label4(:,2));
spine_dendrite_idx4(cellfun(@length, spine_dendrite_idx4) <= 3) = [];
spine_score_4 = cellfun(@(c) (sum(spine_dendrite_type_label4(c,1)== 1)*(1/5) + sum(spine_dendrite_type_label4(c,1)== 3)*(2/5) ...
    + sum(spine_dendrite_type_label4(c,1)== 4)*(3/5) + sum(spine_dendrite_type_label4(c,1)== 5)*(4/5) + sum(spine_dendrite_type_label4(c,1)== 2))/(length(c)), spine_dendrite_idx4);
perimeter_ratio_4 = cellfun(@(c) mean(sinsperimeterRatio4(c),"omitnan"), spine_dendrite_idx4);
preSynapse_ratio_4 = cellfun(@(c) mean(sinspreSynapseTouchingRatio4(c),"omitnan"), spine_dendrite_idx4);
postSynapse_ratio_4 = cellfun(@(c) mean(sinspostSynapseTouchingRatio4(c),"omitnan"), spine_dendrite_idx4);
% headneck_ratio_4 = cellfun(@(c) sum(singleSynHeadNeckTouchingRatio4(c)), spine_dendrite_idx4);

spine_dendrite_idx5 = label2idx(spine_dendrite_type_label5(:,2));
spine_dendrite_idx5(cellfun(@length, spine_dendrite_idx5) <= 3) = [];
spine_score_5 = cellfun(@(c) (sum(spine_dendrite_type_label5(c,1)== 1)*(1/5) + sum(spine_dendrite_type_label5(c,1)== 3)*(2/5) ...
    + sum(spine_dendrite_type_label5(c,1)== 4)*(3/5) + sum(spine_dendrite_type_label5(c,1)== 5)*(4/5) + sum(spine_dendrite_type_label5(c,1)== 2))/(length(c)), spine_dendrite_idx5);
perimeter_ratio_5 = cellfun(@(c) mean(sinsperimeterRatio5(c),"omitnan"), spine_dendrite_idx5);
preSynapse_ratio_5 = cellfun(@(c) mean(sinspreSynapseTouchingRatio5(c),"omitnan"), spine_dendrite_idx5);
postSynapse_ratio_5 = cellfun(@(c) mean(sinspostSynapseTouchingRatio5(c),"omitnan"), spine_dendrite_idx5);
% headneck_ratio_5 = cellfun(@(c) sum(singleSynHeadNeckTouchingRatio5(c)), spine_dendrite_idx5);

spine_score_2_range = linspace(min(spine_score_2), max(spine_score_2), 100);
perimeter_ratio_2_range = linspace(min(perimeter_ratio_2), max(perimeter_ratio_2), 100);
[X2, Y2] = meshgrid(spine_score_2_range, perimeter_ratio_2_range);    
density2 = ksdensity([spine_score_2(:), perimeter_ratio_2(:)], [X2(:), Y2(:)]);
density2 = reshape(density2, size(X2)); 

spine_score_4_range = linspace(min(spine_score_4), max(spine_score_4), 100);
perimeter_ratio_4_range = linspace(min(perimeter_ratio_4), max(perimeter_ratio_4), 100);
[X4, Y4] = meshgrid(spine_score_4_range, perimeter_ratio_4_range);
density4 = ksdensity([spine_score_4(:), perimeter_ratio_4(:)], [X4(:), Y4(:)]);
density4 = reshape(density4, size(X4));

spine_score_5_range = linspace(min(spine_score_5), max(spine_score_5), 100);
perimeter_ratio_5_range = linspace(min(perimeter_ratio_5), max(perimeter_ratio_5), 100);
[X5, Y5] = meshgrid(spine_score_5_range, perimeter_ratio_5_range);
density5 = ksdensity([spine_score_5(:), perimeter_ratio_5(:)], [X5(:), Y5(:)]);
density5 = reshape(density5, size(X5));


figure;
contour(X2, Y2, density2, 'LineWidth', 2, 'LineColor', 'r');
hold on;
contour(X4, Y4, density4, 'LineWidth', 2,'LineColor', 'g');
hold on;
contour(X5, Y5, density5, 'LineWidth', 2, 'LineColor', 'b');

ylim([0,0.6]);
xlabel('spine score')
ylabel('perimeter wrapping ratio')
title('perimeter wrapping ratio vs spine score')

legend({'L2/3', 'L4', 'L5'}, 'Location', 'best')




spine_score_2_range = linspace(min(spine_score_2), max(spine_score_2), 100);
preSynapse_ratio_2_range = linspace(min(preSynapse_ratio_2), max(preSynapse_ratio_2), 100);
[X2, Y2] = meshgrid(spine_score_2_range, preSynapse_ratio_2_range);    
density2 = ksdensity([spine_score_2(:), preSynapse_ratio_2(:)], [X2(:), Y2(:)]);
density2 = reshape(density2, size(X2)); 

spine_score_4_range = linspace(min(spine_score_4), max(spine_score_4), 100);
preSynapse_ratio_4_range = linspace(min(preSynapse_ratio_4), max(preSynapse_ratio_4), 100);
[X4, Y4] = meshgrid(spine_score_4_range, preSynapse_ratio_4_range);
density4 = ksdensity([spine_score_4(:), preSynapse_ratio_4(:)], [X4(:), Y4(:)]);
density4 = reshape(density4, size(X4));

spine_score_5_range = linspace(min(spine_score_5), max(spine_score_5), 100);
preSynapse_ratio_5_range = linspace(min(preSynapse_ratio_5), max(preSynapse_ratio_5), 100);
[X5, Y5] = meshgrid(spine_score_5_range, preSynapse_ratio_5_range);
density5 = ksdensity([spine_score_5(:), preSynapse_ratio_5(:)], [X5(:), Y5(:)]);
density5 = reshape(density5, size(X5));


figure;
contour(X2, Y2, density2, 'LineWidth', 2, 'LineColor', 'r');
hold on;
contour(X4, Y4, density4, 'LineWidth', 2,'LineColor', 'g');
hold on;
contour(X5, Y5, density5, 'LineWidth', 2, 'LineColor', 'b');

ylim([0,0.6]);
xlabel('spine score')
ylabel('presynapse wrapping ratio')
title('presynapse wrapping ratio vs spine score')

legend({'L2/3', 'L4', 'L5'}, 'Location', 'best')




spine_score_2_range = linspace(min(spine_score_2), max(spine_score_2), 100);
postSynapse_ratio_2_range = linspace(min(postSynapse_ratio_2), max(postSynapse_ratio_2), 100);
[X2, Y2] = meshgrid(spine_score_2_range, postSynapse_ratio_2_range);    
density2 = ksdensity([spine_score_2(:), postSynapse_ratio_2(:)], [X2(:), Y2(:)]);
density2 = reshape(density2, size(X2)); 

spine_score_4_range = linspace(min(spine_score_4), max(spine_score_4), 100);
postSynapse_ratio_4_range = linspace(min(postSynapse_ratio_4), max(postSynapse_ratio_4), 100);
[X4, Y4] = meshgrid(spine_score_4_range, postSynapse_ratio_4_range);
density4 = ksdensity([spine_score_4(:), postSynapse_ratio_4(:)], [X4(:), Y4(:)]);
density4 = reshape(density4, size(X4));

spine_score_5_range = linspace(min(spine_score_5), max(spine_score_5), 100);
postSynapse_ratio_5_range = linspace(min(postSynapse_ratio_5), max(postSynapse_ratio_5), 100);
[X5, Y5] = meshgrid(spine_score_5_range, postSynapse_ratio_5_range);
density5 = ksdensity([spine_score_5(:), postSynapse_ratio_5(:)], [X5(:), Y5(:)]);
density5 = reshape(density5, size(X5));


figure;
contour(X2, Y2, density2, 'LineWidth', 2, 'LineColor', 'r');
hold on;
contour(X4, Y4, density4, 'LineWidth', 2,'LineColor', 'g');
hold on;
contour(X5, Y5, density5, 'LineWidth', 2, 'LineColor', 'b');

ylim([0,0.6]);
xlabel('spine score')
ylabel('postsynapse wrapping ratio')
title('postsynapse wrapping ratio vs spine score')

legend({'L2/3', 'L4', 'L5'}, 'Location', 'best')



%% within layer comparison
