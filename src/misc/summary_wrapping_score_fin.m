clear
singleSynHeadVolume = [];
singleSynHeadMeanRadius = [];
singleSynNeckLength = [];
singleSynNeckSection = [];
singleSynNeckMeanRadius = [];
singleDendriteDensity = [];
doubleSynHeadVolume = [];
doubleSynMeanHeadRadius = [];
doubleSynNeckLength = [];
doubleSynNeckSection = [];
doubleSynNeckMeanRadius = [];
doubleDendriteDensity = [];
singleSynapticCleftSize = [];
sinsperimeterRatio = [];
sinsperimeterWeightedWrappingArea = [];
sinspostSynapseTouchingArea = [];
sinspostSynapseTouchingRatio = [];
sinspreSynapseTouchingArea = [];
sinspreSynapseTouchingRatio = [];
singleSynHeadNeckTouchingArea = [];
singleSynHeadNeckTouchingRatio = [];
doubleSynapticCleftSize = [];
dousperimeterRatio = [];
dousperimeterWeightedWrappingArea = [];
douspostSynapseTouchingArea = [];
douspostSynapseTouchingRatio = [];
douspreSynapseTouchingArea = [];
douspreSynapseTouchingRatio = [];
doubleSynHeadNeckTouchingArea = [];
doubleSynHeadNeckTouchingRatio = [];
doubleSynPosition = [];
spineCoordinate = [];
singleSynapse_astroID_neuronID_folder = [];
%[2,4,5,6,7,8,9,10,12,13,14,15,16,17,18,19]
rootFolder = '/work/boyu/EM_astrocyte/astro_11_33_16_16_40/';
neuronListFolder = '/work/boyu/EM_astrocyte/astro_11_33_64_64_80/';
addpath('/work/boyu/EM_astrocyte/EM_dendrite/resources/sigstar-master')
addpath('../resources/Violinplot-Matlab/')

for nf = [23:28]

    disp(nf)
    fullSegfolder_root = [rootFolder,'astro_', num2str(nf), '_minnie65'];
    neuronList = [neuronListFolder, 'astro_', num2str(nf), '_minnie65/', 'top20_neuron_id_no_soma_filtered.txt'];
    opts = delimitedTextImportOptions("NumVariables", 1);

    % Specify range and delimiter
    opts.DataLines = [1, Inf];
    opts.Delimiter = ",";
    
    % Specify column names and types
    opts.VariableNames = "VarName1";
    opts.VariableTypes = "uint64";
    
    % Specify file level properties
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";    
    neuronList_str = table2array(readtable(neuronList, opts));
   for m = 1:length(neuronList_str)
        % create a folder to save the summary.mat for each dendrite ID
        % nf = 26, m = 5
        folderList = zeros(125,3);
        count = 1;
        for ix = 0:4
            for iy = 0:4
                for iz = 0:4
                    folderList(count,:) = [ix, iy, iz];
                    count = count + 1;
                end
            end
        end
        curpsID = neuronList_str(m);
        matSaveFolder = fullfile(fullSegfolder_root, [num2str(curpsID),'_quantify_score_v4']);
        if(exist(fullfile(matSaveFolder, 'summary.mat'), 'file'))
          load(fullfile(matSaveFolder, 'summary.mat'));
        end

        se = strel('sphere', 2);        
        singleSynHeadVolume = [singleSynHeadVolume;summaryStructure.ss1];
        singleSynHeadMeanRadius = [singleSynHeadMeanRadius;summaryStructure.ss2];
        singleSynNeckLength = [singleSynNeckLength;summaryStructure.ss3];
        singleSynNeckSection = [singleSynNeckSection;summaryStructure.ss4];
        singleSynNeckMeanRadius = [singleSynNeckMeanRadius;summaryStructure.ss5];
        singleSynapticCleftSize = [singleSynapticCleftSize;summaryStructure.sa1];
        sinsperimeterRatio = [sinsperimeterRatio;summaryStructure.sa2];
        sinsperimeterWeightedWrappingArea = [sinsperimeterWeightedWrappingArea;summaryStructure.sa3];
        sinspostSynapseTouchingArea = [sinspostSynapseTouchingArea;summaryStructure.sa4];
        sinspostSynapseTouchingRatio = [sinspostSynapseTouchingRatio;summaryStructure.sa5];
        sinspreSynapseTouchingArea = [sinspreSynapseTouchingArea;summaryStructure.sa6];
        sinspreSynapseTouchingRatio = [sinspreSynapseTouchingRatio;summaryStructure.sa7];
        singleSynHeadNeckTouchingArea = [singleSynHeadNeckTouchingArea;summaryStructure.sa8];
        singleSynHeadNeckTouchingRatio = [singleSynHeadNeckTouchingRatio;summaryStructure.sa9];
    
        doubleSynHeadVolume =[doubleSynHeadVolume;summaryStructure.ds1];
        doubleSynMeanHeadRadius = [doubleSynMeanHeadRadius;summaryStructure.ds2];
        doubleSynNeckLength = [doubleSynNeckLength;summaryStructure.ds3];
        doubleSynNeckSection = [doubleSynNeckSection;summaryStructure.ds4];
        doubleSynNeckMeanRadius = [doubleSynNeckMeanRadius;summaryStructure.ds5];
        doubleSynapticCleftSize = [doubleSynapticCleftSize;summaryStructure.da1];
        dousperimeterRatio= [dousperimeterRatio;summaryStructure.da2];
        dousperimeterWeightedWrappingArea = [dousperimeterWeightedWrappingArea;summaryStructure.da3];
        douspostSynapseTouchingArea = [douspostSynapseTouchingArea;summaryStructure.da4];
        douspostSynapseTouchingRatio = [douspostSynapseTouchingRatio;summaryStructure.da5];
        douspreSynapseTouchingArea = [douspreSynapseTouchingArea;summaryStructure.da6];
        douspreSynapseTouchingRatio = [douspreSynapseTouchingRatio;summaryStructure.da7];
        doubleSynHeadNeckTouchingArea = [doubleSynHeadNeckTouchingArea;summaryStructure.da8];
        doubleSynHeadNeckTouchingRatio = [doubleSynHeadNeckTouchingRatio;summaryStructure.da9];
        
        density_score_arr = 0;
%          density_score = density_score_arr(1)/density_score_arr(2);
%         disp([num2str(density_score),' ', num2str(curpsID)])
        
   end

end
layer2_summary.singleSynHeadVolume = singleSynHeadVolume;
layer2_summary.singleSynMeanHeadRadius = singleSynHeadMeanRadius;
layer2_summary.singleSynNeckLength = singleSynNeckLength;
layer2_summary.singleSynNeckSection = singleSynNeckSection;
layer2_summary.singleSynNeckMeanRadius = singleSynNeckMeanRadius;
layer2_summary.doubleSynHeadVolume = doubleSynHeadVolume;
layer2_summary.doubleSynMeanHeadRadius = doubleSynMeanHeadRadius;
layer2_summary.doubleSynNeckLength = doubleSynNeckLength;
layer2_summary.doubleSynNeckSection = doubleSynNeckSection;
layer2_summary.doubleSynNeckMeanRadius = doubleSynNeckMeanRadius;
layer2_summary.singleSynapticCleftSize = singleSynapticCleftSize;
layer2_summary.sinsperimeterRatio = sinsperimeterRatio;
layer2_summary.sinsperimeterWeightedWrappingArea = sinsperimeterWeightedWrappingArea;
layer2_summary.sinspostSynapseTouchingArea = sinspostSynapseTouchingArea;
layer2_summary.sinspostSynapseTouchingRatio = sinspostSynapseTouchingRatio;
layer2_summary.sinspreSynapseTouchingArea = sinspreSynapseTouchingArea;
layer2_summary.sinspreSynapseTouchingRatio = sinspreSynapseTouchingRatio;
layer2_summary.singleSynHeadNeckTouchingArea = singleSynHeadNeckTouchingArea;
layer2_summary.singleSynHeadNeckTouchingRatio = singleSynHeadNeckTouchingRatio;
layer2_summary.doubleSynapticCleftSize = doubleSynapticCleftSize;
layer2_summary.dousperimeterRatio = dousperimeterRatio;
layer2_summary.dousperimeterWeightedWrappingArea = dousperimeterWeightedWrappingArea;
layer2_summary.douspostSynapseTouchingArea = douspostSynapseTouchingArea;
layer2_summary.douspostSynapseTouchingRatio = douspostSynapseTouchingRatio;
layer2_summary.douspreSynapseTouchingArea = douspreSynapseTouchingArea;
layer2_summary.douspreSynapseTouchingRatio = douspreSynapseTouchingRatio;
layer2_summary.doubleSynHeadNeckTouchingArea = doubleSynHeadNeckTouchingArea;
layer2_summary.doubleSynHeadNeckTouchingRatio = doubleSynHeadNeckTouchingRatio;
layer2_summary.doubleSynPosition = doubleSynPosition;
layer2_summary.spineCoordinate = spineCoordinate;
layer2_summary.doubleDendriteDensity = doubleDendriteDensity;
layer2_summary.singleDendriteDensity = singleDendriteDensity;








layer4_summary.singleSynHeadVolume = singleSynHeadVolume;
layer4_summary.singleSynMeanHeadRadius = singleSynHeadMeanRadius;
layer4_summary.singleSynNeckLength = singleSynNeckLength;
layer4_summary.singleSynNeckSection = singleSynNeckSection;
layer4_summary.singleSynNeckMeanRadius = singleSynNeckMeanRadius;
layer4_summary.doubleSynHeadVolume = doubleSynHeadVolume;
layer4_summary.doubleSynMeanHeadRadius = doubleSynMeanHeadRadius;
layer4_summary.doubleSynNeckLength = doubleSynNeckLength;
layer4_summary.doubleSynNeckSection = doubleSynNeckSection;
layer4_summary.doubleSynNeckMeanRadius = doubleSynNeckMeanRadius;
layer4_summary.singleSynapticCleftSize = singleSynapticCleftSize;
layer4_summary.sinsperimeterRatio = sinsperimeterRatio;
layer4_summary.sinsperimeterWeightedWrappingArea = sinsperimeterWeightedWrappingArea;
layer4_summary.sinspostSynapseTouchingArea = sinspostSynapseTouchingArea;
layer4_summary.sinspostSynapseTouchingRatio = sinspostSynapseTouchingRatio;
layer4_summary.sinspreSynapseTouchingArea = sinspreSynapseTouchingArea;
layer4_summary.sinspreSynapseTouchingRatio = sinspreSynapseTouchingRatio;
layer4_summary.singleSynHeadNeckTouchingArea = singleSynHeadNeckTouchingArea;
layer4_summary.singleSynHeadNeckTouchingRatio = singleSynHeadNeckTouchingRatio;
layer4_summary.doubleSynapticCleftSize = doubleSynapticCleftSize;
layer4_summary.dousperimeterRatio = dousperimeterRatio;
layer4_summary.dousperimeterWeightedWrappingArea = dousperimeterWeightedWrappingArea;
layer4_summary.douspostSynapseTouchingArea = douspostSynapseTouchingArea;
layer4_summary.douspostSynapseTouchingRatio = douspostSynapseTouchingRatio;
layer4_summary.douspreSynapseTouchingArea = douspreSynapseTouchingArea;
layer4_summary.douspreSynapseTouchingRatio = douspreSynapseTouchingRatio;
layer4_summary.doubleSynHeadNeckTouchingArea = doubleSynHeadNeckTouchingArea;
layer4_summary.doubleSynHeadNeckTouchingRatio = doubleSynHeadNeckTouchingRatio;
layer4_summary.doubleSynPosition = doubleSynPosition;
layer4_summary.spineCoordinate = spineCoordinate;
layer4_summary.doubleDendriteDensity = doubleDendriteDensity;
layer4_summary.singleDendriteDensity = singleDendriteDensity;









layer5_summary.singleSynHeadVolume = singleSynHeadVolume;
layer5_summary.singleSynMeanHeadRadius = singleSynHeadMeanRadius;
layer5_summary.singleSynNeckLength = singleSynNeckLength;
layer5_summary.singleSynNeckSection = singleSynNeckSection;
layer5_summary.singleSynNeckMeanRadius = singleSynNeckMeanRadius;
layer5_summary.doubleSynHeadVolume = doubleSynHeadVolume;
layer5_summary.doubleSynMeanHeadRadius = doubleSynMeanHeadRadius;
layer5_summary.doubleSynNeckLength = doubleSynNeckLength;
layer5_summary.doubleSynNeckSection = doubleSynNeckSection;
layer5_summary.doubleSynNeckMeanRadius = doubleSynNeckMeanRadius;
layer5_summary.singleSynapticCleftSize = singleSynapticCleftSize;
layer5_summary.sinsperimeterRatio = sinsperimeterRatio;
layer5_summary.sinsperimeterWeightedWrappingArea = sinsperimeterWeightedWrappingArea;
layer5_summary.sinspostSynapseTouchingArea = sinspostSynapseTouchingArea;
layer5_summary.sinspostSynapseTouchingRatio = sinspostSynapseTouchingRatio;
layer5_summary.sinspreSynapseTouchingArea = sinspreSynapseTouchingArea;
layer5_summary.sinspreSynapseTouchingRatio = sinspreSynapseTouchingRatio;
layer5_summary.singleSynHeadNeckTouchingArea = singleSynHeadNeckTouchingArea;
layer5_summary.singleSynHeadNeckTouchingRatio = singleSynHeadNeckTouchingRatio;
layer5_summary.doubleSynapticCleftSize = doubleSynapticCleftSize;
layer5_summary.dousperimeterRatio = dousperimeterRatio;
layer5_summary.dousperimeterWeightedWrappingArea = dousperimeterWeightedWrappingArea;
layer5_summary.douspostSynapseTouchingArea = douspostSynapseTouchingArea;
layer5_summary.douspostSynapseTouchingRatio = douspostSynapseTouchingRatio;
layer5_summary.douspreSynapseTouchingArea = douspreSynapseTouchingArea;
layer5_summary.douspreSynapseTouchingRatio = douspreSynapseTouchingRatio;
layer5_summary.doubleSynHeadNeckTouchingArea = doubleSynHeadNeckTouchingArea;
layer5_summary.doubleSynHeadNeckTouchingRatio = doubleSynHeadNeckTouchingRatio;
layer5_summary.doubleSynPosition = doubleSynPosition;
layer5_summary.spineCoordinate = spineCoordinate;

load('layer2_summary.mat')
load('layer4_summary.mat')


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
% doubleSynPosition4 = layer4_summary.doubleSynPosition;
% spineCoordinate4  =layer4_summary.spineCoordinate;
% singleDendriteDensity4 = layer4_summary.singleDendriteDensity;
% doubleDendriteDensity4 = layer4_summary.doubleDendriteDensity;



%% density 
unique_density2 = unique(singleDendriteDensity2);
ratio_density2 = zeros(length(unique_density2),1);

for i = 1:length(unique_density2)
    ratio_density2(i) = sum(doubleDendriteDensity2 == unique_density2(i))/sum(singleDendriteDensity2 == unique_density2(i));
end
figure;plot(unique_density2, ratio_density2)
unique_density4 = unique(singleDendriteDensity4);
ratio_density4 = zeros(length(unique_density4),1);

for i = 1:length(unique_density4)
    ratio_density4(i) = sum(doubleDendriteDensity4 == unique_density4(i))/sum(singleDendriteDensity4 == unique_density4(i));
end
 figure;plot(unique_density4, ratio_density4)

unique_density5 = unique(singleDendriteDensity5);
ratio_density5 = zeros(length(unique_density5),1);

for i = 1:length(unique_density5)
    ratio_density5(i) = sum(doubleDendriteDensity5 == unique_density5(i))/sum(singleDendriteDensity5 == unique_density5(i));
end

figure;plot(unique_density5, ratio_density5)


figure; scatter(singleDendriteDensity2, singleSynHeadVolume2); hold on; scatter(singleDendriteDensity4, singleSynHeadVolume4);
hold on; scatter(singleDendriteDensity5, singleSynHeadVolume5);
%% head volume 
id2 = find(singleSynHeadVolume2~=0);
id4 = find(singleSynHeadVolume4 ~= 0);
id5 = find(singleSynHeadVolume5 ~= 0);
x = [singleSynHeadVolume2(id2);singleSynHeadVolume4(id4);singleSynHeadVolume5(id5)];
x = x.*16.*16.*40/10^9;
y = [zeros(length(singleSynHeadVolume2(id2)),1);ones(length(singleSynHeadVolume4(id4)),1);ones(length(singleSynHeadVolume5(id5)),1).*2];

figure; violinplot(x,y,'ShowMean',true, 'ShowMedian', true); title(['volume of dendrite spine head'])

set(gca,'xtick',[1,2,3],'xticklabel',{'L2/3', 'L4', 'L5'})
ylabel('volume(\mu m^3)')

% figure; histogram(singleSynHeadVolume2(id2), 'Normalization','probability');
% hold on; histogram(singleSynHeadVolume4(id4), 'Normalization','probability');
% hold on; histogram(singleSynHeadVolume5(id5), 'Normalization','probability');

% id2 = randperm(length(singleSynHeadVolume2), 3000);
% id4 = randperm(length(singleSynHeadVolume4), 3000);
% id5 = randperm(length(singleSynHeadVolume5), 3000);
[h,p1,ci,stats] = ttest2(singleSynHeadVolume2(id2),singleSynHeadVolume4(id4));
[h,p2,ci,stats] = ttest2(singleSynHeadVolume4(id4),singleSynHeadVolume5(id5));
[h,p3,ci,stats] = ttest2(singleSynHeadVolume2(id2),singleSynHeadVolume5(id5));
sigstar({[1,2],[2,3], [1,3]},[p1,p2,p3])
x = [doubleSynHeadVolume2(1:2:end);doubleSynHeadVolume4(1:2:end);doubleSynHeadVolume5(1:2:end)];
y = [zeros(length(doubleSynHeadVolume2(1:2:end)),1);ones(length(doubleSynHeadVolume4(1:2:end)),1);ones(length(doubleSynHeadVolume5(1:2:end)),1).*2];

figure; violinplot(x,y); title(['doubleHeadVolume'])
set(gca,'xtick',[1,2,3],'xticklabel',{'L2', 'L4', 'L5'})



x = [singleSynHeadVolume2(:);doubleSynHeadVolume2(1:2:end)];
y = [zeros(length(singleSynHeadVolume2(:)),1);ones(length(doubleSynHeadVolume2(1:2:end)),1)];
[h,p,ci,stats] = ttest2(singleSynHeadVolume2(:),doubleSynHeadVolume2(1:2:end));

figure; violinplot(x,y); title(['headVolume layer2',' p value:', num2str(p)])
set(gca,'xtick',[1,2],'xticklabel',{'single', 'double'})


x = [singleSynHeadVolume4(:);doubleSynHeadVolume4(1:2:end)];
y = [zeros(length(singleSynHeadVolume4(:)),1);ones(length(doubleSynHeadVolume4(1:2:end)),1)];
[h,p,ci,stats] = ttest2(singleSynHeadVolume4(:),doubleSynHeadVolume4(1:2:end));

figure; violinplot(x,y); title(['headVolume layer4',' p value:', num2str(p)])
set(gca,'xtick',[1,2],'xticklabel',{'single', 'double'})

x = [singleSynHeadVolume5(:);doubleSynHeadVolume5(1:2:end)];
y = [zeros(length(singleSynHeadVolume5(:)),1);ones(length(doubleSynHeadVolume5(1:2:end)),1)];
[h,p,ci,stats] = ttest2(singleSynHeadVolume5(:),doubleSynHeadVolume5(1:2:end));

figure; violinplot(x,y); title(['headVolume layer5',' p value:', num2str(p)])
set(gca,'xtick',[1,2],'xticklabel',{'single', 'double'})

%% binary yes/ no contact ratio
pre_ratio = [sum(sinspreSynapseTouchingRatio2 == 0)/ length(sinspreSynapseTouchingRatio2), sum(sinspreSynapseTouchingRatio4 == 0)/ length(sinspreSynapseTouchingRatio4), ...
    sum(sinspreSynapseTouchingRatio5 == 0)/ length(sinspreSynapseTouchingRatio5)];

peri_ratio = [sum(sinsperimeterRatio2 == 0)/ length(sinsperimeterRatio2), sum(sinsperimeterRatio4 == 0)/ length(sinsperimeterRatio4), ...
    sum(sinsperimeterRatio5 == 0)/ length(sinsperimeterRatio5)];

post_ratio = [sum(sinspostSynapseTouchingRatio2 == 0)/ length(sinspostSynapseTouchingRatio2), sum(sinspostSynapseTouchingRatio4 == 0)/ length(sinspostSynapseTouchingRatio4), ...
    sum(sinspostSynapseTouchingRatio5 == 0)/ length(sinspostSynapseTouchingRatio5)];

figure; plot(pre_ratio); hold on; plot(peri_ratio); hold on; plot(post_ratio);




%% pre synapse touching ratio
sinspreSynapseTouchingRatio2(sinspreSynapseTouchingRatio2 == 0) = [];
sinspreSynapseTouchingRatio4(sinspreSynapseTouchingRatio4 == 0) = [];
sinspreSynapseTouchingRatio5(sinspreSynapseTouchingRatio5 == 0) = [];

x = [sinspreSynapseTouchingRatio2(:);sinspreSynapseTouchingRatio4(:);sinspreSynapseTouchingRatio5(:)];
y = [zeros(length(sinspreSynapseTouchingRatio2(:)),1);ones(length(sinspreSynapseTouchingRatio4(:)),1);ones(length(sinspreSynapseTouchingRatio5(:)),1).*2];

figure; violinplot(x,y,'ShowMean',true,'ShowMedian',true); title(['Pre synapse contact ratio'])
set(gca,'xtick',[1,2,3],'xticklabel',{'L2/3', 'L4', 'L5'})
% id2 = randperm(length(sinspreSynapseTouchingRatio2), 3000);
% id4 = randperm(length(sinspreSynapseTouchingRatio4), 3000);
% id5 = randperm(length(sinspreSynapseTouchingRatio5), 3000);
[h,p1,ci,stats] = ttest2(sinspreSynapseTouchingRatio2(:),sinspreSynapseTouchingRatio4(:));
[h,p2,ci,stats] = ttest2(sinspreSynapseTouchingRatio4(:),sinspreSynapseTouchingRatio5(:));
[h,p3,ci,stats] = ttest2(sinspreSynapseTouchingRatio2(:),sinspreSynapseTouchingRatio5(:));
sigstar({[1,2],[2,3], [1,3]},[p1,p2,p3])
% figure; violinplot(x,y); title(['Pre synapse contact ratio'])
% set(gca,'xtick',[1,2,3],'xticklabel',{'L2/3', 'L4', 'L5'})



douspreSynapseTouchingRatio2(douspreSynapseTouchingRatio2 == 0) = [];
douspreSynapseTouchingRatio4(douspreSynapseTouchingRatio4 == 0) = [];
douspreSynapseTouchingRatio5(douspreSynapseTouchingRatio5 == 0) = [];

x = [douspreSynapseTouchingRatio2(:);douspreSynapseTouchingRatio4(:);douspreSynapseTouchingRatio5(:)];
y = [zeros(length(douspreSynapseTouchingRatio2(:)),1);ones(length(douspreSynapseTouchingRatio4(:)),1);ones(length(douspreSynapseTouchingRatio5(:)),1).*2];

figure; violinplot(x,y); title(['douspreSynapseTouchingRatio'])
set(gca,'xtick',[1,2,3],'xticklabel',{'L2', 'L4', 'L5'})


x = [douspreSynapseTouchingRatio2(:);douspreSynapseTouchingRatio4(:)];
y = [zeros(length(douspreSynapseTouchingRatio2(:)),1);ones(length(douspreSynapseTouchingRatio4(:)),1)];
[h,p,ci,stats] = ttest2(douspreSynapseTouchingRatio2(:),douspreSynapseTouchingRatio4(:));

figure; violinplot(x,y); title(['douspreSynapseTouchingRatio',' p value:', num2str(p)])
set(gca,'xtick',[1,2],'xticklabel',{'L2', 'L4'})



x = [douspreSynapseTouchingRatio4(:);douspreSynapseTouchingRatio5(:)];
y = [zeros(length(douspreSynapseTouchingRatio4(:)),1);ones(length(douspreSynapseTouchingRatio5(:)),1)];
[h,p,ci,stats] = ttest2(douspreSynapseTouchingRatio4(:),douspreSynapseTouchingRatio5(:));

figure; violinplot(x,y); title(['douspreSynapseTouchingRatio',' p value:', num2str(p)])
set(gca,'xtick',[1,2],'xticklabel',{'L4', 'L5'})

x = [sinspreSynapseTouchingRatio2(:);sinspreSynapseTouchingRatio4(:)];
y = [zeros(length(sinspreSynapseTouchingRatio2(:)),1);ones(length(sinspreSynapseTouchingRatio4(:)),1)];
[h,p,ci,stats] = ttest2(sinspreSynapseTouchingRatio2(:),sinspreSynapseTouchingRatio4(:));

figure; violinplot(x,y); title(['sinspreSynapseTouchingRatio',' p value:', num2str(p)])
set(gca,'xtick',[1,2],'xticklabel',{'L2', 'L4'})


x = [sinspreSynapseTouchingRatio4(:);sinspreSynapseTouchingRatio5(:)];
y = [zeros(length(sinspreSynapseTouchingRatio4(:)),1);ones(length(sinspreSynapseTouchingRatio5(:)),1)];
[h,p,ci,stats] = ttest2(sinspreSynapseTouchingRatio4(:),sinspreSynapseTouchingRatio5(:));

figure; violinplot(x,y); title(['sinspreSynapseTouchingRatio',' p value:', num2str(p)])
set(gca,'xtick',[1,2],'xticklabel',{'L4', 'L5'})


x = [sinspreSynapseTouchingRatio2(:);douspreSynapseTouchingRatio2(:)];
y = [zeros(length(sinspreSynapseTouchingRatio2(:)),1);ones(length(douspreSynapseTouchingRatio2(:)),1)];
[h,p,ci,stats] = ttest2(sinspreSynapseTouchingRatio2(:),douspreSynapseTouchingRatio2(:));

figure; violinplot(x,y); title(['layer2 preSynapseTouchingRatio',' p value:', num2str(p)])
set(gca,'xtick',[1,2],'xticklabel',{'single', 'double'})

x = [sinspreSynapseTouchingRatio4(:);douspreSynapseTouchingRatio4(:)];
y = [zeros(length(sinspreSynapseTouchingRatio4(:)),1);ones(length(douspreSynapseTouchingRatio4(:)),1)];
[h,p,ci,stats] = ttest2(sinspreSynapseTouchingRatio4(:),douspreSynapseTouchingRatio4(:));

figure; violinplot(x,y); title(['layer4 preSynapseTouchingRatio',' p value:', num2str(p)])
set(gca,'xtick',[1,2],'xticklabel',{'single', 'double'})


x = [sinspreSynapseTouchingRatio5(:);douspreSynapseTouchingRatio5(:)];
y = [zeros(length(sinspreSynapseTouchingRatio5(:)),1);ones(length(douspreSynapseTouchingRatio5(:)),1)];
[h,p,ci,stats] = ttest2(sinspreSynapseTouchingRatio5(:),douspreSynapseTouchingRatio5(:));

figure; violinplot(x,y); title(['layer5 preSynapseTouchingRatio',' p value:', num2str(p)])
set(gca,'xtick',[1,2],'xticklabel',{'single', 'double'})
%% perimeter touching ratio
sinsperimeterRatio2(sinsperimeterRatio2 == 0) = [];
sinsperimeterRatio4(sinsperimeterRatio4 == 0) = [];
sinsperimeterRatio5(sinsperimeterRatio5 == 0) = [];
x = [sinsperimeterRatio2(:);sinsperimeterRatio4(:);sinsperimeterRatio5(:)];
y = [zeros(length(sinsperimeterRatio2(:)),1);ones(length(sinsperimeterRatio4(:)),1);ones(length(sinsperimeterRatio5(:)),1).*2];
figure; violinplot(x,y,'ShowMean',true,'ShowMedian',true); title(['perimeter ratio'])
set(gca,'xtick',[1,2,3],'xticklabel',{'L2/3', 'L4', 'L5'})
% id2 = randperm(length(sinsperimeterRatio2), 3000);
% id4 = randperm(length(sinsperimeterRatio4), 3000);
% id5 = randperm(length(sinsperimeterRatio5), 3000);
[h,p1,ci,stats] = ttest2(sinsperimeterRatio2,sinsperimeterRatio4);
[h,p2,ci,stats] = ttest2(sinsperimeterRatio4,sinsperimeterRatio5);
[h,p3,ci,stats] = ttest2(sinsperimeterRatio2,sinsperimeterRatio5);
sigstar({[1,2],[2,3], [1,3]},[p1,p2,p3])




figure; violinplot(x,y,'ShowMean',true,'ShowMedian',false); title(['dendrite spine head radius'])
set(gca,'xtick',[1,2,3],'xticklabel',{'L2/3', 'L4', 'L5'}); ylabel('radius(nm)')
% figure; violinplot(x,y,'ShowMean',true,'ShowMedian',false); title(['dendrite spine neck section'])
% set(gca,'xtick',[1,2,3],'xticklabel',{'L2/3', 'L4', 'L5'})
id2 = randperm(length(singleSynMeanHeadRadius2), 3000);
id4 = randperm(length(singleSynMeanHeadRadius4), 3000);
id5 = randperm(length(singleSynMeanHeadRadius5), 3000);
[h,p1,ci,stats] = ttest2(singleSynMeanHeadRadius2(id2),singleSynMeanHeadRadius4(id4));
[h,p2,ci,stats] = ttest2(singleSynMeanHeadRadius4(id4),singleSynMeanHeadRadius5(id5));
[h,p3,ci,stats] = ttest2(singleSynMeanHeadRadius2(id2),singleSynMeanHeadRadius5(id5));
sigstar({[1,2],[2,3], [1,3]},[p1,p2,p3])








x = [sinsperimeterRatio2(:);sinsperimeterRatio4(:);];
y = [zeros(length(sinsperimeterRatio2(:)),1);ones(length(sinsperimeterRatio4(:)),1)];
[h,p,ci,stats] = ttest2(sinsperimeterRatio2(:),sinsperimeterRatio4(:));

figure; violinplot(x,y); title(['sinsperimeterRatio',' p value:', num2str(p)])
set(gca,'xtick',[1,2],'xticklabel',{'L2', 'L4'})


x = [sinsperimeterRatio4(:);sinsperimeterRatio5(:);];
y = [zeros(length(sinsperimeterRatio4(:)),1);ones(length(sinsperimeterRatio5(:)),1)];
[h,p,ci,stats] = ttest2(sinsperimeterRatio4(:),sinsperimeterRatio5(:));

figure; violinplot(x,y); title(['sinsperimeterRatio',' p value:', num2str(p)])
set(gca,'xtick',[1,2],'xticklabel',{'L4', 'L5'})




dousperimeterRatio2(dousperimeterRatio2(:) == 0) = [];
dousperimeterRatio4(dousperimeterRatio4(:) == 0) = [];
dousperimeterRatio5(dousperimeterRatio5(:) == 0) = [];
x = [dousperimeterRatio2(:);dousperimeterRatio4(:);dousperimeterRatio5(:)];
y = [zeros(length(dousperimeterRatio2(:)),1);ones(length(dousperimeterRatio4(:)),1);ones(length(dousperimeterRatio5(:)),1).*2];

figure; violinplot(x,y); title(['dousperimeterRatio'])
set(gca,'xtick',[1,2,3],'xticklabel',{'L2', 'L4', 'L5'})


x = [dousperimeterRatio2(:);dousperimeterRatio4(:);];
y = [zeros(length(dousperimeterRatio2(:)),1);ones(length(dousperimeterRatio4(:)),1)];
[h,p,ci,stats] = ttest2(dousperimeterRatio2(:),dousperimeterRatio4(:));

figure; violinplot(x,y); title(['dousperimeterRatio',' p value:', num2str(p)])
set(gca,'xtick',[1,2],'xticklabel',{'L2', 'L4'})


x = [dousperimeterRatio4(:);dousperimeterRatio5(:);];
y = [zeros(length(dousperimeterRatio4(:)),1);ones(length(dousperimeterRatio5(:)),1)];
[h,p,ci,stats] = ttest2(dousperimeterRatio4(:),dousperimeterRatio5(:));

figure; violinplot(x,y); title(['dousperimeterRatio',' p value:', num2str(p)])
set(gca,'xtick',[1,2],'xticklabel',{'L4', 'L5'})

x = [sinsperimeterRatio2(:);dousperimeterRatio2(:);];
y = [zeros(length(sinsperimeterRatio2(:)),1);ones(length(dousperimeterRatio2(:)),1)];
[h,p,ci,stats] = ttest2(sinsperimeterRatio2(:),dousperimeterRatio2(:));

figure; violinplot(x,y); title(['layer2 perimeterRatio',' p value:', num2str(p)])
set(gca,'xtick',[1,2],'xticklabel',{'single', 'double'})


x = [sinsperimeterRatio4(:);dousperimeterRatio4(:);];
y = [zeros(length(sinsperimeterRatio4(:)),1);ones(length(dousperimeterRatio4(:)),1)];
[h,p,ci,stats] = ttest2(sinsperimeterRatio4(:),dousperimeterRatio4(:));

figure; violinplot(x,y); title(['layer4 perimeterRatio',' p value:', num2str(p)])
set(gca,'xtick',[1,2],'xticklabel',{'single', 'double'})


x = [sinsperimeterRatio5(:);dousperimeterRatio5(:);];
y = [zeros(length(sinsperimeterRatio5(:)),1);ones(length(dousperimeterRatio5(:)),1)];
[h,p,ci,stats] = ttest2(sinsperimeterRatio5(:),dousperimeterRatio5(:));

figure; violinplot(x,y); title(['layer5 perimeterRatio',' p value:', num2str(p)])
set(gca,'xtick',[1,2],'xticklabel',{'single', 'double'})

%% Post synapse touching ratio
douspostSynapseTouchingRatio2(douspostSynapseTouchingRatio2(:) == 0) = [];
douspostSynapseTouchingRatio4(douspostSynapseTouchingRatio4(:) == 0) = [];
douspostSynapseTouchingRatio5(douspostSynapseTouchingRatio5(:) == 0) = [];
x = [douspostSynapseTouchingRatio2(:);douspostSynapseTouchingRatio4(:);douspostSynapseTouchingRatio5(:)];
y = [zeros(length(douspostSynapseTouchingRatio2(:)),1);ones(length(douspostSynapseTouchingRatio4(:)),1);ones(length(douspostSynapseTouchingRatio5(:)),1).*2];

figure; violinplot(x,y,'ShowMean',true,'ShowMedian',false); title(['douspostSynapseTouchingRatio'])
set(gca,'xtick',[1,2,3],'xticklabel',{'L2', 'L4', 'L5'})
% id2 = randperm(length(douspostSynapseTouchingRatio2), 3000);
% id4 = randperm(length(douspostSynapseTouchingRatio4), 3000);
% id5 = randperm(length(douspostSynapseTouchingRatio5), 3000);
[h,p1,ci,stats] = ttest2(douspostSynapseTouchingRatio2(id2),douspostSynapseTouchingRatio4(id4));
[h,p2,ci,stats] = ttest2(douspostSynapseTouchingRatio4(id4),douspostSynapseTouchingRatio5(id5));
[h,p3,ci,stats] = ttest2(douspostSynapseTouchingRatio2(id2),douspostSynapseTouchingRatio5(id5));
sigstar({[1,2],[2,3], [1,3]},[p1,p2,p3])

x = [douspostSynapseTouchingRatio2(:);douspostSynapseTouchingRatio4(:);];
y = [zeros(length(douspostSynapseTouchingRatio2(:)),1);ones(length(douspostSynapseTouchingRatio4(:)),1)];
[h,p,ci,stats] = ttest2(douspostSynapseTouchingRatio2(:),douspostSynapseTouchingRatio4(:));

figure; violinplot(x,y); title(['douspostSynapseTouchingRatio',' p value:', num2str(p)])
set(gca,'xtick',[1,2],'xticklabel',{'L2', 'L4'})

x = [douspostSynapseTouchingRatio4(:);douspostSynapseTouchingRatio5(:);];
y = [zeros(length(douspostSynapseTouchingRatio4(:)),1);ones(length(douspostSynapseTouchingRatio5(:)),1)];
[h,p,ci,stats] = ttest2(douspostSynapseTouchingRatio4(:),douspostSynapseTouchingRatio5(:));

figure; violinplot(x,y); title(['douspostSynapseTouchingRatio',' p value:', num2str(p)])
set(gca,'xtick',[1,2],'xticklabel',{'L4', 'L5'})

sinspostSynapseTouchingRatio2(sinspostSynapseTouchingRatio2 == 0) = [];
sinspostSynapseTouchingRatio4(sinspostSynapseTouchingRatio4 == 0)= [];
sinspostSynapseTouchingRatio5(sinspostSynapseTouchingRatio5 == 0) = [];
x = [sinspostSynapseTouchingRatio2(:);sinspostSynapseTouchingRatio4(:);sinspostSynapseTouchingRatio5(:)];
y = [zeros(length(sinspostSynapseTouchingRatio2(:)),1);ones(length(sinspostSynapseTouchingRatio4(:)),1);ones(length(sinspostSynapseTouchingRatio5(:)),1).*2];

figure; violinplot(x,y,'ShowMean',true,'ShowMedian',true); title(['post-synapse contact ratio'])
set(gca,'xtick',[1,2,3],'xticklabel',{'L2/3', 'L4', 'L5'})
% id2 = randperm(length(sinspostSynapseTouchingRatio2), 3000);
% id4 = randperm(length(sinspostSynapseTouchingRatio4), 3000);
% id5 = randperm(length(sinspostSynapseTouchingRatio5), 3000);
[h,p1,ci,stats] = ttest2(sinspostSynapseTouchingRatio2,sinspostSynapseTouchingRatio4);
[h,p2,ci,stats] = ttest2(sinspostSynapseTouchingRatio4,sinspostSynapseTouchingRatio5);
[h,p3,ci,stats] = ttest2(sinspostSynapseTouchingRatio2,sinspostSynapseTouchingRatio5);
sigstar({[1,2],[2,3], [1,3]},[p1,p2,p3])


figure; histogram(sinspostSynapseTouchingRatio2,50, "Normalization","probability"); 
hold on; histogram(sinspostSynapseTouchingRatio4,50, "Normalization","probability"); 
hold on; histogram(sinspostSynapseTouchingRatio5,50, "Normalization","probability"); 


x = [sinspostSynapseTouchingRatio2(:);sinspostSynapseTouchingRatio4(:);];
y = [zeros(length(sinspostSynapseTouchingRatio2(:)),1);ones(length(sinspostSynapseTouchingRatio4(:)),1)];
[h,p,ci,stats] = ttest2(sinspostSynapseTouchingRatio2(:),sinspostSynapseTouchingRatio4(:));

figure; violinplot(x,y); title(['sinspostSynapseTouchingRatio',' p value:', num2str(p)])
set(gca,'xtick',[1,2],'xticklabel',{'L2', 'L4'})

x = [sinspostSynapseTouchingRatio4(:);sinspostSynapseTouchingRatio5(:);];
y = [zeros(length(sinspostSynapseTouchingRatio4(:)),1);ones(length(sinspostSynapseTouchingRatio5(:)),1)];
[h,p,ci,stats] = ttest2(sinspostSynapseTouchingRatio4(:),sinspostSynapseTouchingRatio5(:));

figure; violinplot(x,y); title(['sinspostSynapseTouchingRatio',' p value:', num2str(p)])
set(gca,'xtick',[1,2],'xticklabel',{'L4', 'L5'})

x = [sinspostSynapseTouchingRatio2(:);sinspostSynapseTouchingRatio4(:);sinspostSynapseTouchingRatio5(:)];
y = [zeros(length(sinspostSynapseTouchingRatio2(:)),1);ones(length(sinspostSynapseTouchingRatio4(:)),1);ones(length(sinspostSynapseTouchingRatio5(:)),1).*2];

figure; violinplot(x,y); title(['sinspostSynapseTouchingRatio'])
set(gca,'xtick',[1,2,3],'xticklabel',{'L2', 'L4', 'L5'})

x = [sinspostSynapseTouchingRatio2(:);douspostSynapseTouchingRatio2(:);];
y = [zeros(length(sinspostSynapseTouchingRatio2(:)),1);ones(length(douspostSynapseTouchingRatio2(:)),1)];
[h,p,ci,stats] = ttest2(sinspostSynapseTouchingRatio2(:),douspostSynapseTouchingRatio2(:));

figure; violinplot(x,y); title(['L2 postSynapseTouchingRatio',' p value:', num2str(p)])

set(gca,'xtick',[1,2],'xticklabel',{'single', 'double'})

x = [sinspostSynapseTouchingRatio5(:);douspostSynapseTouchingRatio5(:);];
y = [zeros(length(sinspostSynapseTouchingRatio5(:)),1);ones(length(douspostSynapseTouchingRatio5(:)),1)];
[h,p,ci,stats] = ttest2(sinspostSynapseTouchingRatio5(:),douspostSynapseTouchingRatio5(:));

figure; violinplot(x,y); title(['L5 postSynapseTouchingRatio',' p value:', num2str(p)])
set(gca,'xtick',[1,2],'xticklabel',{'single', 'double'})


x = [sinspostSynapseTouchingRatio4(:);douspostSynapseTouchingRatio4(:);];
y = [zeros(length(sinspostSynapseTouchingRatio4(:)),1);ones(length(douspostSynapseTouchingRatio4(:)),1)];
[h,p,ci,stats] = ttest2(sinspostSynapseTouchingRatio4(:),douspostSynapseTouchingRatio4(:));

figure; violinplot(x,y); title(['L4 postSynapseTouchingRatio',' p value:', num2str(p)])
set(gca,'xtick',[1,2],'xticklabel',{'single', 'double'})
% plot the bar plot of the contact, conact in head, contact in neck
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
legend({'L2/3', 'L4', 'L5'}, 'Location', 'best');  % Legend

groupNames = {'Contact', 'contact at head', 'contact at neck'}; % Custom names for each group of bars
xticks(1:length(groupNames)); % Set the x-axis ticks to match the number of groups
xticklabels(groupNames);
%% what is the distribution of the contact ratio between head and neck for
%% the dendrite spines which are in contact both in head and neck
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

figure; violinplot(x,y,'ShowMean',true,'ShowMedian',true); title(['contact ratio of head and neck'])
set(gca,'xtick',[1,2,3],'xticklabel',{'L2/3', 'L4', 'L5'})
% id2 = randperm(length(sinspostSynapseTouchingRatio2), 3000);
% id4 = randperm(length(sinspostSynapseTouchingRatio4), 3000);
% id5 = randperm(length(sinspostSynapseTouchingRatio5), 3000);
[h,p1,ci,stats] = ttest2(ratio2,ratio4);
[h,p2,ci,stats] = ttest2(ratio4,ratio5);
[h,p3,ci,stats] = ttest2(ratio2,ratio5);
sigstar({[1,2],[2,3], [1,3]},[p1,p2,p3])





figure; histogram(log2(ratio2), 50,'Normalization','probability');
hold on; histogram(log2(ratio4), 50,'Normalization','probability');
hold on; histogram(log2(ratio5), 50,'Normalization','probability');
xlabel('log-ratio'); ylabel('frequency')
legend('L2/3', 'L4', 'L5')
title('Ratio of the contact area between head and neck')

%% check if there were any 


%% weighted wrapping area
sinsperimeterWeightedWrappingArea2(sinsperimeterWeightedWrappingArea2 == 0) = [];
sinsperimeterWeightedWrappingArea4(sinsperimeterWeightedWrappingArea4 == 0) = [];
sinsperimeterWeightedWrappingArea5(sinsperimeterWeightedWrappingArea5 == 0) = [];
x = [sinsperimeterWeightedWrappingArea2(:);sinsperimeterWeightedWrappingArea4(:);sinsperimeterWeightedWrappingArea5(:)];
y = [zeros(length(sinsperimeterWeightedWrappingArea2(:)),1);ones(length(sinsperimeterWeightedWrappingArea4(:)),1);ones(length(sinsperimeterWeightedWrappingArea5(:)),1).*2];

figure; violinplot(x,y,'ShowMean',true,'ShowMedian',false); title(['sinsWeightedWrappingArea'])
set(gca,'xtick',[1,2,3],'xticklabel',{'L2/3', 'L4', 'L5'})
figure; histogram(sinsperimeterWeightedWrappingArea2, 100,'Normalization','probability');
hold on; histogram(sinsperimeterWeightedWrappingArea5, 100,'Normalization','probability');

x = [sinsperimeterWeightedWrappingArea4(:);sinsperimeterWeightedWrappingArea5(:);];
y = [zeros(length(sinsperimeterWeightedWrappingArea4(:)),1);ones(length(sinsperimeterWeightedWrappingArea5(:)),1)];
[h,p,ci,stats] = ttest2(sinsperimeterWeightedWrappingArea4(:),sinsperimeterWeightedWrappingArea5(:));

figure; violinplot(x,y); title(['WeightedWrappingArea',' p value:', num2str(p)])
set(gca,'xtick',[1,2],'xticklabel',{'L4', 'L5'})


dousperimeterWeightedWrappingArea2(dousperimeterWeightedWrappingArea2(:) == 0) = [];
dousperimeterWeightedWrappingArea4(dousperimeterWeightedWrappingArea4(:) == 0) = [];
dousperimeterWeightedWrappingArea5(dousperimeterWeightedWrappingArea5(:) == 0) = [];
x = [dousperimeterWeightedWrappingArea2(:);dousperimeterWeightedWrappingArea4(:);dousperimeterWeightedWrappingArea5(:)];
y = [zeros(length(dousperimeterWeightedWrappingArea2(:)),1);ones(length(dousperimeterWeightedWrappingArea4(:)),1);ones(length(dousperimeterWeightedWrappingArea5(:)),1).*2];

figure; violinplot(x,y); title(['dousWeightedWrappingArea'])
set(gca,'xtick',[1,2,3],'xticklabel',{'L2', 'L4', 'L5'})


%% cleft size
x = [singleSynapticCleftSize2(id2);singleSynapticCleftSize4(id4);singleSynapticCleftSize5(id5)];
y = [zeros(length(singleSynapticCleftSize2(id2)),1);ones(length(singleSynapticCleftSize4(id4)),1);ones(length(singleSynapticCleftSize5(id5)),1).*2];

figure; violinplot(x,y,'ShowMean',true,'ShowMedian',false); title(['PSD size'])
set(gca,'xtick',[1,2,3],'xticklabel',{'L2/3', 'L4', 'L5'})

% figure; violinplot(x,y,'ShowMean',true,'ShowMedian',false); title(['dendrite spine head radius'])
% set(gca,'xtick',[1,2,3],'xticklabel',{'L2/3', 'L4', 'L5'}); ylabel('radius(nm)')
% figure; violinplot(x,y,'ShowMean',true,'ShowMedian',false); title(['dendrite spine neck section'])
% set(gca,'xtick',[1,2,3],'xticklabel',{'L2/3', 'L4', 'L5'})
id2 = randperm(length(singleSynapticCleftSize2), 3000);
id4 = randperm(length(singleSynapticCleftSize4), 3000);
id5 = randperm(length(singleSynapticCleftSize5), 3000);
[h,p1,ci,stats] = ttest2(singleSynapticCleftSize2(id2),singleSynapticCleftSize4(id4));
[h,p2,ci,stats] = ttest2(singleSynapticCleftSize4(id4),singleSynapticCleftSize5(id5));
[h,p3,ci,stats] = ttest2(singleSynapticCleftSize2(id2),singleSynapticCleftSize5(id5));
sigstar({[1,2],[2,3], [1,3]},[p1,p2,p3])



x = [singleSynapticCleftSize2(:);singleSynapticCleftSize4(:)];
y = [zeros(length(singleSynapticCleftSize2(:)),1);ones(length(singleSynapticCleftSize4(:)),1)];
[h,p,ci,stats] = ttest2(singleSynapticCleftSize2(:),singleSynapticCleftSize4(:));
figure; violinplot(x,y); title(['singleSynapticCleftSize',' p value:', num2str(p)])
set(gca,'xtick',[1,2],'xticklabel',{'L2', 'L4'})

x = [singleSynapticCleftSize4(:);singleSynapticCleftSize5(:)];
y = [zeros(length(singleSynapticCleftSize4(:)),1);ones(length(singleSynapticCleftSize5(:)),1)];
[h,p,ci,stats] = ttest2(singleSynapticCleftSize4(:),singleSynapticCleftSize5(:));
figure; violinplot(x,y); title(['singleSynapticCleftSize',' p value:', num2str(p)])
set(gca,'xtick',[1,2],'xticklabel',{'L4', 'L5'})

x = [doubleSynapticCleftSize2(:);doubleSynapticCleftSize4(:)];
y = [zeros(length(doubleSynapticCleftSize2(:)),1);ones(length(doubleSynapticCleftSize4(:)),1)];
[h,p,ci,stats] = ttest2(doubleSynapticCleftSize2(:),doubleSynapticCleftSize4(:));
figure; violinplot(x,y); title(['doubleSynapticCleftSize',' p value:', num2str(p)])
set(gca,'xtick',[1,2],'xticklabel',{'L2', 'L4'})


x = [doubleSynapticCleftSize4(:);doubleSynapticCleftSize5(:)];
y = [zeros(length(doubleSynapticCleftSize4(:)),1);ones(length(doubleSynapticCleftSize5(:)),1)];
[h,p,ci,stats] = ttest2(doubleSynapticCleftSize4(:),doubleSynapticCleftSize5(:));
figure; violinplot(x,y); title(['doubleSynapticCleftSize',' p value:', num2str(p)])
set(gca,'xtick',[1,2],'xticklabel',{'L4', 'L5'})

x = [doubleSynapticCleftSize2(:);doubleSynapticCleftSize4(:);doubleSynapticCleftSize5(:)];
y = [zeros(length(doubleSynapticCleftSize2(:)),1);ones(length(doubleSynapticCleftSize4(:)),1);ones(length(doubleSynapticCleftSize5(:)),1).*2];

figure; violinplot(x,y); title(['doubleSynapticCleftSize'])
set(gca,'xtick',[1,2,3],'xticklabel',{'L2', 'L4', 'L5'})

%% head volume
singleSynHeadVolume2(singleSynHeadVolume2 == 0) = [];
singleSynHeadVolume4(singleSynHeadVolume4 == 0) = [];
singleSynHeadVolume5(singleSynHeadVolume5 == 0) = [];
x = [singleSynHeadVolume2(:);singleSynHeadVolume4(:);singleSynHeadVolume5(:)];
y = [zeros(length(singleSynHeadVolume2(:)),1);ones(length(singleSynHeadVolume4(:)),1);ones(length(singleSynHeadVolume5(:)),1).*2];

figure; violinplot(x,y); title(['singleSynHeadVolume'])
set(gca,'xtick',[1,2,3],'xticklabel',{'L2', 'L4', 'L5'})



x = [doubleSynHeadVolume2(:);doubleSynHeadVolume4(:);doubleSynHeadVolume5(:)];
y = [zeros(length(doubleSynHeadVolume2(:)),1);ones(length(doubleSynHeadVolume4(:)),1);ones(length(doubleSynHeadVolume5(:)),1).*2];

figure; violinplot(x,y); title(['doubleSynHeadVolume'])
set(gca,'xtick',[1,2,3],'xticklabel',{'L2', 'L4', 'L5'})




x = [doubleSynHeadVolume2(:);doubleSynHeadVolume4(:)];
y = [zeros(length(doubleSynHeadVolume2(:)),1);ones(length(doubleSynHeadVolume4(:)),1)];
[h,p,ci,stats] = ttest2(doubleSynHeadVolume2(:),doubleSynHeadVolume4(:));

figure; violinplot(x,y); title(['doubleSynHeadVolume', 'p value:', num2str(p)])
set(gca,'xtick',[1,2],'xticklabel',{'L2', 'L4'})

x = [doubleSynHeadVolume4(:);doubleSynHeadVolume5(:)];
y = [zeros(length(doubleSynHeadVolume4(:)),1);ones(length(doubleSynHeadVolume5(:)),1)];
[h,p,ci,stats] = ttest2(doubleSynHeadVolume4(:),doubleSynHeadVolume5(:));

figure; violinplot(x,y); title(['doubleSynHeadVolume', 'p value:', num2str(p)])
set(gca,'xtick',[1,2],'xticklabel',{'L4', 'L5'})


%% neck length

x = [singleSynNeckLength2(:);singleSynNeckLength4(:);];
y = [zeros(length(singleSynNeckLength2(:)),1);ones(length(singleSynNeckLength4(:)),1)];
[h,p,ci,stats] = ttest2(singleSynNeckLength2(:),singleSynNeckLength4(:));

figure; violinplot(x,y); title(['singleSynNeckLength',' p value:', num2str(p)])
set(gca,'xtick',[1,2],'xticklabel',{'L2', 'L4'})


x = [singleSynNeckLength4(:);singleSynNeckLength5(:);];
y = [zeros(length(singleSynNeckLength4(:)),1);ones(length(singleSynNeckLength5(:)),1)];
[h,p,ci,stats] = ttest2(singleSynNeckLength4(:),singleSynNeckLength5(:));

figure; violinplot(x,y); title(['singleSynNeckLength','p value:', num2str(p)])
set(gca,'xtick',[1,2],'xticklabel',{'L4', 'L5'})


singleSynNeckLength2(singleSynNeckLength2(:) == 0) = [];
singleSynNeckLength4(singleSynNeckLength4(:) == 0) = [];
singleSynNeckLength5(singleSynNeckLength5(:) == 0) = [];
% singleSynNeckLength2(singleSynNeckLength2 > 4000)  = [];
% singleSynNeckLength4(singleSynNeckLength4 > 4000)  = [];
% singleSynNeckLength5(singleSynNeckLength5 > 4000)  = [];


x = [singleSynNeckLength2(:);singleSynNeckLength4(:);singleSynNeckLength5(:)];
% x = x*sqrt(16^2 + 40^2)/1000;
y = [zeros(length(singleSynNeckLength2(:)),1);ones(length(singleSynNeckLength4(:)),1);ones(length(singleSynNeckLength5(:)),1).*2];

figure; violinplot(x,y,'ShowMean',true,'ShowMedian',true); title(['dendrite spine neck length'])
set(gca,'xtick',[1,2,3],'xticklabel',{'L2/3', 'L4', 'L5'}); ylabel('\mu m')
% id2 = randperm(length(singleSynNeckLength2), 3000);
% id4 = randperm(length(singleSynNeckLength4), 3000);
% id5 = randperm(length(singleSynNeckLength5), 3000);
[h,p1,ci,stats] = ttest2(singleSynNeckLength2(:),singleSynNeckLength4(:));
[h,p2,ci,stats] = ttest2(singleSynNeckLength4(:),singleSynNeckLength5(:));
[h,p3,ci,stats] = ttest2(singleSynNeckLength2(:),singleSynNeckLength5(:));
sigstar({[1,2],[2,3], [1,3]},[p1,p2,p3])

doubleSynNeckLength2(doubleSynNeckLength2(:) == 0) = [];
doubleSynNeckLength4(doubleSynNeckLength4(:) == 0) = [];
doubleSynNeckLength5(doubleSynNeckLength5(:) == 0) = [];
x = [doubleSynNeckLength2(:);doubleSynNeckLength4(:);doubleSynNeckLength5(:)];
y = [zeros(length(doubleSynNeckLength2(:)),1);ones(length(doubleSynNeckLength4(:)),1);ones(length(doubleSynNeckLength5(:)),1).*2];

figure; violinplot(x,y); title(['doubleSynNeckLength'])
set(gca,'xtick',[1,2,3],'xticklabel',{'L2', 'L4', 'L5'})

x = [doubleSynNeckLength4(:);doubleSynNeckLength5(:);];
y = [zeros(length(doubleSynNeckLength4(:)),1);ones(length(doubleSynNeckLength5(:)),1)];
[h,p,ci,stats] = ttest2(doubleSynNeckLength4(:),doubleSynNeckLength5(:));

figure; violinplot(x,y); title(['doubleSynNeckLength','p value:', num2str(p)])
set(gca,'xtick',[1,2],'xticklabel',{'L4', 'L5'})
%% neck section area

singleSynNeckSection2(isinf(singleSynNeckSection2)) = [];
singleSynNeckSection2(singleSynNeckSection2 > 1e5) = [];
singleSynNeckSection4(isinf(singleSynNeckSection4)) = [];
singleSynNeckSection4(singleSynNeckSection4 > 1e5) = [];
singleSynNeckSection5(isinf(singleSynNeckSection5)) = [];
singleSynNeckSection5(singleSynNeckSection5 > 1e5) = [];
x = [singleSynNeckSection2(:);singleSynNeckSection4(:);singleSynNeckSection5(:)];
y = [zeros(length(singleSynNeckSection2(:)),1);ones(length(singleSynNeckSection4(:)),1);ones(length(singleSynNeckSection5(:)),1).*2];
x = x.*16.*16/1000/1000;
figure; violinplot(x,y,'ShowMean',true,'ShowMedian',true); title(['dendrite spine neck section'])
set(gca,'xtick',[1,2,3],'xticklabel',{'L2/3', 'L4', 'L5'}); ylabel('\mu m ^{2}')
% id2 = randperm(length(singleSynNeckSection2), 3000);
% id4 = randperm(length(singleSynNeckSection4), 3000);
% id5 = randperm(length(singleSynNeckSection5), 3000);
[h,p1,ci,stats] = ttest2(singleSynNeckSection2(:),singleSynNeckSection4(:));
[h,p2,ci,stats] = ttest2(singleSynNeckSection4(:),singleSynNeckSection5(:));
[h,p3,ci,stats] = ttest2(singleSynNeckSection2(:),singleSynNeckSection5(:));
sigstar({[1,2],[2,3], [1,3]},[p1,p2,p3])






x = [singleSynNeckSection4(:);singleSynNeckSection5(:);];
y = [zeros(length(singleSynNeckSection4(:)),1);ones(length(singleSynNeckSection5(:)),1)];
[h,p,ci,stats] = ttest2(singleSynNeckSection4(:),singleSynNeckSection5(:));
figure; violinplot(x,y); title(['singleSynNeckSection', 'p value:', num2str(p)])
set(gca,'xtick',[1,2],'xticklabel',{'L4', 'L5'})



x = [doubleSynNeckSection2(:);doubleSynNeckSection4(:);];
y = [zeros(length(doubleSynNeckSection2(:)),1);ones(length(doubleSynNeckSection4(:)),1)];
[h,p,ci,stats] = ttest2(doubleSynNeckSection2(:),doubleSynNeckSection4(:));
figure; violinplot(x,y); title(['doubleSynNeckSection', 'p value:', num2str(p)])
set(gca,'xtick',[1,2],'xticklabel',{'L2', 'L4'})

doubleSynNeckSection2(isinf(doubleSynNeckSection2)) = [];
doubleSynNeckSection2(doubleSynNeckSection2 > 1e5) = [];
doubleSynNeckSection4(isinf(doubleSynNeckSection4)) = [];
doubleSynNeckSection4(doubleSynNeckSection4 > 1e5) = [];
doubleSynNeckSection5(isinf(doubleSynNeckSection5)) = [];
doubleSynNeckSection5(doubleSynNeckSection5 > 1e5) = [];
x = [doubleSynNeckSection2(:);doubleSynNeckSection4(:);doubleSynNeckSection5(:)];
y = [zeros(length(doubleSynNeckSection2(:)),1);ones(length(doubleSynNeckSection4(:)),1);ones(length(doubleSynNeckSection5(:)),1).*2];

figure; violinplot(x,y); title(['doubleSynNeckSection'])
set(gca,'xtick',[1,2,3],'xticklabel',{'L2', 'L4', 'L5'})





x = [doubleSynNeckLength2(:);doubleSynNeckLength4(:);doubleSynNeckLength5(:)];
y = [zeros(length(doubleSynNeckLength2(:)),1);ones(length(doubleSynNeckLength4(:)),1);ones(length(doubleSynNeckLength5(:)),1).*2];

figure; violinplot(x,y); title(['doubleSynNeckLength'])
set(gca,'xtick',[1,2,3],'xticklabel',{'L2', 'L4', 'L5'})



x = [doubleSynNeckLength4(:);doubleSynNeckLength5(:);];
y = [zeros(length(doubleSynNeckLength4(:)),1);ones(length(doubleSynNeckLength5(:)),1)];
[h,p,ci,stats] = ttest2(doubleSynNeckLength4(:),doubleSynNeckLength5(:));
figure; violinplot(x,y); title(['doubleSynNeckLength',' p value:', num2str(p)])
set(gca,'xtick',[1,2],'xticklabel',{'L4', 'L5'})
%% head radius
singleSynMeanHeadRadius2 = singleSynMeanHeadRadius2(singleSynMeanHeadRadius2~=0);
singleSynMeanHeadRadius4 = singleSynMeanHeadRadius4(singleSynMeanHeadRadius4~=0);
singleSynMeanHeadRadius5 = singleSynMeanHeadRadius5(singleSynMeanHeadRadius5~=0);
singleSynMeanHeadRadius2(singleSynMeanHeadRadius2>600) = [];
singleSynMeanHeadRadius4(singleSynMeanHeadRadius4>600) = [];
singleSynMeanHeadRadius5(singleSynMeanHeadRadius5>600) = [];
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
figure; histogram(singleSynMeanHeadRadius2,100, "Normalization","probability"); 
hold on; histogram(singleSynMeanHeadRadius4,100, "Normalization","probability"); 
hold on; histogram(singleSynMeanHeadRadius5,100, "Normalization","probability"); 

legend('L2/3', 'L4', 'L5')
title('histogram of the radius for the max section in head')
xlabel('nm')

x = [doubleSynMeanHeadRadius2(:);doubleSynMeanHeadRadius4(:);doubleSynMeanHeadRadius5(:)];
y = [zeros(length(doubleSynMeanHeadRadius2(:)),1);ones(length(doubleSynMeanHeadRadius4(:)),1);ones(length(doubleSynMeanHeadRadius5(:)),1).*2];

figure; violinplot(x,y); title(['doubleSynMeanHeadRadius'])
set(gca,'xtick',[1,2,3],'xticklabel',{'L2', 'L4', 'L5'})



x = [doubleSynMeanHeadRadius4(:);doubleSynMeanHeadRadius5(:)];
y = [zeros(length(doubleSynMeanHeadRadius4(:)),1);ones(length(doubleSynMeanHeadRadius5(:)),1)];
[h,p,ci,stats] = ttest2(doubleSynMeanHeadRadius4(:),doubleSynMeanHeadRadius5(:));
figure; violinplot(x,y); title(['doubleSynMeanHeadRadius',' p value:', num2str(p)])
set(gca,'xtick',[1,2],'xticklabel',{'L4', 'L5'})



x = [doubleSynHeadVolume2(:);doubleSynHeadVolume4(:);doubleSynHeadVolume5(:)];
y = [zeros(length(doubleSynHeadVolume2(:)),1);ones(length(doubleSynHeadVolume4(:)),1);ones(length(doubleSynHeadVolume5(:)),1).*2];

figure; violinplot(x,y); title(['doubleSynHeadVolume'])
set(gca,'xtick',[1,2,3],'xticklabel',{'L2', 'L4', 'L5'})


x = [doubleSynHeadVolume2(:);doubleSynHeadVolume4(:);doubleSynHeadVolume5(:)];
y = [zeros(length(doubleSynHeadVolume2(:)),1);ones(length(doubleSynHeadVolume4(:)),1);ones(length(doubleSynHeadVolume5(:)),1).*2];

figure; violinplot(x,y); title(['doubleSynHeadVolume'])
set(gca,'xtick',[1,2],'xticklabel',{'single', 'double'})



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
ratio4 = singleSynMeanHeadRadius4(:)./sqrt(singleSynNeckSection4/pi);
ratio4(isinf(ratio4)) = [];
ratio5 = singleSynMeanHeadRadius5(:)./sqrt(singleSynNeckSection5/pi);
ratio5(isinf(ratio5)) = [];
x = [ratio2;ratio4;ratio5];

y = [zeros(length(ratio2),1);ones(length(ratio4),1);ones(length(ratio5),1).*2];

figure; violinplot(x,y,'ShowMean',true,'ShowMedian',true); title(['ratio between the radius of head and neck'])
set(gca,'xtick',[1,2,3],'xticklabel',{'L2/3', 'L4', 'L5'});

[h,p1,ci,stats] = ttest2(ratio2(:),ratio4(:));
[h,p2,ci,stats] = ttest2(ratio4(:),ratio5(:));
[h,p3,ci,stats] = ttest2(ratio2(:),ratio5(:));
sigstar({[1,2],[2,3], [1,3]},[p1,p2,p3])
ratio2(ratio2 >  quantile(ratio2, 0.99)) = [];
ratio4(ratio4 >  quantile(ratio4, 0.99)) = [];
ratio5(ratio5 >  quantile(ratio5, 0.99)) = [];

figure; histogram(ratio2,200, "Normalization","probability"); 
hold on; histogram(ratio4,200, "Normalization","probability"); 
figure; histogram(ratio4,200, "Normalization","probability"); 
% head radius with neck length
singleSynNeckLength2 = singleSynNeckLength2(id2);
singleSynNeckLength4 = singleSynNeckLength4(id4);
singleSynNeckLength5 = singleSynNeckLength5(id5);
figure; scatter(singleSynMeanHeadRadius2, singleSynNeckLength2*16/1000)
hold on; scatter(singleSynMeanHeadRadius4, singleSynNeckLength4*16/1000)
hold on; scatter(singleSynMeanHeadRadius5, singleSynNeckLength5*16/1000)
%% correlation between spine volume and synaptic cleft size
id2 = find(singleSynHeadVolume2 ~=0);
id4 = find(singleSynHeadVolume4 ~= 0);
id5 = find(singleSynHeadVolume5 ~= 0);
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
p1 = polyfit(x1, y1, 1);
p2 = polyfit(x2, y2, 1);
p3 = polyfit(x3, y3, 1);

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


%% correlation between head volume and contact ratio
id2 = find(singleSynHeadVolume2 ~=0);
id4 = find(singleSynHeadVolume4 ~= 0);
id5 = find(singleSynHeadVolume5 ~= 0);
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

    x1 = singleSynapticCleftSize2;
    x2 = singleSynapticCleftSize4;
    x3 = singleSynapticCleftSize5;
    y1 = singleSynHeadNeckTouchingArea2(:,1) + singleSynHeadNeckTouchingArea2(:,2);
    y2 = singleSynHeadNeckTouchingArea4(:,1) + singleSynHeadNeckTouchingArea4(:,2);
    y3 = singleSynHeadNeckTouchingArea5(:,1) + singleSynHeadNeckTouchingArea5(:,2);
    x1(isnan(y1)) = [];
    y1(isnan(y1)) = [];
    x2(isnan(y2)) = [];
    y2(isnan(y2)) = [];
    x3(isnan(y3)) = [];
    y3(isnan(y3)) = [];

    plot_regression_L2L4L5(x1, x2, x3, y1, y2, y3)

figure;
scatter(singleSynapticCleftSize2, singleSynHeadNeckTouchingArea2(:,1) + singleSynHeadNeckTouchingArea2(:,2))
hold on; scatter(singleSynapticCleftSize4, singleSynHeadNeckTouchingArea4(:,1) + singleSynHeadNeckTouchingArea4(:,2))
hold on; scatter(singleSynapticCleftSize5, singleSynHeadNeckTouchingArea5(:,1) + singleSynHeadNeckTouchingArea5(:,2))
