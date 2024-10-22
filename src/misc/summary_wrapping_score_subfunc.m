function layer2_summary = summary_wrapping_score_subfunc(nf_list,rootFolder,neuronListFolder)
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
ttc = 0;
for nf = nf_list

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
        else
                continue;
        end

        se = strel('sphere', 2);        
        ttc = ttc + length(summaryStructure.ss1);
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

disp(ttc)
end
