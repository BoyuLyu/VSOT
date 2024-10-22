function outlabel = classification_dendrite_spine(spineSaveFolder, spine_quantification_folder, spine_classification_folder, resx, resy, resz)
    % based on the features of dendrite spine neck length/ width to classify the dendrite spine
    % input: rootfolder: the root folder of the data
    % the final labels:
    % filopodia: 1, length > 2um
    % mushroom: 2, widht > 0.6um
    % long thin: 3, length > 1um
    % thin: 4, length:width ratio > 1
    % stubby: 5, length: width ratio <= 1
    % output the label of the dendrite spines 
    matSaveFolder = spine_quantification_folder;
    spineSaveFolder = spineSaveFolder;

    load(fullfile(matSaveFolder, 'summary.mat'));
    singleSynHeadVolumeCell = summaryStructure.ss1;
    singleSynHeadMeanRadiusCell = summaryStructure.ss2;
    singleSynNeckLengthCell = summaryStructure.ss3;
    singleSynNeckSectionCell = summaryStructure.ss4;
    singleSynNeckMeanRadiusCell = summaryStructure.ss5;
    singleSynapticCleftSizeCell = summaryStructure.sa1;
    sinsperimeterRatioCell = summaryStructure.sa2;
    sinsperimeterWeightedWrappingAreaCell = summaryStructure.sa3;
    sinspostSynapseTouchingAreaCell = summaryStructure.sa4;
    sinspostSynapseTouchingRatioCell = summaryStructure.sa5;
    sinspreSynapseTouchingAreaCell = summaryStructure.sa6;
    sinspreSynapseTouchingRatioCell = summaryStructure.sa7;
    singleSynHeadNeckTouchingAreaCell = summaryStructure.sa8;
    singleSynHeadNeckTouchingRatioCell = summaryStructure.sa9;

    doubleSynHeadVolumeCell = summaryStructure.ds1;
    doubleSynMeanHeadRadiusCell = summaryStructure.ds2;
    doubleSynNeckLengthCell = summaryStructure.ds3;
    doubleSynNeckSectionCell = summaryStructure.ds4;
    doubleSynNeckMeanRadiusCell = summaryStructure.ds5;
    doubleSynapticCleftSizeCell = summaryStructure.da1;
    dousperimeterRatioCell = summaryStructure.da2;
    dousperimeterWeightedWrappingAreaCell = summaryStructure.da3;
    douspostSynapseTouchingAreaCell = summaryStructure.da4;
    douspostSynapseTouchingRatioCell = summaryStructure.da5;
    douspreSynapseTouchingAreaCell = summaryStructure.da6;
    douspreSynapseTouchingRatioCell = summaryStructure.da7;
    doubleSynHeadNeckTouchingAreaCell = summaryStructure.da8;
    doubleSynHeadNeckTouchingRatioCell = summaryStructure.da9;

    singleSynNeckLengthCell = singleSynNeckLengthCell(:);
    singleSynNeckMeanRadiusCell = singleSynNeckMeanRadiusCell(:);
    doubleSynNeckLengthCell = doubleSynNeckLengthCell(:);
    doubleSynNeckMeanRadiusCell = doubleSynNeckMeanRadiusCell(:);
    singleSynHeadMeanRadiusCell = singleSynHeadMeanRadiusCell(:);
    doubleSynMeanHeadRadiusCell = doubleSynMeanHeadRadiusCell(:);
    singleSynHeadVolumeCell = singleSynHeadVolumeCell(:);
    doubleSynHeadVolumeCell = doubleSynHeadVolumeCell(:);
    neck_length = [singleSynNeckLengthCell + doubleSynNeckLengthCell];
    neck_mean_radius = [singleSynNeckMeanRadiusCell + doubleSynNeckMeanRadiusCell];
    head_mean_radius = [singleSynHeadMeanRadiusCell + doubleSynMeanHeadRadiusCell];
    head_volume = [singleSynHeadVolumeCell + doubleSynHeadVolumeCell];
    head_mean_radius2 = (head_volume.*resx.*resy.*resz/4*3/pi).^(1/3);
    head_mean_radius(head_mean_radius(:) == 1) = head_mean_radius2(head_mean_radius(:) == 1);
    length_width_ratio = neck_length ./ neck_mean_radius;
    outlabel = zeros(length(length_width_ratio), 1);
    for i = 1:length(length_width_ratio)
        if head_volume(i) > 0
            if((neck_length(i)+ head_mean_radius(i))> 2000 && head_mean_radius(i)/neck_mean_radius(i) < 1.5)
                outlabel(i) = 1; % filopodia
            elseif(head_mean_radius(i) > 300)
                outlabel(i) = 2; % mushroom
            elseif((neck_length(i) + head_mean_radius(i)) > 1000)
                outlabel(i) = 3; %longthin
            elseif(length_width_ratio(i) > 2)
                outlabel(i) = 4; % thin
            elseif(length_width_ratio(i) <= 2)
                outlabel(i) = 5; % stubby
            end

        end
    end

    filoFolder = fullfile(spine_classification_folder, 'filopodia');
    if(exist(filoFolder, 'dir'))
        rmdir(filoFolder, 's');
        mkdir(filoFolder);
    elseif(~exist(filoFolder, 'dir'))
        mkdir(filoFolder);
    end 
    mushFolder = fullfile(spine_classification_folder, 'mushroom');
    if(exist(mushFolder, 'dir'))
        rmdir(mushFolder, 's');
        mkdir(mushFolder);
    elseif(~exist(mushFolder, 'dir'))
        mkdir(mushFolder);
    end
    longthinFolder = fullfile(spine_classification_folder, 'longthin');
    if(exist(longthinFolder, 'dir'))
        rmdir(longthinFolder, 's');
        mkdir(longthinFolder);
    elseif(~exist(longthinFolder, 'dir'))
        mkdir(longthinFolder);
    end
    thinFolder = fullfile(spine_classification_folder, 'thin');
    if(exist(thinFolder, 'dir'))
        rmdir(thinFolder, 's');
        mkdir(thinFolder);
    elseif(~exist(thinFolder, 'dir'))
        mkdir(thinFolder);
    end
    stubbyFolder = fullfile(spine_classification_folder, 'stubby');
    if(exist(stubbyFolder, 'dir'))
        rmdir(stubbyFolder, 's');
        mkdir(stubbyFolder);
    elseif(~exist(stubbyFolder, 'dir'))
        mkdir(stubbyFolder);
    end

    for i = 1:length(outlabel)
        if outlabel(i) == 1
            copyfile(fullfile(spineSaveFolder, [num2str(i), '.tif']), filoFolder);
            % save_volshow_2_jpg(filoFolder, [num2str(i), '.tif']);
        elseif outlabel(i) == 2
            copyfile(fullfile(spineSaveFolder, [num2str(i), '.tif']), mushFolder);
            % save_volshow_2_jpg(mushFolder, [num2str(i), '.tif']);
        elseif outlabel(i) == 3
            copyfile(fullfile(spineSaveFolder, [num2str(i), '.tif']), longthinFolder);
            % save_volshow_2_jpg(longthinFolder, [num2str(i), '.tif']);
        elseif outlabel(i) == 4
            copyfile(fullfile(spineSaveFolder, [num2str(i), '.tif']), thinFolder);
            % save_volshow_2_jpg(thinFolder, [num2str(i), '.tif']);
        elseif outlabel(i) == 5
            copyfile(fullfile(spineSaveFolder, [num2str(i), '.tif']), stubbyFolder);
            % save_volshow_2_jpg(stubbyFolder, [num2str(i), '.tif']);
        end
    end





end
