function quantify_dendrite_spine_score_single_psd(curpsID, lenx, leny, lenz,fullSegfolder_root,spine_save_folder, spine_head_neck_save_folder, matSaveFolder, resx, resy, resz)
    % the quantification result is saved in a mat file
    % the variable contained in the mat file is a struct with the following fields:
    % s1: head volume of the dendrite spine head
    % s2: mean radius of the dendrite spine head
    % s3: neck length of the dendrite spine
    % s4: neck section of the dendrite spine
    % s5: mean radius of the dendrite spine neck
    % s6: size of the synaptic cleft that is in contact with the dendrite spine
    % s7: perimeter ratio of astrocyte
    % s8: weighted wrapping area of astrocyte
    % s9: post-synapse touching area
    % s10: post-synapse touching ratio
    % s11: pre-synapse touching area
    % s12: pre-synapse touching ratio
    % s13: head-neck touching area
    % s14: head-neck touching ratio
    % s15: std of the neck radius 
    spine_mask_coordinate_list = table2array(readtable(fullfile(spine_save_folder,'spine_coordinate.csv')));
    synpatic_double_single_indicator = zeros(size(spine_mask_coordinate_list,1), 1);%s0
    singleSynHeadVolume = zeros(length(spine_mask_coordinate_list),1);%s1
    singleSynMeanHeadRadius = zeros(length(spine_mask_coordinate_list),1);%s2
    singleSynNeckLength = zeros(length(spine_mask_coordinate_list),1);%s3
    singleSynNeckSection = zeros(length(spine_mask_coordinate_list),1);%s4   
    singleSynNeckMeanRadius = zeros(length(spine_mask_coordinate_list),1);%s5
    
    singleSynapticCleftSize = zeros(length(spine_mask_coordinate_list),1);%s6
    sinsperimeterRatio = zeros(length(spine_mask_coordinate_list),1);%s7
    sinsperimeterWeightedWrappingArea = zeros(length(spine_mask_coordinate_list),1);%s8
    sinspostSynapseTouchingArea = zeros(length(spine_mask_coordinate_list),1);%s9
    sinspostSynapseTouchingRatio = zeros(length(spine_mask_coordinate_list),1);%s10
    sinspreSynapseTouchingArea = zeros(length(spine_mask_coordinate_list),1);%s11
    sinspreSynapseTouchingRatio = zeros(length(spine_mask_coordinate_list),1);%s12
    singleSynHeadNeckTouchingArea = zeros(length(spine_mask_coordinate_list),2);%s13
    singleSynHeadNeckTouchingRatio = zeros(length(spine_mask_coordinate_list),2);%s14

    singleSynNeckRadiusStd = zeros(length(spine_mask_coordinate_list),1);%s15



    spineCoordinate = zeros(length(spine_mask_coordinate_list),3);
    folder_group = zeros(length(spine_mask_coordinate_list), 6); % assign the group information to each dendrite spine to avoid the repeated calculation of the same chunk
    %         spine_mask_larger_cell = cell(length(spine_mask_coordinate_list), 1);
    %         seg_mask_all_cell = cell(length(spine_mask_coordinate_list), 1);
    %         astro_mask_cell = cell(length(spine_mask_coordinate_list), 1);
    %         cleft_mask_cell = cell(length(spine_mask_coordinate_list), 1);
    spine_mask_bbx_coordinate = zeros(size(spine_mask_coordinate_list,1), 6);
    for i = 1:size(spine_mask_coordinate_list,1)
        if(exist(fullfile(spine_head_neck_save_folder, [num2str(i),'.tif']),'file'))
            spineMask = tiffreadVolume(fullfile(spine_head_neck_save_folder, [num2str(i),'.tif']));
            [lsx, lsy, lsz] = size(spineMask);
            x0 = spine_mask_coordinate_list(i,1);
            y0 = spine_mask_coordinate_list(i,2);
            z0 = spine_mask_coordinate_list(i,3);
            % extract the region with three times size in each
            % direction of the extracted region
            x0_min = max(x0-lsx, 1);
            x0_max = min(x0 + 2*lsx, 5*lenx);
            y0_min = max(y0 - lsy, 1);
            y0_max = min(y0 + 2*lsy,5*leny);
            z0_min = max(z0 - lsz, 1);
            z0_max = min(z0 + 2*lsz, 5*lenz);
            spine_mask_bbx_coordinate(i,1) = x0_min;
            spine_mask_bbx_coordinate(i,2) = y0_min;
            spine_mask_bbx_coordinate(i,3) = z0_min;
            spine_mask_bbx_coordinate(i,4) = x0_max;
            spine_mask_bbx_coordinate(i,5) = y0_max;
            spine_mask_bbx_coordinate(i,6) = z0_max;
    
            ix_min = floor((x0_min-1)/lenx);
            iy_min = floor((y0_min-1)/leny);
            iz_min = floor((z0_min-1)/lenz);
            ix_max = floor((x0_max-1)/lenx);
            iy_max = floor((y0_max-1)/leny);
            iz_max = floor((z0_max-1)/lenz);
            folder_group(i,1) = max(ix_min,0);
            folder_group(i,2) = max(iy_min,0);
            folder_group(i,3) = max(iz_min,0);
            folder_group(i,4) = min(ix_max,4);
            folder_group(i,5) = min(iy_max,4);
            folder_group(i,6) = min(iz_max,4);
        else
            folder_group(i,1) = nan;
            folder_group(i,2) = nan;
            folder_group(i,3) = nan;
            folder_group(i,4) = nan;
            folder_group(i,5) = nan;
            folder_group(i,6) = nan;               
        end
    end
    [unique_rows, ia, ic] = unique(folder_group, 'rows');
        % any rows in spine_mask_bbx_coordinate that are all zeros will be
        % discarded at last (denoted as nan)
    
    
    %Find the dendrite spines that have the same group of chunks and work
    %on them together
    tic;
        for j = 1:size(unique_rows, 1)
            % disp([j, size(unique_rows, 1)])
            if(~isnan(unique_rows(j,1)))
                astro_mask = false((unique_rows(j,4) - unique_rows(j,1) + 1)*lenx, (unique_rows(j,5) - unique_rows(j,2) + 1)*leny,...
                    (unique_rows(j,6) - unique_rows(j,3) + 1)*lenz);
                cleft_mask = false((unique_rows(j,4) - unique_rows(j,1) + 1)*lenx, (unique_rows(j,5) - unique_rows(j,2) + 1)*leny,...
                    (unique_rows(j,6) - unique_rows(j,3) + 1)*lenz);
                seg_mask_all = zeros((unique_rows(j,4) - unique_rows(j,1) + 1)*lenx, (unique_rows(j,5) - unique_rows(j,2) + 1)*leny,...
                    (unique_rows(j,6) - unique_rows(j,3) + 1)*lenz);
                keys_all = [];
                values_all = [];
                max_label = 0;
                for ix = unique_rows(j,1):unique_rows(j,4)
                    for iy = unique_rows(j,2):unique_rows(j,5)
                        for iz = unique_rows(j,3):unique_rows(j,6)        
                                targetFolder = [fullSegfolder_root,'/',num2str(ix),'/', num2str(iy),'/',num2str(iz)];
                                fid = fopen(fullfile(targetFolder, 'seg_mapping.txt')); % Opening the file
                                raw = fread(fid,inf); % Reading the contents
                                str = char(raw'); % Transforma  tion
                                fclose(fid); % Closing the file
                                datax = jsondecode(str); % Using the jsondecode function to parse JSON from string
                                if(isfield(datax,'x0'))
                                    datax = rmfield(datax, 'x0');
                                end
                                segNewID = struct2array(datax);
                                % % Read .h5 file
                                % 
                                % h5FilePath = fullfile(targetFolder, 'segMaskFull.h5');
                                % if exist(h5FilePath, 'file')
                                %     segMaskFull_tmp = h5read(h5FilePath, '/seg_uint16');
                                % else
                                %     error('H5 file not found: %s', h5FilePath);
                                % end
                                segMaskFull_tmp = tiffreadVolume(fullfile(targetFolder,'segMaskFull.tif'));
                                seg_mask_all(((ix - unique_rows(j,1))*lenx + 1):((ix - unique_rows(j,1) + 1)*lenx),...
                                    ((iy - unique_rows(j,2))*leny + 1):((iy - unique_rows(j,2) + 1)*leny)...
                                    ,((iz - unique_rows(j,3))*lenz + 1):((iz - unique_rows(j,3) + 1)*lenz)) = double(segMaskFull_tmp) + double(segMaskFull_tmp>0).*max_label;
                                segRootIDfields = fieldnames(datax);
                                keys = cell(length(segRootIDfields),1);
                                values = cell(length(segRootIDfields),1);
                                for k = 1:length(segRootIDfields)
                                    keys{k} = datax.(segRootIDfields{k}) + max_label;
                                    values{k} = segRootIDfields{k};
                                end
                                keys_all = [keys_all;keys];
                                values_all = [values_all;values];
                                max_label = double(max(segMaskFull_tmp(:))) + max_label; 
                                cleft_tmp = tiffreadVolume(fullfile(targetFolder,['cleft.tif'])) > 0;
                                % h5FilePath = fullfile(targetFolder, 'cleft.h5');
                                % if exist(h5FilePath, 'file')
                                %     cleft_tmp = h5read(h5FilePath, '/seg_uint16');
                                % else
                                %     error('H5 file not found: %s', h5FilePath);
                                % end
                                cleft_mask(((ix - unique_rows(j,1))*lenx + 1):((ix - unique_rows(j,1) + 1)*lenx),...
                                    ((iy - unique_rows(j,2))*leny + 1):((iy - unique_rows(j,2) + 1)*leny)...
                                    ,((iz - unique_rows(j,3))*lenz + 1):((iz - unique_rows(j,3) + 1)*lenz)) = cleft_tmp;
                                astro_tmp = tiffreadVolume(fullfile(targetFolder,['new_astrocyte_seg.tif'])) > 0;
                                astro_mask(((ix - unique_rows(j,1))*lenx + 1):((ix - unique_rows(j,1) + 1)*lenx),...
                                    ((iy - unique_rows(j,2))*leny + 1):((iy - unique_rows(j,2) + 1)*leny)...
                                    ,((iz - unique_rows(j,3))*lenz + 1):((iz - unique_rows(j,3) + 1)*lenz)) = astro_tmp;
                        end
                    end
                end
                id_selected = find(ic == j);
                datax_big = containers.Map(keys_all,values_all);
                for k = 1:length(id_selected)
                    kk = id_selected(k); % the real index to write to the array
                    x0_min = spine_mask_bbx_coordinate(kk,1);
                    y0_min = spine_mask_bbx_coordinate(kk,2) ;
                    z0_min = spine_mask_bbx_coordinate(kk,3);
                    x0_max = spine_mask_bbx_coordinate(kk,4);
                    y0_max = spine_mask_bbx_coordinate(kk,5);
                    z0_max = spine_mask_bbx_coordinate(kk,6);
                    spine_mask = zeros((x0_max - x0_min + 1), (y0_max - y0_min + 1),(z0_max - z0_min + 1));
                    spineMask_small = tiffreadVolume(fullfile(spine_head_neck_save_folder, [num2str(kk),'.tif']));
                    [lsx, lsy, lsz] = size(spineMask_small);
                    spine_mask((spine_mask_coordinate_list(kk,1) - x0_min + 1):(spine_mask_coordinate_list(kk,1) - x0_min + lsx),...
                        (spine_mask_coordinate_list(kk,2) - y0_min + 1):(spine_mask_coordinate_list(kk,2) - y0_min + lsy),...
                        (spine_mask_coordinate_list(kk,3) - z0_min + 1):(spine_mask_coordinate_list(kk,3) - z0_min + lsz)) = spineMask_small;
                    cleft_mask_small = cleft_mask((x0_min - unique_rows(j,1)*lenx):(x0_max - unique_rows(j,1)*lenx),...
                        (y0_min - unique_rows(j,2)*leny ):(y0_max - unique_rows(j,2)*leny), ...
                        (z0_min - unique_rows(j,3)*lenz):(z0_max - unique_rows(j,3)*lenz));
                    astro_mask_small = astro_mask((x0_min - unique_rows(j,1)*lenx):(x0_max - unique_rows(j,1)*lenx),...
                        (y0_min - unique_rows(j,2)*leny):(y0_max - unique_rows(j,2)*leny), ...
                        (z0_min - unique_rows(j,3)*lenz):(z0_max - unique_rows(j,3)*lenz));
                    seg_mask_small = seg_mask_all((x0_min - unique_rows(j,1)*lenx):(x0_max - unique_rows(j,1)*lenx),...
                        (y0_min - unique_rows(j,2)*leny):(y0_max - unique_rows(j,2)*leny), ...
                        (z0_min - unique_rows(j,3)*lenz):(z0_max - unique_rows(j,3)*lenz));
                    [seg_mask_small, datax_small] = structQuant.re_organize_segMask(seg_mask_small, datax_big);
                    curpsIDNew = ['x',curpsID];
                    [single_features_all] = structQuant.gen_quantification_spine_sub_function_single_psd(spine_mask,cleft_mask_small,astro_mask_small,seg_mask_small,datax_small(curpsIDNew), resx, resy, resz);
                    if(~isempty(single_features_all))
                        synpatic_double_single_indicator(kk) = single_features_all.synpatic_double_single_indicator;%s0
                        singleSynHeadVolume(kk) = single_features_all.headVolumex;%s1
                        singleSynMeanHeadRadius(kk) = single_features_all.headMeanRadiusx;%s2
                        singleSynNeckLength(kk) = single_features_all.neckLengthx;%s3
                        singleSynNeckSection(kk) = single_features_all.neckSectionx;%s4   
                        singleSynNeckMeanRadius(kk) = single_features_all.neckMeanRadiusx;%s5
                        singleSynapticCleftSize(kk) = single_features_all.singleSynapticCleftSize;%s6
                        % sinsperimeterRatio(kk) = single_features_all.sinsperimeterRatio;%s7
                        % sinsperimeterWeightedWrappingArea(kk) = single_features_all.sinsperimeterWeightedWrappingArea;%s8
                        % sinspostSynapseTouchingArea(kk) = single_features_all.sinspostSynapseTouchingArea;%s9
                        % sinspostSynapseTouchingRatio(kk) = single_features_all.sinspostSynapseTouchingRatio;%s10
                        % sinspreSynapseTouchingArea(kk) = single_features_all.sinspreSynapseTouchingArea;%s11
                        % sinspreSynapseTouchingRatio(kk) = single_features_all.sinspreSynapseTouchingRatio;%s12
                        % singleSynHeadNeckTouchingArea(kk,:) = single_features_all.singleSynHeadNeckTouchingArea;%s13
                        % singleSynHeadNeckTouchingRatio(kk,:) = single_features_all.singleSynHeadNeckTouchingRatio;%s14
                        spineCoordinate(kk,:) = single_features_all.spineCoordinate;
                        singleSynNeckRadiusStd(kk) = single_features_all.neckRadiusSTD;%s15
                    end

                end
            end
        end
    toc;
    
            summaryStructure.synpatic_double_single_indicator = synpatic_double_single_indicator;
            summaryStructure.ss1 = singleSynHeadVolume;
            summaryStructure.ss2 = singleSynMeanHeadRadius;
            summaryStructure.ss3 = singleSynNeckLength;
            summaryStructure.ss4 = singleSynNeckSection;
            summaryStructure.ss5 = singleSynNeckMeanRadius;
            summaryStructure.sa1 = singleSynapticCleftSize;
            % summaryStructure.sa2 = sinsperimeterRatio;
            % summaryStructure.sa3 = sinsperimeterWeightedWrappingArea;
            % summaryStructure.sa4 = sinspostSynapseTouchingArea;
            % summaryStructure.sa5 = sinspostSynapseTouchingRatio;
            % summaryStructure.sa6 = sinspreSynapseTouchingArea;
            % summaryStructure.sa7 = sinspreSynapseTouchingRatio;
            % summaryStructure.sa8 = singleSynHeadNeckTouchingArea;
            % summaryStructure.sa9 = singleSynHeadNeckTouchingRatio;
            summaryStructure.spineCoordinate = spineCoordinate;
            summaryStructure.ss15 = singleSynNeckRadiusStd;
            if(~exist(matSaveFolder, 'dir'))
                mkdir(matSaveFolder);
            end
            save(fullfile(matSaveFolder, 'summary_spine.mat'), 'summaryStructure')
    
    
    end