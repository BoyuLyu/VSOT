function dendrite_feature  = quantify_dendrite_sub_func(curpsID, lenx, leny, lenz,resx, resy, resz, fullSegfolder_root, spine_save_folder, neuron_ds_folder)
    mask_spine = false(5*lenx, 5*leny,5*lenz);
    mask_dendrite = false(5*lenx, 5*leny, 5*lenz);
    % spine_save_folder = fullfile(fullSegfolder_root, curpsID);
    neuron_x_removal = tiffreadVolume(fullfile(neuron_ds_folder, [num2str(curpsID),'_dendrite_soma.tif'])) > 0;
    neuron_x_removal_roi = bwlabeln(neuron_x_removal);
    neuron_x_removal_roi_idx = label2idx(neuron_x_removal_roi);
    neuron_x_removal_roi_idx(cellfun(@length, neuron_x_removal_roi_idx) < 1000) = [];
    neuron_x_removal_roi_idx = neuron_x_removal_roi_idx(:);
    neuron_x_removal = false(size(neuron_x_removal));
    neuron_x_removal(cell2mat(neuron_x_removal_roi_idx)) = 1;
    skel_x = tiffreadVolume(fullfile(neuron_ds_folder, [num2str(curpsID),'_skel.tif']));
    skel_x = skel_x > 0;
    se = strel('sphere', 1);
    neuron_x_removal = imdilate(neuron_x_removal, se);
    % split the skeleton into sub-branches and extract each sub branch
    skel_x = skel_x.*(1 - neuron_x_removal) > 0;
    skel_x_branchpoints = bwmorph3(skel_x, 'branchpoints');
    skel_x_parts = skel_x - imdilate(skel_x_branchpoints, se) >0;
    skel_x_parts_roi = bwlabeln(skel_x_parts);
    skel_x_parts_roi_idx = label2idx(skel_x_parts_roi);
    skel_x_parts_roi = skel_x_parts_roi.*0;
    skel_x_parts_roi_idx(cellfun(@length, skel_x_parts_roi_idx) < 20) = [];
    for i = 1:length(skel_x_parts_roi_idx)
        skel_x_parts_roi(skel_x_parts_roi_idx{i}) = i;
    end
    % the length of each skeleton branch(80x80x200)
    len_skel_parts = comSeg.check_skel_length(skel_x_parts_roi_idx, lenx, leny, lenz, 5*resx,5*resy,5*resz);    
    for ix = 0:4
        for iy = 0:4
            for iz = 0:4
                fullSegfolder = fullfile(fullSegfolder_root, num2str(ix), num2str(iy),num2str(iz));
                if(exist(fullfile(fullSegfolder,['spine_new_no_soma_revised_',num2str(curpsID),'.tif']),'file') )
                    tmp = tiffreadVolume(fullfile(fullSegfolder,['spine_new_no_soma_revised_',num2str(curpsID),'.tif']));
                    mask_spine((1 + ix*lenx):(ix+1)*lenx, (1 + iy*leny):(iy+1)*leny, (1 + iz*lenz):(iz+1)*lenz) = tmp == 2;
                    mask_dendrite((1 + ix*lenx):(ix+1)*lenx, (1 + iy*leny):(iy+1)*leny, (1 + iz*lenz):(iz+1)*lenz) = tmp == 1;
    %                 tifwrite(uint8(tmp),fullfile(fullSegfolder,['spine_new_no_soma_revised_',num2str(curpsID)]))                        
                elseif(exist(fullfile(fullSegfolder,['spine_new_no_soma_',num2str(curpsID),'.tif']),'file'))
                    tmp = tiffreadVolume(fullfile(fullSegfolder,['spine_new_no_soma_',num2str(curpsID),'.tif']));
                    mask_spine((1 + ix*lenx):(ix+1)*lenx, (1 + iy*leny):(iy+1)*leny, (1 + iz*lenz):(iz+1)*lenz) = tmp == 2;
                    mask_dendrite((1 + ix*lenx):(ix+1)*lenx, (1 + iy*leny):(iy+1)*leny, (1 + iz*lenz):(iz+1)*lenz) = tmp == 1;
    %               
                end
            end
        end
    end 
    mask_spine_ds_all = imresize3(mask_spine, [lenx, leny, lenz], 'Method', 'nearest'); % 80x80x200
    mask_dendrite_ds_all = imresize3(mask_dendrite, [lenx, leny, lenz], 'Method', 'nearest') + mask_spine_ds_all; % 80x80x200
    mask_dendrite_ds = mask_dendrite_ds_all - neuron_x_removal > 0;
    mask_dendrite_ds_1d = mask_dendrite_ds(:);
    mask_dendrite_ds_rm_edt = edt_mex(mask_dendrite_ds_1d,lenx,leny, lenz, 5*resx,5*resy,5*resz);

    mask_dendrite_dt_2D = zeros(lenx, leny, lenz);
    for i = 1:lenz
        mask_dendrite_dt_2D(:,:,i) = bwdist(1 - mask_dendrite_ds(:,:,i));
    end
    %refine the skeleton to force it closer to the centerline 
    curDendriteDist2 = max(mask_dendrite_dt_2D(:)) - mask_dendrite_dt_2D;
    nodeMap = zeros(lenx, leny, lenz);
    totalID = find(mask_dendrite_ds(:) > 0);
    nodeMap(totalID) = 1:length(totalID);
    nodeMap(1,:,:) = 0;
    nodeMap(size(nodeMap, 1),:,:) = 0;
    nodeMap(:,1,:) = 0;
    nodeMap(:,size(nodeMap, 2),:) = 0; 
    nei26 = regionGrow3D(totalID, double(lenx), double(leny), double(lenz), xxshift, yyshift);
    score = sqrt(curDendriteDist2(nei26(:,1)).*curDendriteDist2(nei26(:,2)));
    nodex = [nodeMap(nei26),score];
    nodex(nodex(:,1) == 0 | nodex(:,2) ==0,:) = [];
    G = graph(nodex(:,1), nodex(:,2), nodex(:,3));
    skel_x_parts_roi_idx_dist_center = skel_x_parts_roi_idx;
    for i = 1:length(skel_x_parts_roi_idx)
        curID = skel_x_parts_roi_idx{i};
        curID(nodeMap(curID) ==0) = [];
        if(~isempty(curID))
            [curIDx, curIDy, curIDz] = ind2sub([lenx, leny, lenz], curID);
            dist_matrix = squareform(pdist([curIDx(:), curIDy(:),curIDz(:)]));
            [max_dist, idx] = max(dist_matrix(:));
            [furthest_point1, furthest_point2] = ind2sub(size(dist_matrix), idx);
            furthest_points = [curID(furthest_point1), curID(furthest_point2)];
            ssPath = shortestpath(G,nodeMap(furthest_points(1)),nodeMap(furthest_points(2)));
            pathID = totalID(ssPath);
            skel_x_parts_roi_idx_dist_center{i} = pathID;
        end
    end
    % thickness_skel_parts = cellfun(@(c) mean(mask_dendrite_dt_2D(c)), skel_x_parts_roi_idx_dist_center)*80;



    thickness_skel_parts = cellfun(@(c) mean(mask_dendrite_ds_rm_edt(c)), skel_x_parts_roi_idx_dist_center)*5*lenx; 
    len_skel_parts = comSeg.check_skel_length(skel_x_parts_roi_idx_dist_center, lenx, leny, lenz, 5*resx,5*resy,5*resz); 

    % create a score map so that the center of the dendrite shaft is the lowest. Wateshed algorithm grows from the center to the boundary.
    mask_dendrite_ds_rm_edt = - mask_dendrite_ds_rm_edt;
    mask_dendrite_ds_rm_edt = imgaussfilt3(mask_dendrite_ds_rm_edt, 1);
    mask_dendrite_ds_rm_edt(skel_x_parts(:) > 0) = -inf;
    mask_dendrite_ds_dist = reshape(mask_dendrite_ds_rm_edt, size(mask_dendrite_ds));
    L = watershed(mask_dendrite_ds_dist);
    L(~mask_dendrite_ds_all) = 0;
    L_idx = label2idx(L);
    L_idx = L_idx(:);
    L_idx(cellfun(@length,L_idx) == 0) = [];
    final_mask = skel_x_parts_roi;
    se2 = strel('sphere',2);
    while ~isempty(L_idx)
    % check each segments from the watershed

        ws_with_skel = cellfun(@(c) max(final_mask(c)), L_idx);
        if(max(ws_with_skel) == 0)
            break;
        else
            for i = 1:length(skel_x_parts_roi_idx)
                final_mask(cell2mat(L_idx(ws_with_skel == i))) = i;
            end
            final_mask = imdilate(final_mask, se2).*double(mask_dendrite_ds);
            L_idx(ws_with_skel > 0) = [];
        end

    end
    spine_coordinates = table2array(readtable(fullfile(spine_save_folder, 'spine_coordinate.csv')));
    dendrite_feature = zeros(size(spine_coordinates,1), 4); 
    % store the saved spine ID, saved dendrite ID, the length of the corresponding dendrite, the thickness of the corresponding dendrite
    mask_spine_roi = bwlabeln(mask_spine);
    mask_spine_roi_idx = label2idx(mask_spine_roi);
    final_mask_original_scale = imresize3(uint8(final_mask), [5*lenx, 5*leny, 5*lenz], 'Method', 'nearest');
    final_mask_original_scale2 = final_mask_original_scale;
    for i = 1:size(spine_coordinates,1)
        spine_small = tiffreadVolume(fullfile(spine_save_folder, [num2str(i), '.tif']));
        [lenx_small, leny_small, lenz_small] = size(spine_small);
        spine_labels = mask_spine_roi(spine_coordinates(i,1):(spine_coordinates(i,1) + lenx_small - 1)...
            , spine_coordinates(i, 2): (spine_coordinates(i, 2) + leny_small - 1), ...
            spine_coordinates(i, 3):(spine_coordinates(i, 3)+lenz_small - 1));
        spine_labels(spine_labels(:) == 0) = [];
        spine_label = mode(spine_labels);
        if(spine_label ~= 0 && ~isnan(spine_label))
            dendrite_branch_label = max(final_mask_original_scale(mask_spine_roi_idx{spine_label}));
            final_mask_original_scale2(mask_spine_roi_idx{spine_label}) = dendrite_branch_label;
            dendrite_feature(i,1) = i;
            dendrite_feature(i,2) = dendrite_branch_label;
            if(dendrite_branch_label ~= 0)
                dendrite_feature(i,3) = len_skel_parts(dendrite_branch_label);
                dendrite_feature(i,4) = thickness_skel_parts(dendrite_branch_label);
            end
        end
    end
