function level_2_preprocessing_remove_branch(fullSegfolder_root, spine_save_folder,...
    branched_spine_save_folder, curpsID,lenx,leny, lenz, resx, resy, resz, second_path_length)

    % second_path_length: the length of the second longest branch in the spine (default 20)
mask_spine = false(5*lenx, 5*leny,5*lenz);
mask_dendrite = false(5*lenx, 5*leny, 5*lenz);


xxshift3D = zeros(3,3,3);
yyshift3D = zeros(3,3,3);
zzshift3D = zeros(3,3,3);
for i = -1:1
    for j =-1:1
        for k = -1:1
            xxshift3D((i+2), (j+2), (k+2)) = i;
            yyshift3D((i+2), (j+2), (k+2)) = j;
            zzshift3D((i+2), (j+2), (k+2)) = k;
        end
    end
end


% check if the object is in contact with dendrite shaft in more than two
% regions
% put together the spine mask                     
for ix = 0:4
    for iy = 0:4
        for iz = 0:4
            fullSegfolder = fullfile(fullSegfolder_root, num2str(ix), num2str(iy),num2str(iz));
            if(exist(fullfile(fullSegfolder,['spine_new_no_soma_revised_',num2str(curpsID),'.tif']),'file') )
                tmp = tiffreadVolume(fullfile(fullSegfolder,['spine_new_no_soma_revised_',num2str(curpsID),'.tif']));
                mask_spine((1 + ix*lenx):(ix+1)*lenx, (1 + iy*leny):(iy+1)*leny, (1 + iz*lenz):(iz+1)*lenz) = tmp == 2;
                mask_dendrite((1 + ix*lenx):(ix+1)*lenx, (1 + iy*leny):(iy+1)*leny, (1 + iz*lenz):(iz+1)*lenz) = tmp == 1;           
            end
        end
    end
end 

mask_spine_roi = bwlabeln(mask_spine);
mask_spine_roi_idx = label2idx(mask_spine_roi);
mask_spine_roi_idx(cellfun(@length, mask_spine_roi_idx) < 500) = []; % this is a very relaxed threshold simply to remove any small noisy structures
% save each spine region into a separate file
if(~exist(spine_save_folder, 'dir'))
    mkdir(spine_save_folder)
end
mask_spine_cell = cell(length(mask_spine_roi_idx),1);
mask_spine_bbx = zeros(length(mask_spine_roi_idx),3);
[lx, ly, lz] = size(mask_spine);
se = strel('sphere', 1);
for i = 1:length(mask_spine_roi_idx)
    id_tmp = mask_spine_roi_idx{i};
    [id_tmpx, id_tmpy, id_tmpz] = ind2sub([lx, ly, lz], id_tmp);
    bbx = [max(min(id_tmpx) - 2, 1), min(max(id_tmpx) +2, lx);max(min(id_tmpy) - 2, 1), min(max(id_tmpy) +2, ly); max(min(id_tmpz) - 2, 1), min(max(id_tmpz) +2, lz)];
    mask_spine_bbx(i,:) = [bbx(1,1), bbx(2,1), bbx(3,1)];
    mask_dendrite_tmp = mask_dendrite(bbx(1,1):bbx(1,2), bbx(2,1):bbx(2,2), bbx(3,1):bbx(3,2));
    mask_dendrite_tmp2 = imopen(mask_dendrite_tmp,se); % conversion from 32x32x40 to 16x16x40 makes some dendritic region one layer larger than the denrite spine
    tmp_spine = 2*(mask_spine(bbx(1,1):bbx(1,2), bbx(2,1):bbx(2,2), bbx(3,1):bbx(3,2)))  + mask_dendrite_tmp2;
    mask_spine_cell{i} = tmp_spine;
end
mask_spine = [];
mask_spine_roi = [];
mask_dendrite = [];
clear mask_spine_roi mask_spine mask_dendrite
rmid = [];
branch_spine_id = [];
for i = 1:length(mask_spine_cell)
    tmp_spine = mask_spine_cell{i};
    [lsx, lsy, lsz] = size(tmp_spine);
    if(lsx < 10 || lsy < 10 || lsz < 10)
        rmid = [rmid;i];
        continue;
    else
        % check if the candidate spine is part of the dendrite shaft (the majority is within 1 layer neighbor of the shaft)
        % or if there are multiple touch point at the dendrite shaft
        tmp_spine_roi = bwlabeln(tmp_spine == 2);
        tmp_spine_idx = label2idx(tmp_spine_roi);
        tmp_spine_idx_length = cellfun(@length, tmp_spine_idx);
        id_selected = find(tmp_spine_idx_length == max(tmp_spine_idx_length),1);
        tmp_spine = tmp_spine_roi == id_selected;
        tmp_spine_dilated = imdilate(tmp_spine , se);
        tmp_dendrite_roi = bwlabeln(mask_spine_cell{i} == 1);
        tmp_dendrite_idx = label2idx(tmp_dendrite_roi);
        indi_touch_spine = cellfun(@(c) max(tmp_spine_dilated(c)), tmp_dendrite_idx);
        ratiox = sum(cellfun(@(c) sum(tmp_spine_dilated(c)), tmp_dendrite_idx))/ sum(tmp_spine(:));
        if(sum(indi_touch_spine) >= 2 || sum(indi_touch_spine) == 0 ||ratiox > 0.1)
            % check if the spine region is in contact with the shaft in more than two parts
            % Also check if the spine region is too small to be considered as a spine
            rmid = [rmid;i];
            continue;
        else
            % check if there are branches in the spine 
            skelx = bwskel(tmp_spine, 'MinBranchLength',20);
%                 skelx = imclose(skelx, se);
            bp_point = bwmorph3(skelx, "branchpoints");
            bp_point = imdilate(bp_point, se);
            skelx2 = (skelx - bp_point) > 0;
            skelx_roi = bwlabeln(skelx2);
            skelx_roi_idx = label2idx(skelx_roi);
            if(sum(skelx2(:)) < 10)
                continue;
            elseif(~isempty(skelx_roi_idx))
                % find the longest path between the point( the point cloest to the dendrite shaft)
                tmp_dendrite = tmp_dendrite_roi == find(indi_touch_spine == 1);
                mask_spine_cell{i} = 2*tmp_spine + tmp_dendrite;
                tmp_dendrite_1d = tmp_dendrite(:);
                tmp_den_spine_1d = tmp_dendrite_1d(:) + tmp_spine(:) > 0;
                dist_2_shaft = imchaferDist3D(tmp_den_spine_1d, tmp_dendrite_1d, lsx, lsy, lsz, resx,resy,resz);
                dist_2_shaft = reshape(dist_2_shaft, size(tmp_spine));
                nodeID = find(skelx(:) == 1);
                nearest_skel_point = nodeID(find(dist_2_shaft(nodeID) == min(dist_2_shaft(nodeID)),1)); % the point cloest to the dendrite shaft
                node_map = zeros(lsx, lsy, lsz);
                node_map(nodeID) = 1:length(nodeID);
                nearest_skel_label = node_map(nearest_skel_point);
                [idX0, idY0, idZ0] = ind2sub([lsx, lsy, lsz], nodeID);
                lsx = double(lsx);
                lsy = double(lsy);
                lsz = double(lsz);
                idX = min(max(idX0(:) + xxshift3D(:)',1),lsx);
                idY = min(max(idY0(:) + yyshift3D(:)',1),lsy);
                idZ = min(max(idZ0(:) + zzshift3D(:)',1),lsz);
                outputID_Nei = sub2ind([lsx, lsy, lsz],idX(:),idY(:), idZ(:));
                graphx = [repmat(nodeID(:), 27, 1), outputID_Nei(:)];
                graphx(isnan(graphx(:,1))|isnan(graphx(:,2)), :) = [];
                graphx_label = node_map(graphx);
                graphx_label(graphx_label(:,1) == 0| graphx_label(:,2) ==0,:) = [];
                graphx_label(graphx_label(:,1) == graphx_label(:,2),:) = [];
                graphx_label = sort(graphx_label , 2);
                graphx_label = unique(graphx_label, 'rows');
                g0 = graph(graphx_label(:,1), graphx_label(:,2));
                dist_2_shaft_skel = distances(g0, nearest_skel_label); % the point on the skeleton that is farthest from the dendrite shaft
                farthest_label = find(dist_2_shaft_skel == max(dist_2_shaft_skel),1);
                sspath = shortestpath(g0, nearest_skel_label, farthest_label);
                % remove any small non-branch segments that are overlapped
                % with the shortest path
                sspath_region = false(size(node_map));
                sspath_region(nodeID(sspath)) = 1;
                overlapped_seg = cellfun(@(c) max(sspath_region(c)), skelx_roi_idx);
                skelx_roi_idx(overlapped_seg == 1) = [];
                len_skelx_seg = cellfun(@length, skelx_roi_idx);
                if(~isempty(len_skelx_seg) && sum(len_skelx_seg > second_path_length) >=1)
                    % if after removing one path, there is still one long path remained.
                    rmid = [rmid;i];
                    branch_spine_id = [branch_spine_id;i];
                    continue;      
                end
            end
        end
    end
end

mask_spine_branched_cell = mask_spine_cell(branch_spine_id);
mask_spine_bbx_branched = mask_spine_bbx(branch_spine_id,:);
for i = 1:length(mask_spine_branched_cell)
    tifwrite(uint8(mask_spine_branched_cell{i}), fullfile(branched_spine_save_folder, num2str(i)));
end
writematrix(mask_spine_bbx_branched, fullfile(branched_spine_save_folder,'spine_coordinate.csv'))



mask_spine_cell(rmid) = [];
mask_spine_bbx(rmid,:) = [];
for i = 1:length(mask_spine_cell)
    tifwrite(uint8(mask_spine_cell{i}), fullfile(spine_save_folder, num2str(i)));
end
writematrix(mask_spine_bbx, fullfile(spine_save_folder,'spine_coordinate.csv'))
