function [singleSynHeadVolume,singleSynHeadMeanRadius,singleSynNeckLength, singleSynNeckSection, singleSynNeckMeanRadius] =...
    genDendriteSpineScore_nocleft(spineID, lenx, leny, lenz,spineHNROI, resx, resy, resz)

    % addpath('/home/boyu/Documents/edt_mex/edt_mex/edt_mex')
    xxshift = zeros(3,3);
    yyshift = zeros(3,3);
    for i = -1:1
        for j = -1:1
            xxshift((i+2), (j+2)) = i;
            yyshift((i+2), (j+2)) = j;
        end
    end
    spineHNROIx = spineHNROI;
    % spineHNROIx(spineID) = spineHNROI(spineID);
    [allidx, allidy, allidz] = ind2sub([lenx, leny, lenz], [find(spineHNROIx(:) >0)]);
    bbx = [min(allidx), max(allidx)];
    bby = [min(allidy), max(allidy)];
    bbz = [min(allidz), max(allidz)];
    spineHNROIx = spineHNROIx(bbx(1):bbx(2), bby(1):bby(2), bbz(1):bbz(2));
    spineHNROIx = padarray(spineHNROIx, [1,1,1], 0, 'both');
    spineID = find(spineHNROIx(:) > 1);
    % check if the spine contains only single 3 or single 2
    if(sum(spineHNROIx(:) == 3) == 0 || sum(spineHNROIx(:) == 2) == 0)
        spineHNROI_head = double(spineHNROIx > 1);
        spineHNROI_head_roi = bwlabeln(spineHNROI_head);
        spineHNROI_head_roi_idx = label2idx(spineHNROI_head_roi);
        len_roi = cellfun(@length, spineHNROI_head_roi_idx);
        hid = spineHNROI_head_roi_idx{len_roi == max(len_roi)};
        spineHNROI_neck = zeros(size(spineHNROIx));
        singleSynHeadVolume = length(hid);
        se  = strel('sphere', 1);
        singleSynHeadMeanRadius = (singleSynHeadVolume*resx*resy*resz/4*3/pi)^(1/3);
        % cleftID = find(cleftx(:) == 1);
        [lenx, leny, lenz] = size(spineHNROIx);
        % [cleftIDx, cleftIDy,cleftIDz] = ind2sub([lenx, leny, lenz], cleftID);
        % cleftIDxyz = [cleftIDx(:)*16, cleftIDy(:)*16,cleftIDz(:)*40];
        % [coeff,score_w,latent] = pca(cleftIDxyz);  
        % 
        % [hidx, hidy, hidz] = ind2sub([lenx, leny, lenz], hid);
        % hidxyz = [hidx(:)*16, hidy(:)*16, hidz(:)*40];
        % hidxyz = (hidxyz - mean(cleftIDxyz,1))*coeff;
        % % check along the z direction, group all the points within the range of
        % % the same bin
        spineID = find(spineHNROIx(:) > 0);
        % maxRadiusx = zeros(length(max(round(hidxyz(:,3)))), 1);
        % for i = 0:max(round(hidxyz(:,3))) - 1
        %     tmpx = hidxyz(hidxyz(:,3) >= i & hidxyz(:,3) < (i+10),:);
        %     centerxy = [mean(tmpx(:,1)), mean(tmpx(:,2))];
        %     distAll = sqrt((tmpx(:,1) - centerxy(1)).^2 + (tmpx(:,2) - centerxy(2)).^2);
        %     maxRadiusx(i + 1) = mean(distAll); 
        % end
        % singleSynHeadMeanRadius = max(maxRadiusx);
    else
        spineHNROI_head = double(spineHNROIx == 3);
        spineHNROI_head_roi = bwlabeln(spineHNROI_head);
        spineHNROI_head_roi_idx = label2idx(spineHNROI_head_roi);
        len_roi = cellfun(@length, spineHNROI_head_roi_idx);
        hid = spineHNROI_head_roi_idx{len_roi == max(len_roi)};
        spineHNROI_neck = double(spineHNROIx == 2);
        singleSynHeadVolume = length(hid);
        se  = strel('sphere', 1);
        singleSynHeadMeanRadius = (singleSynHeadVolume*resx*resy*resz/4*3/pi)^(1/3);
        % cleftID = find(cleftx(:) == 1);
        [lenx, leny, lenz] = size(spineHNROIx);
        % [cleftIDx, cleftIDy,cleftIDz] = ind2sub([lenx, leny, lenz], cleftID);
        % cleftIDxyz = [cleftIDx(:)*16, cleftIDy(:)*16,cleftIDz(:)*40];
        % [coeff,score_w,latent] = pca(cleftIDxyz);  
        % 
        % [hidx, hidy, hidz] = ind2sub([lenx, leny, lenz], hid);
        % hidxyz = [hidx(:)*16, hidy(:)*16, hidz(:)*40];
        % hidxyz = (hidxyz - mean(cleftIDxyz,1))*coeff;
        % % check along the z direction, group all the points within the range of
        % % the same bin
        % spineID = find(spineHNROIx(:) > 0);
        % maxRadiusx = zeros(length(max(round(hidxyz(:,3)))), 1);
        % for i = 0:max(round(hidxyz(:,3))) - 1
        %     tmpx = hidxyz(hidxyz(:,3) >= i & hidxyz(:,3) < (i+10),:);
        %     centerxy = [mean(tmpx(:,1)), mean(tmpx(:,2))];
        %     distAll = sqrt((tmpx(:,1) - centerxy(1)).^2 + (tmpx(:,2) - centerxy(2)).^2);
        %     maxRadiusx(i + 1) = mean(distAll); 
        % end
        % singleSynHeadMeanRadius = max(maxRadiusx);
        
        % neckShaftJctid = spineHNROI_head_roi_idx{len_roi == min(len_roi)};
        % neckShaftJctid = sort(neckShaftJctid);
        % id1 = neckShaftJctid(round(length(neckShaftJctid)/2));
        % hid_growed = regionGrow3D(hid, lenx, leny, lenz, xxshift, yyshift);
        % headNeckJctid  = hid_growed(spineHNROI_neck(hid_growed) == 1);
        % headNeckJctid = sort(headNeckJctid);
        % id2 = headNeckJctid(round(length(headNeckJctid)/2));
        % headNeckTotal = spineHNROI_head + spineHNROI_neck;
        % spineHNROI_neckx = spineHNROI_neck;
        % spineHNROI_neckx(headNeckJctid) = 2;
        % spineHNROI_neckx(neckShaftJctid) = 2;
    end

%         tifwrite(uint8(spineHNROIx), '../test_spineHNROI' )
%     tifwrite(uint8(spineHNROI_neckx), '../test_neck' )
    if(sum(spineHNROI_neck(:)) > 100)   
        totalID = find(spineHNROIx(:) > 0);
        headNeckTotal = spineHNROIx > 0;
        mask_dendrite_1D = logical(headNeckTotal(:));
        dist_dendrite_1D = edt_mex(mask_dendrite_1D, lenx, leny, lenz, resx,resy,resz);
        mask_dendrite_1D = [];
        dist_dendrite = reshape(dist_dendrite_1D, lenx, leny, lenz);
        curDendriteDist = dist_dendrite;
        curID = totalID;
        curDendriteDist(curID) = dist_dendrite(curID);
        curDendriteDist2 = max(curDendriteDist(:)) - curDendriteDist;
        nodeMap = zeros(lenx, leny, lenz);
        nodeMap(totalID) = 1:length(totalID);
        nodeMap(1,:,:) = 0;
        nodeMap(size(nodeMap, 1),:,:) = 0;
        nodeMap(:,1,:) = 0;
        nodeMap(:,size(nodeMap, 2),:) = 0; 
        nei26 = regionGrow3D(curID, lenx, leny, lenz, xxshift, yyshift);
        score = sqrt(curDendriteDist2(nei26(:,1)).*curDendriteDist2(nei26(:,2)));
        nodex = [nodeMap(nei26),score];
        nodex(nodex(:,1) == 0 | nodex(:,2) ==0,:) = [];
        G = graph(nodex(:,1), nodex(:,2), nodex(:,3));
        id1 = find(spineHNROIx(:) == 1, 1);
        if(isempty(id1))
            id1 = setdiff(find(spineHNROIx(:) == 3), hid);
        end
        if(~isempty(id1))
            id1 = id1(1);
            id2 = hid(1);
            ssPath = shortestpath(G,nodeMap(id1),nodeMap(id2));
            pathID = curID(ssPath);
            pathID = pathID(spineHNROIx(pathID) == 2);
            if(length(pathID) >=2)
                linex = zeros(lenx, leny, lenz);
                linex(pathID) = 1;
                linex_output = imdilate(linex, se);
        %         tifwrite(uint8(linex_output), '../test_centerline' )
                [pathIDx, pathIDy, pathIDz] = ind2sub([lenx, leny, lenz], pathID);
                dist_pair = sqrt((pathIDx(2:end) - pathIDx(1:end-1)).^2.*resx^2 + (pathIDy(2:end) - pathIDy(1:end-1)).^2.*resy^2+ (pathIDz(2:end) - pathIDz(1:end-1)).^2.*resz^2);
                singleSynNeckLength = sum(dist_pair);
            else
                singleSynNeckLength = 0;
            end
        else
            singleSynNeckLength = 0;
        end

          
    else
        singleSynNeckLength = 0;
    end
    singleSynNeckSection = sum(spineHNROI_neck(:))*resx*resy*resz/ singleSynNeckLength;
    singleSynNeckMeanRadius = sqrt(singleSynNeckSection/ pi);

end
