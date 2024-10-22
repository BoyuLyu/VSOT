function [out_skel_small_scale, bbx, bbx_position] = get_bbx_from_split_skeleton_large_cluster(skel_x, inputRegion, xxshift3D, yyshift3D, zzshift3D)
    % each cell in the bbx contains the volume data, larger than 0,
    % dendrite, larger than 1, centerline/ skeleton
    % split the skeleton and then pick out each part to split the space
    % connnect the skeleton into graph and then check the degree of
    % each node
    [lenx, leny, lenz] = size(skel_x);
    [outLenx, outLeny, outLenz] = size(inputRegion);
    out_skel_small_scale = false(lenx, leny, lenz);
    nodeID = find(skel_x(:) == 1);
    node_map = zeros(lenx, leny, lenz);
    node_map(nodeID) = 1:length(nodeID);
    [idX0, idY0, idZ0] = ind2sub([lenx, leny, lenz], nodeID);
    lenx = double(lenx);
    leny = double(leny);
    lenz = double(lenz);
    idX = min(max(idX0(:) + xxshift3D(:)',1), lenx);
    idY = min(max(idY0(:) + yyshift3D(:)',1),leny);
    idZ = min(max(idZ0(:) + zzshift3D(:)',1),lenz);
    outputID_Nei = sub2ind([lenx,leny, lenz],idX(:),idY(:), idZ(:));
    graphx = [repmat(nodeID(:), 27, 1), outputID_Nei(:)];
    graphx(isnan(graphx(:,1))|isnan(graphx(:,2)), :) = [];
    graphx_label = node_map(graphx);
    graphx_label(graphx_label(:,1) == 0| graphx_label(:,2) ==0,:) = [];
    graphx_label(graphx_label(:,1) == graphx_label(:,2),:) = [];
    graphx_label = sort(graphx_label , 2);
    graphx_label = unique(graphx_label, 'rows');
    g0 = graph(graphx_label(:,1), graphx_label(:,2));
    degree_all = g0.degree;
    dist_all = distances(g0);
    [bins,binsizes] = conncomp(g0);
    bins_idx = label2idx(bins);
    bins_idx = bins_idx(:);
    % search between each two furthest points and then go along the path
    % between the two points 
    saved_nodes = [];
    bins_long_line = [];
    while(~isempty(bins_idx))
         ccid_label = bins_idx{1};
         bins_idx(1) = [];
         dist_ccid = dist_all(ccid_label, ccid_label);
         [maxPaira, maxPairb] = find(dist_ccid == max(dist_ccid(:)));
         pointA = ccid_label(maxPaira(1));
         pointB = ccid_label(maxPairb(1));
         ssPath = shortestpath(g0, pointA, pointB);
         % the only two possible wrong segments might be at the two
         % terminals check the last branch and first branch and then
         % compare their direction with the direction of the path betwee
         % the two branch points (remained to be checked)
         graphx_label_2 = graphx_label(ismember(graphx_label(:,1), ccid_label) & ismember(graphx_label(:,2), ccid_label),:);
         degree_path = degree_all(ssPath);%!!!!!!!remained to be checked
         id_label_bp = find(degree_all >=3);
         graphx_label_2(ismember(graphx_label_2(:,1), id_label_bp)|ismember(graphx_label_2(:,2),id_label_bp),:) = [];
         if(~isempty(graphx_label_2))
             g1 = graph(graphx_label_2(:,1), graphx_label_2(:,2));
             [bins2,binsizes2] = conncomp(g1);
             bins2_idx = label2idx(bins2);
             for k = 1:length(bins2_idx)
                 id_label_short = bins2_idx{k};
                 len_path = sum(sqrt((idX0(id_label_short(2:end)) - idX0(id_label_short(1:end-1))).^2.*32^2 + ...
                     (idY0(id_label_short(2:end)) - idY0(id_label_short(1:end-1))).^2.*32^2 +...
                     (idZ0(id_label_short(2:end)) - idZ0(id_label_short(1:end-1))).^2.*40^2));
                 if(len_path > 1000 && ~any(ismember(id_label_short, ssPath)))
                     bins_idx = [bins_idx;{[id_label_short(:)]}];
                 end
             end
         end
         bins_long_line = [bins_long_line;{ssPath(:)}];
    end
    bins_long_line = bins_long_line(:);
    out_skel_small_scale(nodeID(cell2mat(bins_long_line))) = 1;
    inputRegion_roi = bwlabeln(inputRegion);
    inputRegion_roi_idx = label2idx(inputRegion_roi);
    out_skel_large_scale = imresize3(out_skel_small_scale, [outLenx, outLeny, outLenz], 'Method','nearest');
    inputRegion_roi_idx(cellfun(@(c) max(out_skel_large_scale(c)), inputRegion_roi_idx) == 0) = [];
    combinedRegion = inputRegion + out_skel_large_scale;
    bbx = cell(length(inputRegion_roi_idx),1);
    bbx_position = zeros(length(inputRegion_roi_idx),6);
    for k = 1:length(inputRegion_roi_idx)
        [idx, idy, idz] = ind2sub([outLenx, outLeny, outLenz],inputRegion_roi_idx{k});
        tmp_combined = false([outLenx, outLeny, outLenz]);
        tmp_combined(inputRegion_roi_idx{k}) = 1;
        tmp_combined = double(tmp_combined).*combinedRegion ;
        tmpx = tmp_combined(min(idx):max(idx), min(idy):max(idy), min(idz):max(idz));
        bbx_position(k,:) = [min(idx),max(idx),min(idy),max(idy),min(idz),max(idz)];
%         tmpx = imresize3(tmpx, [outLenx, outLeny, outLenz], 'Method','nearest');
        bbx{k} = tmpx;
    end
    


% 
%     
% 
%     
% % constrain the length of each branch to be longer than certain threshold
%     id_label_bp = find(degree_all >=3);
%     graphx_label_2 = graphx_label;
%     graphx_label_2(ismember(graphx_label_2(:,1), id_label_bp)|ismember(graphx_label_2(:,2),id_label_bp),:) = [];
%     g1 = graph(graphx_label_2(:,1), graphx_label_2(:,2));
%     [bins,binsizes] = conncomp(g1);
%     id_label_short = find(ismember(bins, find(binsizes < 100)));
%     id_label_short = setdiff(id_label_short, id_label_bp);
%     nodeID2 = nodeID;
%     nodeID2(id_label_short) = [];
%     node_map2 = false(lenx, leny, lenz);
%     node_map2(nodeID2) =1;
%     figure; volshow(imdilate(node_map2, se))




    
















end
