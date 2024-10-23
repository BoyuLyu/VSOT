function [out_skel_small_scale, bbx, bbx_position] = get_bbx_from_split_skeleton(skel_x, inputRegion, xxshift3D, yyshift3D, zzshift3D)
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
                 if(length(id_label_short) > 100)
                     bins_long_line = [bins_long_line;{[id_label_short(:)]}];
                 end
             end
         end
         % bins_long_line = [bins_long_line;{ssPath(:)}];
    end
    bins_long_line = bins_long_line(:);
    % check all the geodesic distance between the points to the centerline to cluster the points
    input_region_small = imresize3(inputRegion, [lenx, leny, lenz], 'Method','nearest');
    input_region_small = imdilate(input_region_small, strel('sphere',3)); % to maintain the connectivity between parts
    tic;
    input_region_small_1d = input_region_small(:);
    output_region_roi = uint8(zeros(lenx, leny, lenz));
    dist_2_centerline = uint32(zeros(lenx, leny, lenz)) + 100000;
    for k = 1:length(bins_long_line)
        line_id = nodeID(bins_long_line{k});
        out_skel_small_scale(line_id) = k;
        tmp_centerline = false(lenx, leny, lenz);
        tmp_centerline(line_id) = 1;
        tmp_centerline_1d = tmp_centerline(:);
        dist_1d = imchaferDist3D(input_region_small_1d, tmp_centerline_1d, lenx, leny, lenz, 32,32,40);
        output_region_roi(uint32(dist_1d) < dist_2_centerline(:)) = k;
        dist_2_centerline = min(dist_2_centerline, reshape(uint32(dist_1d), lenx, leny, lenz));
    end
    toc;

    dist_2_centerline = [];
    bins_long_line = [];
    clear dist_2_centerline bins_long_line
    output_region_roi = output_region_roi.*uint8(input_region_small);
    out_skel_large_scale = imresize3(skel_x, [outLenx, outLeny, outLenz], 'Method','nearest');
    output_region_roi_large_scale = imresize3(output_region_roi, [outLenx, outLeny, outLenz], 'Method','nearest');
    output_region_roi_large_scale = imdilate(output_region_roi_large_scale, strel('sphere', 3));
    output_region_roi_large_scale = output_region_roi_large_scale.*uint8(inputRegion);
    output_region_roi_idx = label2idx(output_region_roi_large_scale);
    combinedRegion = inputRegion + out_skel_large_scale;
    bbx = cell(length(output_region_roi_idx),1);
    bbx_position = zeros(length(output_region_roi_idx),6);
    for k = 1:length(output_region_roi_idx)
        [idx, idy, idz] = ind2sub([outLenx, outLeny, outLenz],output_region_roi_idx{k});
        tmp_combined = false([outLenx, outLeny, outLenz]);
        tmp_combined(output_region_roi_idx{k}) = 1;
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
