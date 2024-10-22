function maskRemoveAll = extract_shaft_v5(mask_dendrite,skel_x, dist_dendrite, xxshift, yyshift, resx, resy, resz)
%one parameter to decide the size of the dendrite shaft
    [lenx, leny, lenz] = size(mask_dendrite);
    maskRemoveAll = false(lenx,leny,lenz);
    maskRemoveAll(:,:,1) = mask_dendrite(:,:,1);
    maskRemoveAll(:,:,lenz) = mask_dendrite(:,:,lenz);
    maskRemoveAll(:,leny,:) = mask_dendrite(:,leny,:);
    maskRemoveAll(:,1,:) = mask_dendrite(:,1,:);
    maskRemoveAll(lenx,:,:) = mask_dendrite(lenx,:,:);
    maskRemoveAll(1,:,:) = mask_dendrite(1,:,:);
    mask_dendrite(mask_dendrite > 0) = 1;
    se = strel('sphere',3);
    skel_x = imdilate(skel_x, se);
    skel_x = skel_x.*mask_dendrite > 0;
    skel_x_id = skel_x(:);
    pathID = find(skel_x(:) == 1);
    distLine = dist_dendrite(pathID);
    radius = min(distLine(distLine > resx));
    mask_target_whole_Area1d = logical(mask_dendrite(:));
    dist_2_center = imchaferDist3D(mask_target_whole_Area1d,skel_x_id, lenx, leny, lenz, resx,resy,resz);
    dist_2_center = reshape(dist_2_center, size(mask_dendrite));
    maskRemove = dist_2_center <= radius;
    se = strel('sphere',1);
    boundary_mask = xor(imdilate(mask_dendrite,se) , mask_dendrite);
    boundary_mask = imdilate(boundary_mask, se) & mask_dendrite;
    graphMapping = zeros(lenx, leny, lenz);
    ccID = find(mask_dendrite(:) == 1);
    graphMapping(ccID) = 1:length(ccID);
    source_cc = graphMapping(maskRemove(:));
    sink_cc = graphMapping(boundary_mask(:));
    sink_cc = setdiff(sink_cc, source_cc);
    rest_cc = setdiff([1:length(ccID)], [source_cc(:);sink_cc(:)]);
    % construct the graph with only six connectivity
    nei6 = utils.regionGrow3D_conn6(ccID, lenx, leny, lenz);
    mask_target_whole_Area1d = logical(mask_dendrite(:));
    maskRemove1d = logical(maskRemove(:));
    scoreIJ = imchaferDist3D(mask_target_whole_Area1d,maskRemove1d, lenx, leny, lenz,resx,resy,resz);
    scoreIJ(scoreIJ == 1e10) = inf;
    scoreIJ = reshape(scoreIJ, size(mask_dendrite));
    score_graph_1 =sqrt(scoreIJ(nei6(:,1)).*scoreIJ(nei6(:,2)));
    score_graph_1(isnan(score_graph_1)) = 0;
    scale_factor = 8;
    score_graph_1 = score_graph_1./scale_factor;
    graph_mat_1 = graphMapping(nei6);
    psSourceID = length(ccID) + 1;
    psSinkID = length(ccID) + 2;
    graph_mat_2 = [[source_cc(:);rest_cc(:)],repmat(psSourceID, length(source_cc) + length(rest_cc), 1)]; % link to source, this relates to the surface area of the final extracted shaft region
    score_graph_2 = [inf(length(source_cc), 1); ones(length(rest_cc),1)];
    lambda = 1/4;
    score_graph_2 = score_graph_2*(resx*resy*resz)^(2/3)/scale_factor*lambda;
    graph_mat_3 = [sink_cc(:), repmat(psSinkID, length(sink_cc), 1)]; % link to sink
    score_graph_3 = inf(length(sink_cc), 1);
    graph_matx = [graph_mat_1, score_graph_1;graph_mat_2,score_graph_2; graph_mat_3,  score_graph_3];
    graph_matx(graph_matx(:,1) == 0| graph_matx(:,2) == 0,:) = [];
    graph_matx((ismember(graph_matx(:,1), source_cc) & ismember(graph_matx(:,2), sink_cc))...
     | (ismember(graph_matx(:,2), source_cc) & ismember(graph_matx(:,1), sink_cc)),:) = [];
    if(~isempty(graph_matx))
        G2 = graph(graph_matx(:,1), graph_matx(:,2), graph_matx(:,3));
        [mf, Gf, cs, ct] = maxflow(G2, psSourceID, psSinkID);
        cs(cs == psSourceID) = [];
        maskRemoveAll(ccID(cs)) = 1;
        maskRemove21d = logical(maskRemoveAll(:));
        dist_new = imchaferDist3D(mask_target_whole_Area1d,maskRemove21d, lenx, leny, lenz, resx,resy,resz);
        dist_new = reshape(dist_new,size(mask_dendrite));
        maskRemoveAll = dist_new <= (resx + resy);
    end
    maskRemoveAll_roi = bwlabeln(maskRemoveAll);
    maskRemoveAll_roi_idx = label2idx(maskRemoveAll_roi);
    len_roi = cellfun(@length, maskRemoveAll_roi_idx);
    maskRemoveAll_roi_idx(len_roi<20000) = [];
    maskRemoveAll_roi_idx = maskRemoveAll_roi_idx(:);
    maskRemoveAll = false(size(maskRemoveAll));
    maskRemoveAll(cell2mat(maskRemoveAll_roi_idx)) = 1;
end