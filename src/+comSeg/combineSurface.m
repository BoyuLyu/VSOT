function [outSurface, outVertices] = combineSurface(cutFace, cutVertices, headFace, headVertices, bbNodes)
    

    cutFace_uniqueID = unique(cutFace(:));
    mapping1 = zeros(max(cutFace_uniqueID), 1);
    mapping1(cutFace_uniqueID) = 1:length(cutFace_uniqueID);
    cutVertices = cutVertices(cutFace_uniqueID,:);
    newCutFace = mapping1(cutFace);
    headFace_uniqueID = unique(headFace(:));
    mapping2 = zeros(max(headFace_uniqueID), 1);
    mapping2(headFace_uniqueID) = 1:length(headFace_uniqueID);
    headVertices = headVertices(headFace_uniqueID,:);
    newHeadFace = mapping2(headFace);
    % all the normals of the upper half triangles and the directions of
    % these points towards the cut surface should be within 
%     meanCutVex = mean(cutVertices, 1);
%     
%     headVertices_vector = meanCutVex - (headVertices(newHeadFace(:,1),:) + headVertices(newHeadFace(:,2),:) + headVertices(newHeadFace(:,3),:))/3;
%     normals = cross(headVertices(newHeadFace(:,2),:) - headVertices(newHeadFace(:,1),:), headVertices(newHeadFace(:,3),:) - headVertices(newHeadFace(:,2),:),2);
%     dot_pro = dot(normals, headVertices_vector,2);
%     newHeadFace(dot_pro<0,:) = [newHeadFace(dot_pro<0,3), newHeadFace(dot_pro<0,2),newHeadFace(dot_pro<0,1)];


    ccVertex = [cutVertices;headVertices];
    ccVertex = unique(ccVertex, 'rows');
    [Lia,Locb] = ismember(ccVertex,cutVertices, 'rows');
    mapping1_1 = zeros(size(cutVertices,1),1);
    mapping1_1(Locb(Lia)) = find(Lia);
    newCutFace = mapping1_1(newCutFace);
    [Lia,Locb] = ismember(ccVertex,headVertices, 'rows');
    mapping2_1 = zeros(size(headVertices,1),1);
    mapping2_1(Locb(Lia)) = find(Lia);
    newHeadFace = mapping2_1(newHeadFace);
% use the winding order of the head part to re-organize the cut surface
    bbEdge1 = find(ismember(ccVertex, bbNodes(2:3,:), 'rows')); % pick one edge on the boundary
    bbEdge1 = sort(bbEdge1);
    bbEdge1 = bbEdge1(:)';
    % find the two triangles that shares this edge, one on the head part
    % the other on the cut surface
    edge_matx = [newHeadFace(:,1:2), [1:size(newHeadFace,1)]', zeros(size(newHeadFace,1), 1); newHeadFace(:,2:3), [1:size(newHeadFace,1)]',zeros(size(newHeadFace,1), 1); newHeadFace(:,[3,1]), [1:size(newHeadFace,1)]',zeros(size(newHeadFace,1), 1)];
    edge_matx = [edge_matx; ...
        newCutFace(:,1:2), [1:size(newCutFace,1)]', ones(size(newCutFace,1), 1); newCutFace(:,2:3), [1:size(newCutFace,1)]',ones(size(newCutFace,1), 1); newCutFace(:,[3,1]), [1:size(newCutFace,1)]',ones(size(newCutFace,1), 1)];
    pairID = find(ismember(sort(edge_matx(:,1:2), 2), bbEdge1, 'rows'));
    pair_edge_mat = edge_matx(pairID,:);
    pair_triangle_head = pair_edge_mat(pair_edge_mat(:,4) == 0, 3);
    tmpSurface = [newHeadFace(pair_triangle_head,:); newCutFace];
    edge_x = [tmpSurface(:,1:2), [1:size(tmpSurface,1)]';tmpSurface(:,2:3), [1:size(tmpSurface,1)]';tmpSurface(:,[3,1]), [1:size(tmpSurface,1)]'];
    edge_x(:,1:2) = sort(edge_x(:,1:2), 2);
    edge_x = sortrows(edge_x, [1,2]);
    [uedge, ia_edge, ic_edgex] = unique(edge_x(:,1:2), 'rows');
    count_pair = accumarray(ic_edgex, 1);
    count_pair_id = find(count_pair == 2);
    graph_face = graph(edge_x(ia_edge(count_pair_id), 3), edge_x(ia_edge(count_pair_id) + 1, 3));
    %start  from face 1
    tmpSurface = comSeg.checkWindingOrderWithGraph(tmpSurface, graph_face, 1); % correct the winding order for all the vertices
    outSurface = [newHeadFace; tmpSurface(2:end,:)];
    
%     outSurface = checkWindingOrder(outSurface);
%     calVolfromMesh(outSurface2, ccVertex)
%     surfvolume(ccVertex, outSurface)
%     
%     outSurface = [newCutFace;newHeadFace];
    outVertices = ccVertex;




end