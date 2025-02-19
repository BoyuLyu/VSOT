function cutvexCell = searchSSpathCutGraph_v2_speedup(Vnieuw, id_selected, mF2,headpoint,neckpoint, ccHead, ccNeck)
    % start from a vertex
%only calculate the selected vertex
graphVertex = zeros(size(mF2,1)*3, 3);

for j = 1:size(mF2,1)
    graphVertex(3*(j-1)+1,1) = mF2(j,1);
    graphVertex(3*(j-1)+1,2) = mF2(j,2);
    graphVertex(3*(j-1)+2,1) = mF2(j,1);
    graphVertex(3*(j-1)+2,2) = mF2(j,3);
    graphVertex(3*(j-1)+3,1) = mF2(j,2);
    graphVertex(3*(j-1)+3,2) = mF2(j,3);
end
distxVertex = vecnorm(Vnieuw(graphVertex(:,1),:) - Vnieuw(graphVertex(:,2),:),2,2);
graphVertex(:,3) = distxVertex;
cutvexCell = cell(size(Vnieuw,1), 1);
% vexInd = [1:size(Vnieuw,1)];
ns = max(mF2(:)) + 1;
ts = ns + 1;
graphVertex2 = [graphVertex; [repmat(ns, length(ccHead),1), ccHead, zeros(length(ccHead),1)];...
    [repmat(ts, length(ccNeck),1), ccNeck, zeros(length(ccNeck),1)]];
G_v = graph(graphVertex2(:,1), graphVertex2(:,2), graphVertex2(:,3));
ssPath = shortestpath(G_v, ns, ts);
ssPath = ssPath(2:end-1);
% cut through the shortest path
cutPath = [ssPath(1:end-1), ssPath(2:end)];
cutPathx = setdiff(ssPath(:),[ccHead;ccNeck]);
neiNodes = mF2(ismember(mF2(:,1), cutPathx(:))| ismember(mF2(:,2), cutPathx(:))|ismember(mF2(:,3), cutPathx(:)), :);
mF2_t = mF2;
Vnieuw_tmp = Vnieuw;
Vnieuw_tmp = [Vnieuw_tmp;Vnieuw_tmp(ssPath(:),:)]; % the vertex with newly added nodes
newNodes = (size(Vnieuw,1) + 1) : (size(Vnieuw,1) + length(unique(ssPath(:))));
mappingx = zeros(max(ssPath(:)),1);
mappingx(ssPath(:)) = newNodes;
graphSmall = [neiNodes(:,1), neiNodes(:,2);neiNodes(:,3), neiNodes(:,2);neiNodes(:,1), neiNodes(:,3)];
graphSmall = sort(graphSmall, 2);
graphSmall = unique(graphSmall,'rows');
graphSmall(ismember(graphSmall(:,1), cutPath(:))|ismember(graphSmall(:,2),cutPath(:)),:) = [];
G_x = graph(graphSmall(:,1), graphSmall(:,2));
[bins,binsizes] = conncomp(G_x);
[~,idx] = sort(binsizes,'descend');
vexLeft = find(bins == idx(1));
vexRight = find(bins == idx(2));
% link all the vexRight to the new nodes
tface = neiNodes(ismember(neiNodes(:,1), vexRight)| ismember(neiNodes(:,2), vexRight) | ismember(neiNodes(:,3), vexRight), :);
remappedFace = tface;
remappedFace(tface > max(cutPath(:)))  = 0;
remappedFace(tface <= max(cutPath(:)))  = mappingx(tface(tface <= max(cutPath(:))));
tface(remappedFace > 0) = remappedFace(remappedFace > 0);
neiNodes(ismember(neiNodes(:,1), vexRight)| ismember(neiNodes(:,2), vexRight) | ismember(neiNodes(:,3), vexRight), :) = tface;
mF2_t(ismember(mF2(:,1), cutPathx(:))| ismember(mF2(:,2), cutPathx(:))|ismember(mF2(:,3), cutPathx(:)), :) = neiNodes;
% aa = zeros(size(Vnieuw_tmp,1),1);
% aa(cutPathx) = 1;
% aa(5372) = 2;
% figure;trisurf(mF2,Vnieuw_tmp(:,1),Vnieuw_tmp(:,2),Vnieuw_tmp(:,3),'FaceVertexCData',aa,'Facecolor','interp');
graphNode_t = [mF2_t(:,1),mF2_t(:,2);mF2_t(:,3),mF2_t(:,2); mF2_t(:,3),mF2_t(:,1)];
graphNode_t = sort(graphNode_t, 2);
graphNode_t = unique(graphNode_t, 'rows');
distx = vecnorm(Vnieuw_tmp(graphNode_t(:,1),:) - Vnieuw_tmp(graphNode_t(:,2),:), 2, 2);
distx(ismember(graphNode_t(:,1), ccHead(:)) | ismember(graphNode_t(:,2), ccHead(:))) = inf;
distx(ismember(graphNode_t(:,1), ccNeck(:)) | ismember(graphNode_t(:,2), ccNeck(:))) = inf;
G_cut = graph(graphNode_t(:,1), graphNode_t(:,2), distx);
for i = 1:length(cutPathx)
    [ss_cycle] = shortestpath(G_cut, cutPathx(i) ,mappingx(cutPathx(i)));
    ss_cycle(ss_cycle == mappingx(cutPathx(i))) = [];
    [Lia,Locb] = ismember(ss_cycle(ss_cycle >= min(newNodes)),newNodes);
    ss_cycle(ss_cycle >= min(newNodes)) = ssPath(Locb);
    cutvexCell{cutPathx(i)} = ss_cycle(:); 
    % aa(ss_cycle(:)) = 1;

end
% 5372
% newSource = max(graphNode_t(:)) + 1;
% newEdges = [repmat(newSource,length(cutPathx), 1),cutPathx];
% distNew = zeros(length(cutPathx), 1);
% G_cut_new = graph([graphNode_t(:,1);newEdges(:,1)], [graphNode_t(:,2);newEdges(:,2)],[distx(:); distNew(:)] );
id_selected = setdiff(id_selected, cutPathx);
cutvexCell_small = cell(length(id_selected), 1);
headpp = mF2(ismember(mF2(:,1),ccHead) | ismember(mF2(:,2),ccHead) | ismember(mF2(:,3),ccHead) ,:);
headpp = unique(headpp(:));
neckpp = mF2(ismember(mF2(:,1),ccNeck) | ismember(mF2(:,2),ccNeck) | ismember(mF2(:,3),ccNeck) ,:);
neckpp = unique(neckpp(:));
for j = 1:length(id_selected)
    i = id_selected(j);
    try
    if(any(headpp == i)|| any(neckpp == i))
        continue;
    else
        d1 = distances(G_cut,i,ssPath);
        d2 = distances(G_cut, i, newNodes);
        sumx = d1 + d2;
        labelleft = ssPath(sumx == min(sumx));
        sspathtmpLeft = shortestpath(G_cut, labelleft(1), i);
        labelRight = newNodes(sumx == min(sumx));
        sspathtmpRight = shortestpath(G_cut, i, labelRight(1));
        sspathtmpRight = sspathtmpRight(:);
        ss_cycle = [sspathtmpLeft(:);sspathtmpRight(2:end-1)];
        redundant_points= ss_cycle(ss_cycle >= min(newNodes));
        ss_cycle(ss_cycle >= min(newNodes)) = ssPath(redundant_points - min(newNodes) + 1);       
        cutvexCell_small{j} = ss_cycle;
    end
    
    catch ME
     error('problem happens at %s', num2str(i))
    continue;
    end
end
cutvexCell(id_selected) = cutvexCell_small;
end
        
                
                
                
  