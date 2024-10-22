function [newCurvature] = genCurvatureFace4Scale_speedup(mF2, Vnieuw,distrange,idSelected)
% group all the surface points within the range of the target pixel
% build a graph of vertex
graphNode = [[mF2(:,1), mF2(:,2)];[mF2(:,2), mF2(:,3)]; [mF2(:,1), mF2(:,3)]];
graphNode = sort(graphNode,2);
graphNode = unique(graphNode,'rows');
G0 = graph(graphNode(:,1), graphNode(:,2));
newCurvature  = zeros(size(Vnieuw,1),1) + nan;
% [Vnieuw2, index] = sortrows(Vnieuw, [1,2,3]);
newCurvature_sub = zeros(length(idSelected),1);
Vnieuw2 = Vnieuw(idSelected, :);
for i = 1:length(idSelected)
  % disp(i)
    curDist = 0;
    distAll = vecnorm(Vnieuw - Vnieuw2(i,:),2, 2);
    curNode = idSelected(i);
    neiNew = [];
    while curDist < distrange      
        neiNode = [mF2(ismember(mF2(:,1), curNode),:);mF2(ismember(mF2(:,2), curNode),:);mF2(ismember(mF2(:,3), curNode),:)];
        neiNew = setdiff(neiNode(:), curNode);
        curDist = max(distAll(neiNew));
        curNode = unique(neiNode(:));
    end
    %redo triangulation using the vertices
    if(~isempty(neiNew))
        H = subgraph(G0, neiNew);
        cyclesCell = cyclebasis(H);
        lengthCycle = cellfun(@length, cyclesCell);
        cycleSelected = cyclesCell{find(lengthCycle == max(lengthCycle))};
        cycleSelected = [cycleSelected(:);cycleSelected(1)];
        % uniformly split into 6 groups
        if(length(cycleSelected) > 6)
            splitChunks = zeros(6,1);
            splitChunks = splitChunks + floor(length(cycleSelected)/6);
            remainx = rem(length(cycleSelected),6);
            splitChunks(1:remainx) = splitChunks(1:remainx) + 1;
        else
            splitChunks = ones(length(cycleSelected));
        end
        splitChunks = cumsum(splitChunks);
        binNodex = neiNew([cycleSelected(1);cycleSelected(splitChunks(:))]);
        mFtmp = zeros(min(6, length(cycleSelected)), 3);
        for j = 1:size(mFtmp,1)
            mFtmp(j,:) = [idSelected(i),binNodex(j),binNodex(j+1)];
        end
        % figure;trisurf(mFtmp,Vnieuw(:,1),Vnieuw(:,2),Vnieuw(:,3),'FaceVertexCData',xxc,'Facecolor','interp');
        FV2 = [];
        [unqiueNode,ia, ic] = unique(mFtmp(:));
        mFtmp = reshape(ic, [], 3);
        FV2.faces = mFtmp;
        FV2.vertices = Vnieuw(unqiueNode,:);
        getderivatives = 0;
        [PrincipalCurvatures2,PrincipalDir1,PrincipalDir2,FaceCMatrix,VertexCMatrix,Cmagnitude, vertexNormals]= GetCurvatures(FV2 ,getderivatives);
        newCurvature_sub(i) = PrincipalCurvatures2(1,1)*PrincipalCurvatures2(2,1);
    end


end

newCurvature(idSelected) = newCurvature_sub;



end