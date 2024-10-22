function outF = checkWindingOrder(inputF)

edge_all = [inputF(:,1:2),[1:size(inputF,1)]'; inputF(:,2:3),[1:size(inputF,1)]';inputF(:,[3,1]),[1:size(inputF,1)]'];
edge_all(:,1:2) = sort(edge_all(:,1:2), 2);
edge_all = sortrows(edge_all, [1,2]);
graph_faces = graph(edge_all(1:2:end-1,3), edge_all(2:2:end,3));
bins = conncomp(graph_faces);
for i = 1:max(bins(:))
    ccFaces = find(bins == i, 1);
    indic = zeros(size(inputF,1),1);
    while(~isempty(ccFaces))
        indic(ccFaces) = 1;
        s0 = ccFaces(1);
        ccFaces(1) = [];
        neix = neighbors(graph_faces, s0);
        neix(indic(neix) == 1) = [];
        if(~isempty(neix))
            ccFaces = [ccFaces;neix(:)];
            for j = 1:length(neix)
                [indit, newNeiF] = comSeg.wdOrderChange(inputF(s0, :), inputF(neix(j),:));
                if(indit)
                    inputF(neix(j),:) = newNeiF;
                else
                    continue;
                end


            end
        end
    end
end
outF = inputF;








end