function outF = checkWindingOrderWithGraph(inputF, graph_faces, id)


ccFaces = id;
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

outF = inputF;

end