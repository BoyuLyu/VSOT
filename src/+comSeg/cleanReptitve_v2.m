function newCell = cleanReptitve_v2(inputCell)
    lengthAll = cellfun(@length, inputCell);
    lengthAll2 = cellfun(@(c) length(unique(c)),inputCell );
    xxx = lengthAll - lengthAll2;
    ccID =find(xxx~= 0);
    newCell = inputCell;
    for i =1:length(ccID)
        tmpNode = newCell{ccID(i)};
        [u_label, ia, ic] = unique(tmpNode);
        a_counts = accumarray(ic,1);
        repeti_ID = u_label(a_counts == 2);
        rm_idx  = [];
        for j = 1:length(repeti_ID)
            idx = find(tmpNode ==repeti_ID(j));
            if(idx(2) - idx(1) > 0.5*length(tmpNode))
                rm_idx = [rm_idx; [1:idx(1)]'; [idx(2)+1:length(tmpNode)]'];
            else
                rm_idx = [rm_idx;[idx(1)+1:idx(2)]'];
            end
        end
        rm_idx = unique(rm_idx);
        tmpNode(rm_idx) = [];
        newCell{ccID(i)} = tmpNode;
    end

end