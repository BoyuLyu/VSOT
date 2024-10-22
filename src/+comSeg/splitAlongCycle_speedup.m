function outputSurface = splitAlongCycle_speedup(mF2, cutVex,facep)
    edgeGraph = [mF2(:,1), mF2(:,2), [1:size(mF2, 1)]'; mF2(:,2), mF2(:,3), [1:size(mF2, 1)]';mF2(:,1), mF2(:,3), [1:size(mF2, 1)]'];
    cutVex = cutVex(:);
    cutVex_path = [cutVex, [cutVex(2:end); cutVex(1)]];
    edgeGraph(:,1:2) = sort(edgeGraph(:,1:2),2);
    cutVex_path = sort(cutVex_path, 2);
    edgeGraph(ismember(edgeGraph(:,1:2), cutVex_path, 'rows'),:) = [];
    edgeGraph = sortrows(edgeGraph, [1,2]);
    graphx = graph(edgeGraph(1:2:end-1,3), edgeGraph(2:2:end,3));
    bins = conncomp(graphx);
    bins_count = accumarray(bins(:),1);
    [~, sortedID] = sort(bins_count, 'descend');
    if(length(sortedID) == 1)
        % warning('unseparable surface');
    else
        out_binlabel = unique(bins(facep));
        outputSurface = mF2(bins == out_binlabel,:);
    end
end