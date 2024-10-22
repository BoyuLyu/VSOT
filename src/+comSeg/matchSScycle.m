function real_id_point = matchSScycle(Vnieuwx,cutVexfineCell, Vnieuw3, coarse_cycle_path)
    meanDist = inf(length(cutVexfineCell),1);
    cutVexfineCell = cutVexfineCell(:);
    corseNNode = Vnieuw3(coarse_cycle_path,:);
    for i = 1:length(cutVexfineCell)
        findCutVex = Vnieuwx(cutVexfineCell{i},:);
        middlePoint = mean(findCutVex,1);
        disttmp = sqrt((middlePoint(:,1) - corseNNode(:,1)).^2 + (middlePoint(:,2) - corseNNode(:,2)).^2 + (middlePoint(:,3) - corseNNode(:,3)).^2);
        meanDist(i) = mean(disttmp);
    end

    real_id_point = coarse_cycle_path(meanDist == min(meanDist));













end