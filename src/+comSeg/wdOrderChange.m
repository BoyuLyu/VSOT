function [indit, newNeiF] = wdOrderChange(refF, targetF)
% refF is the triangle used as reference, targetF is the newly input
% triangle
refEdges = [refF(1), refF(2);refF(2), refF(3);refF(3), refF(1)];
targetFEdges = [targetF(1), targetF(2);targetF(2),targetF(3);targetF(3),targetF(1)];
diff = abs(refEdges(:,1) - targetFEdges(:,1)') + abs(refEdges(:,2) - targetFEdges(:,2)');
if(any(diff(:) == 0))
    indit = 1;
    newNeiF = [targetF(3),targetF(2),targetF(1)];
else
    indit = 0;
    newNeiF = [];
end
end