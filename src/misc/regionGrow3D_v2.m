function [marker_xyz, growedID] = regionGrow3D_v2(id0, lenx, leny, lenz, xxshift, yyshift)
% marker_xyz marks the whether they are at the same xy space or z space
    idxy = rem(id0, lenx*leny);
    idz = floor(id0/(lenx*leny)) + 1;
    idy = floor(idxy/ lenx) + 1;
    idx = rem(idxy, lenx);
    shiftedID = (idx(:) + xxshift(:)') + (idy(:) + yyshift(:)' - 1).*lenx;
    shiftedID = [shiftedID + repmat(idz(:) - 1, [1,9])*lenx*leny, shiftedID + repmat(idz(:) - 2, [1,9])*lenx*leny, shiftedID + repmat(idz(:), [1,9])*lenx*leny];
    marker_xyz = zeros(length(id0),27);
    marker_xyz(:, 1:18) = 1;
    shiftedID(shiftedID(:) <= 0 ) = nan;
    shiftedID(shiftedID(:) >= lenx*leny*lenz) = nan;
    growedID = [repmat(id0, 27, 1), shiftedID(:)];
    marker_xyz = marker_xyz(:);
    marker_xyz(isnan(growedID(:,2))) = [];
    growedID(isnan(growedID(:,2)), :) = [];
    growedID(growedID(:,1) == growedID(:,2),:) = [];
    growedID = sort(growedID, 2);
    [growedID, ia, ic] = unique(growedID, 'rows');
    marker_xyz = marker_xyz(ia);


end
