function growedID = regionGrow3D(id0, lenx, leny, lenz, xxshift, yyshift)

    idxy = rem(id0, lenx*leny);
    idz = floor(id0/(lenx*leny)) + 1;
    idy = floor(idxy/ lenx) + 1;
    idx = rem(idxy, lenx);
    shiftedID = (idx(:) + xxshift(:)') + (idy(:) + yyshift(:)' - 1).*lenx;
    shiftedID = [shiftedID + repmat(idz(:) - 1, [1,9])*lenx*leny, shiftedID + repmat(idz(:) - 2, [1,9])*lenx*leny, shiftedID + repmat(idz(:), [1,9])*lenx*leny];
    shiftedID(shiftedID(:) <= 0 ) = nan;
    shiftedID(shiftedID(:) >= lenx*leny*lenz) = nan;
    growedID = [repmat(id0, 27, 1), shiftedID(:)];
    growedID(isnan(growedID(:,2)), :) = [];
    growedID(growedID(:,1) == growedID(:,2),:) = [];
    growedID = sort(growedID, 2);
    growedID = unique(growedID, 'rows');


end
