function growedID = regionGrow3D_conn6(id0, lenx, leny, lenz)
    xxshift = [0,1,0,-1];
    yyshift = [1,0,-1,0];
    idxy = rem(id0, lenx*leny);
    idz = floor(id0/(lenx*leny)) + 1;
    idy = floor(idxy/ lenx) + 1;
    idx = rem(idxy, lenx);
    shiftedID = (idx(:) + xxshift(:)') + (idy(:) + yyshift(:)' - 1).*lenx;
    shiftedID = shiftedID + repmat((idz - 1).*lenx.*leny, 1, size(shiftedID, 2));
    shiftedID = [shiftedID, id0(:) + [lenx*leny, -lenx*leny]];
    shiftedID(shiftedID(:) <= 0 ) = nan;
    shiftedID(shiftedID(:) >= lenx*leny*lenz) = nan;
    growedID = [repmat(id0, 6, 1), shiftedID(:)];
    growedID(isnan(growedID(:,2)), :) = [];
    growedID(growedID(:,1) == growedID(:,2),:) = [];
    growedID = sort(growedID, 2);
    growedID = unique(growedID, 'rows');


end
