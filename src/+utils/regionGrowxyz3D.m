function outputID = regionGrowxyz3D(inputID,iters, lx,ly,lz)
   %%% input the numbers of the circles to grow
    [idX, idY, idZ] = ind2sub([lx,ly, lz],inputID);
    xxshift = zeros(2*iters+1, 2*iters+1);
    yyshift = zeros(2*iters+1, 2*iters+1);
    zzshift = zeros(2*iters+1, 2*iters+1);
    for i = -iters:iters
        for j = -iters:iters
            for k = -iters:iters
                xxshift(i+ iters+1, j + iters +1, k +iters + 1) = i;
                
                yyshift(i+ iters+1, j + iters +1, k +iters + 1) = j;
                zzshift(i+ iters+1, j + iters +1, k +iters + 1) = k;
            end
        end
    end
    
    idX = min(max(idX + xxshift(:)',1), lx);
    idY = min(max(idY + yyshift(:)',1),ly);
    idZ = min(max(idZ + zzshift(:)',1),lz);
    outputID = sub2ind([lx,ly, lz],idX(:),idY(:), idZ(:));
    outputID = unique(outputID);
    [idX, idY] = ind2sub([lx, ly,lz], outputID);

    
end