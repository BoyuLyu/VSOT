function outputID = regionGrowWithBBconstraint(regionIDx, regionIDy, regionIDz, xxshift, yyshift, zzshift, lenx, leny, lenz, refinedBoudary, iters, bbx)

    count = 1;
    while count <= iters
        
        regionIDx = regionIDx + xxshift(:)';
        regionIDy = regionIDy + yyshift(:)';
        regionIDz = regionIDz + zzshift(:)';
        regionIDssXYZ = [regionIDx(:), regionIDy(:), regionIDz(:)];
        regionIDssXYZ = unique(regionIDssXYZ, 'rows');
        regionIDss = sub2ind([lenx, leny, lenz], regionIDssXYZ(:,1), regionIDssXYZ(:,2), regionIDssXYZ(:,3));
        regionIDss = setdiff(regionIDss, refinedBoudary);
        regionIDss(bbx(regionIDss) == 0) = [];
        [regionIDx, regionIDy, regionIDz] = ind2sub([lenx, leny, lenz], regionIDss);
        count = count + 1;
    end
    outputID = regionIDss;
end