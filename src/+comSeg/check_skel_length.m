function len_all = check_skel_length(skel_x_parts_roi_idx, lenx, leny, lenz, resx, resy, resz)
    len_all = zeros(length(skel_x_parts_roi_idx), 1);
    for i = 1:length(skel_x_parts_roi_idx)
        idtmp = skel_x_parts_roi_idx{i};
        tmp_mask = zeros(lenx, leny, lenz);
        tmp_mask(idtmp) = 1;
        [idtmpx, idtmpy, idtmpz] = ind2sub([lenx, leny, lenz], idtmp);
        skel3Dline = tmp_mask(min(idtmpx):max(idtmpx), min(idtmpy):max(idtmpy), min(idtmpz):max(idtmpz));
        skel3Dline = padarray(skel3Dline, [1,1,1], 0, 'both');
        skel3Dline_1d = logical(skel3Dline(:));
        startID1_region_1d = zeros(size(skel3Dline_1d),'logical');
        id0 = find(skel3Dline(:) == 1);
        startID1_region_1d(id0(1)) = 1;
        [lenxtmp, lenytmp, lenztmp] = size(skel3Dline);
        dist0 = imchaferDist3D(skel3Dline_1d,startID1_region_1d, lenxtmp, lenytmp, lenztmp, resx,resy,resz);
        dist0(dist0 == 1e10) = 0;
        id1 = find(dist0(:) == max(dist0(:)), 1);
        startID1_region_1d = zeros(size(skel3Dline_1d),'logical');
        startID1_region_1d(id1) = 1;
        dist1 = imchaferDist3D(skel3Dline_1d,startID1_region_1d,lenxtmp, lenytmp, lenztmp, resx,resy,resz);
        dist1(dist1 == 1e10) = 0;
        len_all(i) = max(dist1(:));



    end
    



















end