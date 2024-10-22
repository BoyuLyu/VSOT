function [pathidxyz, width] = calWidthfromVol(vol, point_shaft_xyz, point_head_xyz, resx, resy, resz)

    xxshift = zeros(3,3);
    yyshift = zeros(3,3);
    for i = -1:1
        for j = -1:1
            xxshift((i+2), (j+2)) = i;
            yyshift((i+2), (j+2)) = j;
        end
    end
    [lenx, leny, lenz] = size(vol);
    point_shaft = sub2ind([lenx, leny, lenz], point_shaft_xyz(1), point_shaft_xyz(2), point_shaft_xyz(3));
    point_head = sub2ind([lenx, leny, lenz], point_head_xyz(1), point_head_xyz(2), point_head_xyz(3));
    curID = find(vol(:) >0);
    mask_spine_1D = logical(vol(:));
    dist_spine_1D = edt_mex(mask_spine_1D, lenx, leny, lenz, resx,resy,resz);
    dist_spine = reshape(dist_spine_1D, lenx, leny, lenz);
    dist_spine2 = max(dist_spine(:)) - dist_spine;
    nodeMap = zeros(lenx, leny, lenz);
    nodeMap(curID) = 1:length(curID);
    nodeMap(1,:,:) = 0;
    nodeMap(size(nodeMap, 1),:,:) = 0;
    nodeMap(:,1,:) = 0;
    nodeMap(:,size(nodeMap, 2),:) = 0; 
    nei26 = utils.regionGrow3D(curID, lenx, leny, lenz, xxshift, yyshift);
    score = sqrt(dist_spine2(nei26(:,1)).*dist_spine2(nei26(:,2)));
    nodex = [nodeMap(nei26),score];
    nodex(nodex(:,1) == 0 | nodex(:,2) ==0,:) = [];
    G = graph(nodex(:,1), nodex(:,2), nodex(:,3));
    ssPath = shortestpath(G,nodeMap(point_shaft),nodeMap(point_head));
    pathID = curID(ssPath);
    pathidxyz = zeros(length(pathID),3);
    [pathidxyz(:,1), pathidxyz(:,2), pathidxyz(:,3)] = ind2sub([lenx, leny, lenz], pathID);
    pathidxyz(:,1) = pathidxyz(:,1)*resx;
    pathidxyz(:,2) = pathidxyz(:,2)*resy;
    pathidxyz(:,3) = pathidxyz(:,3)*resz;
    width = dist_spine(pathID);









end