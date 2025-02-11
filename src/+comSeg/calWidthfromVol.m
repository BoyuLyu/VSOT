function [pathidxyz, width] = calWidthfromVol(vol, point_shaft_xyz, point_head_xyz, point_non_neck, resx,resy, resz)
%   point_non_neck is with the actual scale
%   point_shaft_xyz, point_head_xyz are under the volume raw scale

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
    % mask_spine_1D = logical(vol(:));
    % dist_spine_1D = edt_mex(mask_spine_1D, lenx, leny, lenz, resx,resy,resz);
    % dist_spine = reshape(dist_spine_1D, lenx, leny, lenz);
    % dist_spine2 = max(dist_spine(:)) - dist_spine;
    dist_spine = zeros(lenx, leny, lenz);
    for i = 1:lenz
        dist_spine(:,:,i) = bwdist(1 - vol(:,:,i));

    end
    nodeMap = zeros(lenx, leny, lenz);
    nodeMap(curID) = 1:length(curID);
    nodeMap(1,:,:) = 0;
    nodeMap(size(nodeMap, 1),:,:) = 0;
    nodeMap(:,1,:) = 0;
    nodeMap(:,size(nodeMap, 2),:) = 0; 
    nei26 = utils.regionGrow3D(curID, lenx, leny, lenz, xxshift, yyshift);
    score = sqrt(dist_spine(nei26(:,1)).*dist_spine(nei26(:,2)));
    nodex = [nodeMap(nei26),1./score];
    nodex(nodex(:,1) == 0 | nodex(:,2) ==0,:) = [];
    G = graph(nodex(:,1), nodex(:,2), nodex(:,3));
    ssPath = shortestpath(G,nodeMap(point_shaft),nodeMap(point_head));
    pathID = curID(ssPath);
    pathidxyz = zeros(length(pathID),3);
    [pathidxyz(:,1), pathidxyz(:,2), pathidxyz(:,3)] = ind2sub([lenx, leny, lenz], pathID);
    pathidxyz(:,1) = pathidxyz(:,1)*resx;
    pathidxyz(:,2) = pathidxyz(:,2)*resy;
    pathidxyz(:,3) = pathidxyz(:,3)*resz;
    width = zeros(size(pathidxyz,1),1);
    for i = 1:size(pathidxyz,1)
        distall = vecnorm((point_non_neck - pathidxyz(i,:)),2,2);
        width(i) = min(distall);
    end









end