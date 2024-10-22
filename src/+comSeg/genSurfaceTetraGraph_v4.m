function [face_graph0, uniqueXYZ,tera_x] = genSurfaceTetraGraph_v4(test_region2, xxshift, yyshift, zzshift)
    %automatically avoid one edge from being shared by more than 2 triangles
    %fix the non-manifold vertex in v2 
    % also ensure that no dimension can only have one layer (if only one layer, patch in two directions)
    
    test_region2ID = find(test_region2(:) == 1);
    [lenx, leny,lenz] = size(test_region2);
    [allPointx,allPointy,allPointz] = ind2sub([lenx, leny, lenz],test_region2ID);
    
    [X_tmp,Y_tmp,Z_tmp] = meshgrid(min(allPointx):2:max(allPointx),min(allPointy):2:max(allPointy),min(allPointz):2:max(allPointz));
    id_subsample = sub2ind([lenx, leny, lenz], X_tmp(:), Y_tmp(:), Z_tmp(:));
    id_subsample = intersect(id_subsample, test_region2ID);
    [allPointx2,allPointy2,allPointz2] = ind2sub([lenx, leny, lenz],id_subsample);
    
    allPointxSS = allPointx2(:)' + xxshift(:);
    allPointySS = allPointy2(:)' + yyshift(:);
    allPointzSS = allPointz2(:)' + zzshift(:);
    allPointXYZ = [allPointxSS(:), allPointySS(:), allPointzSS(:)]; % 2, 6, 10 rows contains the newly added integer coordinates
    idRM = find(allPointXYZ(:,1) < 1 | allPointXYZ(:,1) > lenx | allPointXYZ(:,2) < 1| allPointXYZ(:,2) > leny | allPointXYZ(:,3) < 1| allPointXYZ(:,3) > lenz);
    groupx = unique(floor((idRM - 1)/4) + 1);
    idRM_groupID = zeros(4*length(groupx), 1);
    for i = 1:length(groupx)
        idRM_groupID((4*(i-1)+1) : 4*i) = [((groupx(i) - 1)*4 + 1):(4*groupx(i))];
    end
    
    allPointXYZ(idRM_groupID,:) = [];
    
    allPointid = sub2ind([lenx, leny, lenz], allPointXYZ(:,1), allPointXYZ(:,2), allPointXYZ(:,3));
    intergerID = id_subsample;
    % obtain the newly added integer point and then remove them
    
    newlyadded_intid = find(~ismember(allPointid, intergerID));
    newlyadded_intid = intersect(newlyadded_intid, [2:4:length(allPointid)]);
    groupx_newlyadded_intid = unique(floor((newlyadded_intid - 1)/4 + 1));
    idRM_groupID = zeros(4*length(groupx_newlyadded_intid), 1);
    for i = 1:length(groupx_newlyadded_intid)
        idRM_groupID((4*(i-1)+1) : 4*i) = [((groupx_newlyadded_intid(i) - 1)*4 + 1):(4*groupx_newlyadded_intid(i))];
    end
    allPointid(idRM_groupID) = [];
    
    % check the center, remove any volumetric center that are having fewer than
    % 4 surrounding grid points
    
    [uniqueID, iax, icx] = unique(allPointid);
    
    count_uniqueID = accumarray(icx, 1);
    
    idrm_iax = find(~ismember(uniqueID,intergerID) & count_uniqueID <= 4);
    idrm_center_icx = find(ismember(icx, idrm_iax));
    idrm_center_icx = idrm_center_icx(rem(idrm_center_icx,4) ==3);
    groupx_centerid = unique(floor((idrm_center_icx - 1)/4) + 1);
    idRM_groupID = zeros(4*length(groupx_centerid), 1);
    for i = 1:length(groupx_centerid)
        idRM_groupID((4*(i-1)+1) : 4*i) = [((groupx_centerid(i) - 1)*4 + 1):(4*groupx_centerid(i))];
    end
    allPointid(idRM_groupID) = [];
    
    [uniqueid,ia, ic] = unique(allPointid);
    [uniqueX, uniqueY, uniqueZ] = ind2sub([lenx, leny, lenz], uniqueid);
    uniqueXYZ = [uniqueX(:), uniqueY(:), uniqueZ(:)];
    % uniqueXYZ = uniqueXYZ./2;
    tera_x = reshape(ic, 4,[]);
    tera_x = tera_x';
    tera_x = sort(tera_x, 2);
    tera_x = unique(tera_x, 'rows');
    tera_x(tera_x(:,1) == 0| tera_x(:,2) == 0| tera_x(:,3) == 0| tera_x(:,4) ==0,:) = [];
    facesx = [tera_x(:, [1,2,3]), [1:size(tera_x,1)]'; tera_x(:,[1,3,4]), [1:size(tera_x,1)]';tera_x(:, [1,2,4]), [1:size(tera_x,1)]';tera_x(:, [2,3,4]), [1:size(tera_x,1)]'];
    facesx(:,1:3) = sort(facesx(:,1:3), 2);
    facesx = sortrows(facesx, [1:3]);
    [uniqueface, ia_face,ic_face] = unique(facesx(:, 1:3), 'rows');
    countsx = accumarray(ic_face, 1);
    id_single_tri = ia_face(countsx == 1);
    face_graph0 = facesx(id_single_tri, 1:3);% the outer surface of the tetrahedron graph
    %  figure;trisurf(face_graph0,uniqueXYZ(:,1),uniqueXYZ(:,2),uniqueXYZ(:,3),'Facecolor','red','FaceAlpha',0.1);
    end