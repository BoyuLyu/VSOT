function segmentation_peer_method_dorkenwald(surfaceOutFolder, volumeGTOutfolder,skelfolder,dorkenwalk_result_folder)
    addpath('../resources/curvatures/')
    addpath('../resources/ImprovedSurfaceSmooth/')
    addpath('../resources/CurvatureEstimation')
    addpath('/home/boyu/Documents/edt_mex/edt_mex/edt_mex')
    addpath('../resources/TAUBIN/TAUBIN/')
    addpath('/home/boyu/Documents/tinymesh_mex/tiny_mesh_mex')
    addpath('/home/boyu/Documents/iso2mesh/')
    addpath('../resources/data_to_off/')
    addpath('../resources/data_to_off/')
    addpath('/home/boyu/Documents/graph_related/graph_mex/')
    xxshift = zeros(3,3);
    yyshift = zeros(3,3);
    for i = -1:1
        for j = -1:1
            xxshift((i+2), (j+2)) = i;
            yyshift((i+2), (j+2)) = j;
        end
    end


    % surfaceOutFolder ='/work/boyu/EM_astrocyte/test_segmentation_samples/gt_300/surface_off_300/';
    % volumeGTOutfolder = '/work/boyu/EM_astrocyte/test_segmentation_samples/gt_300/volume_300/'; % the ground truth volume with shaft labeled as 1 and spine labeled as 2
    % skelfolder = '/work/boyu/EM_astrocyte/dendrite_segmentation_peer_methods/dorkenwald/skeleton_300';
    % dorkenwalk_result_folder = '/work/boyu/EM_astrocyte/dendrite_segmentation_peer_methods/dorkenwald/result';
    listx = dir(fullfile(volumeGTOutfolder, '*.tif'));
    id_all = zeros(length(listx), 6); %astrocyte ID, ix, iy, iz, neuronID, dendrite spine ID
    xxshift = zeros(3,3);
    yyshift = zeros(3,3);
    for i = -1:1
        for j = -1:1
            xxshift((i+2), (j+2)) = i;
            yyshift((i+2), (j+2)) = j;
        end
    end
    resx = 16;
    resy = 16;
    resz = 40;
    listx_name = {listx.name};
    for i = 1:length(listx_name)
        try
            namex = listx_name{i};
            namex_split = strsplit(namex, '.');
            id_str = namex_split{1};
            disp(id_str);
            vol_target = tiffreadVolume(fullfile(volumeGTOutfolder, [id_str, '.tif']));
            shaft_region = vol_target == 1;
            vol_target = double(vol_target == 2); % the dendrite spine is labeled as 2 and the shaft is labeled as 1
            skel_target = tiffreadVolume(fullfile(skelfolder, [id_str, '_skel.tif']));
            skel_target = double(skel_target >0);
            id_skel = find(skel_target(:) == 1);
    
            [lenx, leny, lenz] = size(vol_target);
            [id_skel_x, id_skel_y, id_skel_z] = ind2sub([lenx, leny, lenz], id_skel);
            label_skel = zeros(lenx, leny, lenz);
            label_skel(id_skel) = [1:length(id_skel)];
            idxy = rem(id_skel, lenx*leny);
            idz = floor(id_skel/(lenx*leny)) + 1;
            idy = floor(idxy/ lenx) + 1;
            idx = rem(idxy, lenx);
            shiftedID = (idx(:) + xxshift(:)') + (idy(:) + yyshift(:)' - 1).*lenx;
            shiftedID = [shiftedID + repmat(idz(:) - 1, [1,9])*lenx*leny, shiftedID + repmat(idz(:) - 2, [1,9])*lenx*leny, shiftedID + repmat(idz(:), [1,9])*lenx*leny];
            shiftedID(shiftedID(:) <= 0 ) = nan;
            shiftedID(shiftedID(:) >= lenx*leny*lenz) = nan;
            growedID = [repmat(id_skel, 27, 1), shiftedID(:)];
            graphx = label_skel(growedID);
        %     graphx_dist = sqrt((id_skel_x(growedID(:,1)) - id_skel_x(growedID(:,2))).^2 + (id_skel_x(growedID(:,1)) - id_skel_x(growedID(:,2))).^2 ...
        %     +(id_skel_x(growedID(:,1)) - id_skel_x(growedID(:,2))).^2);
            graphx(graphx(:,1) == 0 | graphx(:,2) ==0,:) = [];
            graph_sspath = graph(graphx(:,1), graphx(:,2)); % build the graph of the centerline
            d = distances(graph_sspath);
            [ss,tt] = find(d == max(d(:)));
            spath = shortestpath(graph_sspath, ss(1), tt(1));
            dist_2_boundary = zeros(length(id_skel));
            id_vol_target = find(vol_target(:) == 1);
            [id_vol_x, id_vol_y, id_vol_z] = ind2sub([lenx, leny, lenz], id_vol_target);
            vol_target_1d = logical(vol_target(:));
            tmp_dist_1d = edt_mex(vol_target_1d, lenx,leny, lenz, resx,resy,resz);
            tmp_dist = reshape(tmp_dist_1d, size(vol_target));
            skel_dist = tmp_dist(id_skel);
    
            % order the skeleton id by checking which one is closer to the shaft!!!!!
            skel_dist_ordered = skel_dist(spath);
            id_shaft = find(shaft_region(:) == 1);
            [id_shaft_x, id_shaft_y, id_shaft_z] = ind2sub([lenx, leny, lenz], id_shaft);
            s1_2_shaft = min(sqrt((id_skel_x(spath(1)) - id_shaft_x).^2 + (id_skel_y(spath(1)) - id_shaft_y).^2 + (id_skel_z(spath(1)) - id_shaft_z).^2));
            s2_2_shaft = min(sqrt((id_skel_x(spath(end)) - id_shaft_x).^2 + (id_skel_y(spath(end)) - id_shaft_y).^2 + (id_skel_z(spath(end)) - id_shaft_z).^2));
            if(s1_2_shaft > s2_2_shaft)
                spath = spath(end:-1:1);
                skel_dist_ordered = skel_dist_ordered(end:-1:1);
            end
            % if(mean(skel_dist_ordered(1:round(length(skel_dist_ordered)/2))) > mean(skel_dist_ordered(round(length(skel_dist_ordered)/2):end)))
            %     spath = spath(end:-1:1);
            %     skel_dist_ordered = skel_dist_ordered(end:-1:1);
            % end
            id_small_half = [(1:round(length(skel_dist_ordered)/2))];
            id_large_half = [(round(length(skel_dist_ordered)/2)+ 1):length(skel_dist_ordered)];
            anchor1 = id_small_half(find(skel_dist_ordered(id_small_half) == min(skel_dist_ordered(id_small_half)),1));
            db1 = min(skel_dist_ordered(id_small_half));
            anchor2 = id_large_half(find(skel_dist_ordered(id_large_half) == max(skel_dist_ordered(id_large_half)),1));
            db2 = max(skel_dist_ordered(id_large_half));
            for j = anchor2:-1:anchor1
                if(skel_dist(j) < (1/5*db1 + (4/5)*db2))
                    anchor1 = j;
                    break;
                end

            end
            cutPoint = anchor1;
    
            for j = anchor1 : anchor2
    
                if(skel_dist_ordered(j) > ((1/3)*db1 + (2/3)*db2))
                    cutPoint = j;
                    break;
                end
    
            end
    
            head_skeleton = spath(cutPoint:end);
            neck_skeleton = spath(1:(cutPoint - 1));
            [Pts,Tri] = read_off(fullfile(surfaceOutFolder, [id_str,'.off']));
            Tri  =Tri';
            Pts = Pts';
            face_center_idxyz = (Pts(Tri(:,1),:)+ Pts(Tri(:,2),:)+ Pts(Tri(:,3),:))/3;
            % check which face center should belong to the head part which should
            % belong to the neck part
            %surface coordinates are multiplied by [2,2,5]
            face_label = zeros(size(face_center_idxyz,1),1);
            face_center_idxyz = face_center_idxyz;
            for j = 1:size(face_center_idxyz,1)
                dist2head = min(vecnorm(face_center_idxyz(j,:) - [id_skel_x(head_skeleton).*2, id_skel_y(head_skeleton).*2, id_skel_z(head_skeleton).*5], 2,2));
                dist2neck = min(vecnorm(face_center_idxyz(j,:) - [id_skel_x(neck_skeleton).*2, id_skel_y(neck_skeleton).*2, id_skel_z(neck_skeleton).*5], 2,2));
                if(dist2head < dist2neck)
                    face_label(j) = 2;
                else
                    face_label(j) = 1;
                end
            end
            node_colorMap = zeros(length(Pts),4);
            node_colorMap2 = zeros(length(Pts),1);
            colormapx = [255, 0, 0, 255;0, 255, 0, 255;0, 0, 255, 255];
            edge_all = [[Tri(:,1), Tri(:,2), [1:size(Tri,1)]']; [Tri(:,1), Tri(:,3), [1:size(Tri,1)]'];[Tri(:,3), Tri(:,2), [1:size(Tri,1)]']];
            edge_all(:,1:2) = sort(edge_all(:,1:2), 2);
            edge_all = sortrows(edge_all,[1,2]);
            edge_all(:,3) = face_label(edge_all(:,3));
            a1 = [1:2:size(edge_all,1)];
            a2 = [2:2:size(edge_all,1)];
            id0 = find(edge_all(a1,3) ~= edge_all(a2,3));
            node_list = edge_all(a1(id0),1:2);
            node_list = unique(node_list(:));
            writematrix(node_list, fullfile(dorkenwalk_result_folder,[id_str,'.cut.txt']));
            head_neck_label = face_label;
            save(fullfile(dorkenwalk_result_folder,[id_str,'.face_label.mat']), 'head_neck_label');
        catch ME
            warning('error happened at %s', id_str);
            continue;
        end



    end

end