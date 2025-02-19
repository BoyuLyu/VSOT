function gen_surface_cut_gt(offFolder, tifFolder, cutFolder)
    listx = dir([cutFolder, '/*.cut.txt']);
    listx_name = {listx.name};
    for i = 1:length(listx_name)
        namex = listx_name{i};
        namex = strsplit(namex, '.');
        namex = namex{1};
        disp(namex)
        [Pts,Tri] = read_off(fullfile(offFolder, [namex , '.off']));
        Tri = Tri';
        Pts = Pts';
        Tri_center = (Pts(Tri(:,1),:) + Pts(Tri(:,2),:) + Pts(Tri(:,3),:))/3;
        vol_target = tiffreadVolume(fullfile(tifFolder, [namex, '.tif']));
        shaft_region = vol_target == 1;
        [lenx, leny, lenz] = size(vol_target);
        [id_shaft_x, id_shaft_y, id_shaft_z] = ind2sub([lenx, leny, lenz], find(shaft_region(:) == 1));
        id_shaft_x = id_shaft_x(:)*2;
        id_shaft_y = id_shaft_y(:)*2;
        id_shaft_z = id_shaft_z(:)*5;
        % split the surface along the cutvex 
        cutID = readtable(fullfile(cutFolder, [namex, '.cut.txt']));
        cutID = table2array(cutID);
        curve_vertex = cutID(:);
        edge_graph = [Tri(:,1), Tri(:,2), [1:size(Tri,1)]';Tri(:,1), Tri(:,3),[1:size(Tri,1)]';Tri(:,2), Tri(:,3),[1:size(Tri,1)]'];
        edge_graph(:,1:2) = sort(edge_graph(:,1:2), 2);
        edge_graph = sortrows(edge_graph, [1,2]);
        cut_curve_edges = [curve_vertex(1:end),[curve_vertex(2:end);curve_vertex(1)]];
        cut_curve_edges = sort(cut_curve_edges, 2);
        edge_graph_2 = edge_graph;
        edge_graph_2(ismember(edge_graph_2(:,1:2), cut_curve_edges, 'rows'),:) = [];
        edge_G = graph(edge_graph_2(1:2:end,3), edge_graph_2(2:2:end,3));
        [bins, binsizes] = conncomp(edge_G);
        face_label_idx = label2idx(bins);
        face_label = bins;
        if(length(binsizes) == 1)
            warning('incomplete separation %s',namex);
            cut_vertex = [];
            face_fin_label = [];
        else
            [~, sortedID] = sort(cellfun(@length, face_label_idx), 'descend');
            face_label_idx = face_label_idx(sortedID(1:2));
            part1 = face_label_idx{1};
            part2 = face_label_idx{2};
            face_center_label_1 = Tri_center(part1,:);
            face_center_label_2 = Tri_center(part2,:);
            dist_all_1 = sqrt((id_shaft_x - face_center_label_1(:,1)').^2 + ...
                (id_shaft_y - face_center_label_1(:,2)').^2 + ...
                (id_shaft_z - face_center_label_1(:,3)').^2);
            dist_all_2 = sqrt((id_shaft_x - face_center_label_2(:,1)').^2 + ...
                (id_shaft_y - face_center_label_2(:,2)').^2 + ...
                (id_shaft_z - face_center_label_2(:,3)').^2);
            face_label = face_label.*0;
            if(min(dist_all_1(:)) > min(dist_all_2(:)))
                face_label(part1) = 2;
                face_label(part2) = 1;
            else
                face_label(part2) = 2;
                face_label(part1) = 1;
            end
            head_neck_label = face_label;
            save(fullfile(cutFolder,[namex,'.face_label.mat']), 'head_neck_label');
        end



    end







end