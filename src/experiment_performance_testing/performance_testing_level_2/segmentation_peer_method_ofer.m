function segmentation_peer_method_ofer(offFolder, ofer_output_folder,volumeGTOutfolder)
% ofer_output_folder = '/work/boyu/EM_astrocyte/dendrite_segmentation_peer_methods/ofer/result';
% offFolder = '/work/boyu/EM_astrocyte/test_segmentation_samples/gt_300/surface_off_300';
% annotationFolder = '/work/boyu/EM_astrocyte/test_segmentation_samples/gt_300/annotation_json_300';
% annotationFolder_gt_curves = '/work/boyu/EM_astrocyte/test_segmentation_samples/gt_300/cut_cycle_ID_300';
listx = dir(fullfile(volumeGTOutfolder, '*.tif'));
% run segmentation 
listx_name = {listx.name};
% for j = 1:length(listx_name)
% 
%     namex = listx_name{j};
%     namex_split = strsplit(namex, '.');
%     namex = namex_split{1};
%     disp(namex)
%     input_file_path = fullfile(offFolder, [namex, '.off']);
%     output_file_path = fullfile(ofer_output_folder, [namex, '.CGAL.txt']);    
%     system(['./ofer_method_CGAL/Segmentation_SDF_Skeleton_Radius',' ',...
%         input_file_path, '>> ', output_file_path]);
% 
% end
% organize the results
for j = 1:length(listx_name)
    namex = listx_name{j};
    namex_split = strsplit(namex, '.');
    namex = namex_split{1};
    disp(namex)

    vol_target = tiffreadVolume(fullfile(volumeGTOutfolder, [namex, '.tif']));
    shaft_region = vol_target == 1;
    [lenx, leny, lenz] = size(vol_target);
    [id_shaft_x, id_shaft_y, id_shaft_z] = ind2sub([lenx, leny, lenz], find(shaft_region(:) == 1));
    id_shaft_x = id_shaft_x(:)*2;
    id_shaft_y = id_shaft_y(:)*2;
    id_shaft_z = id_shaft_z(:)*5;

    [Pts,Tri] = read_off(fullfile(offFolder, [namex,'.off']));
    Tri = Tri';
    Pts = Pts';
    CGAL_output = readtable(fullfile(ofer_output_folder, [namex, '.CGAL.txt']));
    output_mat = table2array(CGAL_output);
    % output_mat = strsplit(output_mat{1},' ');
    face_label = output_mat(2,:);
    face_label = face_label(:);
    face_label = face_label + 1; % the original label is 0, 1 => change to 1, 2
    face_label(isnan(face_label(:))) = 0;
    % find the min-distance from each group to the shaft region
    face_label_idx_cell = label2idx(face_label);
    if(length(face_label_idx_cell) == 1)
        warning('single segmentation ofer method %s', namex)
        else
        if(length(face_label_idx_cell) > 2)
            warning('multiple separation ofer method %s',namex);
            continue;
        end
        face_label_length = cellfun(@length, face_label_idx_cell);
        [~,sortedID] = sort(face_label_length, 'descend');
        face_label_idx_cell = face_label_idx_cell(sortedID(1:2));
        % pick the two largest clusters
        Tri_center = (Pts(Tri(:,1),:) + Pts(Tri(:,2),:) + Pts(Tri(:,3),:))/3;
        face_center_label_1 = Tri_center(face_label_idx_cell{1},:);
        face_center_label_2 = Tri_center(face_label_idx_cell{2},:);
        dist_all_1 = sqrt((id_shaft_x - face_center_label_1(:,1)').^2 + ...
            (id_shaft_y - face_center_label_1(:,2)').^2 + ...
            (id_shaft_z - face_center_label_1(:,3)').^2);
        dist_all_2 = sqrt((id_shaft_x - face_center_label_2(:,1)').^2 + ...
            (id_shaft_y - face_center_label_2(:,2)').^2 + ...
            (id_shaft_z - face_center_label_2(:,3)').^2);
        face_label = face_label.*0;
        if(min(dist_all_1(:)) > min(dist_all_2(:)))
            face_label(face_label_idx_cell{1}) = 2;
            face_label(face_label_idx_cell{2}) = 1;
        else
            face_label(face_label_idx_cell{2}) = 2;
            face_label(face_label_idx_cell{1}) = 1;
        end
        edge_all = [[Tri(:,1), Tri(:,2), [1:size(Tri,1)]']; [Tri(:,1), Tri(:,3), [1:size(Tri,1)]'];[Tri(:,3), Tri(:,2), [1:size(Tri,1)]']];
        edge_all(:,1:2) = sort(edge_all(:,1:2), 2);
        edge_all = sortrows(edge_all,[1,2]);
        edge_all(:,3) = face_label(edge_all(:,3));
        a1 = [1:2:size(edge_all,1)]; % check the label for the triangle on both sides of each edge
        a2 = [2:2:size(edge_all,1)];
        id0 = find(edge_all(a1,3) ~= edge_all(a2,3));
        node_list = edge_all(a1(id0),1:2);
        node_list = node_list(:);
        writematrix(unique(node_list(:)), fullfile(ofer_output_folder,[namex,'.cut.txt']));
        %
        head_neck_label = face_label;
        save(fullfile(ofer_output_folder,[namex,'.face_label.mat']), 'head_neck_label');
    
    
    
    
    %     colormapx = [255, 0, 0, 255;0, 255, 0, 255;0, 0, 255, 255];
        % find the edges that are shared by both 
    
    
    
        % node_list_out
    
    
    %     face_link_can = [edge_all(1:2:end, 3), edge_all(2:2:end, 3)];
    %     face_link_can(output_mat_num(face_link_can(:,1))~=output_mat_num(face_link_can(:,2)),:) = [];
    %     graphx_face = graph(face_link_can(:,1),face_link_can(:,2));
    %     [bins, binsizes] = conncomp(graphx_face);
        % normally the neck part and part of the head part will be wrongly
        % clustered together => use the neck part to decide the right cut cycle
        
    
        % gt_cut = readtable(fullfile(annotationFolder_gt_curves, [namex,'_cut.txt']));
        % gt_cut = table2array(gt_cut);
        % gt_cut_coor = Pts(gt_cut, :);
        % edge_all(:,3) = face_label(edge_all(:,3));
        % % 
        % a1 = [1:2:size(edge_all,1)]; % check the label for the triangle on both sides of each edge
        % a2 = [2:2:size(edge_all,1)];
        % id0 = find(edge_all(a1,3) ~= edge_all(a2,3));
        % node_list = edge_all(a1(id0),1:2);
        % node_list = node_list(:);
        % edge_all_2 = edge_all(ismember(edge_all(:,1), node_list) & ismember(edge_all(:,2), node_list),1:2);
        % graph_edge = graph(edge_all_2(:,1), edge_all_2(:,2));
        % [bins, binsizes] = conncomp(graph_edge);
        % [binsizes_sorted, sortedID] = sort(binsizes, 'descend');
        % sortedID(binsizes_sorted < 5) = [];
        % % since this method might generate multiple cut cycles, we will only output the one that is separated between by the two largest clusters
        
        % min_length = 10000000;
        % selected_id = 1;
        % if(length(sortedID) ~= 1)
        %     for m = 1:length(sortedID)
        %         vertx_tmp = find(bins == sortedID(m));
        %         vertx_tmp_coor = Pts(vertx_tmp,:);
        %         dist_all = sqrt((gt_cut_coor(:,1) - vertx_tmp_coor(:,1)').^2 + ...
        %             (gt_cut_coor(:,2) - vertx_tmp_coor(:,2)').^2 + ...
        %             (gt_cut_coor(:,3) - vertx_tmp_coor(:,3)').^2);
        %         if(min(dist_all(:)) < min_length)
        %             min_length = min(dist_all(:));
        %             selected_id = sortedID(m);
        %         end
        %     end
    
        % else
        %     selected_id = sortedID(1);
        %     face_label = face_label(:);
        % end
        % node_list_out = find(bins == selected_id);
    end


end
    
