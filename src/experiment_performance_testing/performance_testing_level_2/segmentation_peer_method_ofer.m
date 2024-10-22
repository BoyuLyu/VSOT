clear
ofer_output_folder = '/work/boyu/EM_astrocyte/dendrite_segmentation_peer_methods/ofer/result';
offFolder = '/work/boyu/EM_astrocyte/test_segmentation_samples/gt_300/surface_off_300';
annotationFolder = '/work/boyu/EM_astrocyte/test_segmentation_samples/gt_300/annotation_json_300';
annotationFolder_gt_curves = '/work/boyu/EM_astrocyte/test_segmentation_samples/gt_300/cut_cycle_ID_300';
listx = dir([annotationFolder, '/*.json']);
% run segmentation 
for j = 1:length(listx)

    namex = listx(j).name;
    namex = namex(1:end-9);
    disp(namex)
    input_file_path = fullfile(offFolder, [namex, '.off']);
    output_file_path = fullfile(ofer_output_folder, [namex, '.txt']);    
    system(['./ofer_method_CGAL/Segmentation_SDF_Skeleton_Radius',' ',...
        input_file_path, '>> ', output_file_path])


end
% organize the results
for j = 1:length(listx)
    namex = listx(j).name;
    namex = namex(1:end-9);
    disp(namex)


    [Pts,Tri] = read_off(fullfile('/work/boyu/EM_astrocyte/test_segmentation_samples/annotation/', [namex,'.off']));
    Tri = Tri';
    Pts = Pts';
    output = readtable(fullfile(ofer_output_folder, [namex, '.txt']));
    output_mat = table2array(output);
    output_mat = strsplit(output_mat{1},' ');
    output_mat_num = cellfun(@(c) str2double(c), output_mat);
    output_mat_num(isnan(output_mat_num)) = 0;

    node_colorMap = zeros(length(Pts),4);
    node_colorMap2 = zeros(length(Pts),1);
%     colormapx = [255, 0, 0, 255;0, 255, 0, 255;0, 0, 255, 255];
    % find the edges that are shared by both 
    edge_all = [[Tri(:,1), Tri(:,2), [1:size(Tri,1)]']; [Tri(:,1), Tri(:,3), [1:size(Tri,1)]'];[Tri(:,3), Tri(:,2), [1:size(Tri,1)]']];
    edge_all(:,1:2) = sort(edge_all(:,1:2), 2);
    edge_all = sortrows(edge_all,[1,2]);
%     face_link_can = [edge_all(1:2:end, 3), edge_all(2:2:end, 3)];
%     face_link_can(output_mat_num(face_link_can(:,1))~=output_mat_num(face_link_can(:,2)),:) = [];
%     graphx_face = graph(face_link_can(:,1),face_link_can(:,2));
%     [bins, binsizes] = conncomp(graphx_face);
    % normally the neck part and part of the head part will be wrongly
    % clustered together => use the neck part to decide the right cut cycle
    

    gt_cut = readtable(fullfile(annotationFolder_gt_curves, [namex,'_cut.txt']));
    gt_cut = table2array(gt_cut);
    gt_cut_coor = Pts(gt_cut, :);
    edge_all(:,3) = output_mat_num(edge_all(:,3));
    % 
    a1 = [1:2:size(edge_all,1)]; % check the label for the triangle on both sides of each edge
    a2 = [2:2:size(edge_all,1)];
    id0 = find(edge_all(a1,3) ~= edge_all(a2,3));
    node_list = edge_all(a1(id0),1:2);
    node_list = node_list(:);
    edge_all_2 = edge_all(ismember(edge_all(:,1), node_list) & ismember(edge_all(:,2), node_list),1:2);
    graph_edge = graph(edge_all_2(:,1), edge_all_2(:,2));
    [bins, binsizes] = conncomp(graph_edge);
    [binsizes_sorted, sortedID] = sort(binsizes, 'descend');
    sortedID(binsizes_sorted < 5) = [];
    min_length = 10000000;
    selected_id = 1;
    if(length(sortedID) ~= 1)
        for m = 1:length(sortedID)
            vertx_tmp = find(bins == sortedID(m));
            vertx_tmp_coor = Pts(vertx_tmp,:);
            dist_all = sqrt((gt_cut_coor(:,1) - vertx_tmp_coor(:,1)').^2 + ...
                (gt_cut_coor(:,2) - vertx_tmp_coor(:,2)').^2 + ...
                (gt_cut_coor(:,3) - vertx_tmp_coor(:,3)').^2);
            if(min(dist_all(:)) < min_length)
                min_length = min(dist_all(:));
                selected_id = sortedID(m);
            end
        end
    else
        selected_id = sortedID(1);
    end
    node_list_out = find(bins == selected_id);

    writematrix(unique(node_list_out(:)), fullfile(ofer_output_folder,[namex,'_cut.txt']));
    



end
    
