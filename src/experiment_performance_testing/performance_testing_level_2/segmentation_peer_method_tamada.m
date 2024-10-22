clear
offFolder = '/work/boyu/EM_astrocyte/test_segmentation_samples/gt_300/surface_off_300';
annotationFolder = '/work/boyu/EM_astrocyte/test_segmentation_samples/gt_300/annotation_json_300';
annotationFolder_gt_curves = '/work/boyu/EM_astrocyte/test_segmentation_samples/gt_300/cut_cycle_ID_300';
listx = dir([offFolder, '/*.off']);
tamada_result_folder = '/work/boyu/EM_astrocyte/dendrite_segmentation_peer_methods/tamada/results';
for j = 1:length(listx)
    
    namex = listx(j).name;
    namex = namex(1:end-4);
    disp(namex)


    [Pts,Tri] = read_off(fullfile(offFolder, [namex , '.off']));
    Tri = Tri';
    Pts = Pts';

    fid = fopen(fullfile(tamada_result_folder, [namex,'.json'])); % Opening the file
    raw = fread(fid,inf); % Reading the contents
    str = char(raw'); % Transformation
    fclose(fid); % Closing the file
    datax = jsondecode(str); % Using the jsondecode function to parse JSON from string
    coor_cut_tmp = datax;
    cut_vex_to_compare = zeros(size(coor_cut_tmp,1), 1);
    
    for i = 1:size(coor_cut_tmp,1)
        dist = sqrt( (Pts(:,1) - coor_cut_tmp(i,1)).^2 + (Pts(:,2) - coor_cut_tmp(i,2)).^2 + (Pts(:,3) - coor_cut_tmp(i,3)).^2);
        cut_vex_to_compare(i) = find(dist == min(dist),1);
    end
    cut_vex_to_compare2 = unique(cut_vex_to_compare,'stable');
    
    v_graph = [Tri(:,1), Tri(:,2);Tri(:,1), Tri(:,3);Tri(:,2), Tri(:,3)];
    v_graph = sort(v_graph, 2);
    v_graph = unique(v_graph, 'rows');
    e_weight = sqrt((Pts(v_graph(:,1), 1) - Pts(v_graph(:,2), 1)).^2 + (Pts(v_graph(:,1), 2) - Pts(v_graph(:,2), 2)).^2 ...
        + (Pts(v_graph(:,1), 3) - Pts(v_graph(:,2), 3)).^2);
    G_vertex = graph(v_graph(:,1), v_graph(:,2), e_weight);
    
    index_ccPoints_loop = [cut_vex_to_compare2(:);cut_vex_to_compare2(1)];
    curve_vertex_to_compare3 = [cut_vex_to_compare2(1)];
    % the points along the cycle are just sparsely positioned link them
    % into a complete cycle
    rest_of_index = setdiff(cut_vex_to_compare2, curve_vertex_to_compare3);
    while(~isempty(rest_of_index))
       dist = distances(G_vertex, curve_vertex_to_compare3(end), rest_of_index);
       next_index = rest_of_index(dist == min(dist));
       pathx = shortestpath(G_vertex, curve_vertex_to_compare3(end), next_index);
       pathx = pathx(:);
       curve_vertex_to_compare3 = [curve_vertex_to_compare3;pathx(2:end)];
       rest_of_index = setdiff(cut_vex_to_compare2, curve_vertex_to_compare3);
    end
    pathx_fin = shortestpath(G_vertex, curve_vertex_to_compare3(end), cut_vex_to_compare2(1));
    pathx_fin = pathx_fin(:);
    curve_vertex_to_compare3 = [curve_vertex_to_compare3;pathx_fin(2:end-1)];
    curve_vertex_to_compare3 = curve_vertex_to_compare3(:);
    % remove any repetitve segment in the part
    edges_all = [curve_vertex_to_compare3(1:end), [curve_vertex_to_compare3(2:end);curve_vertex_to_compare3(1)]];
    edges_all = sort(edges_all, 2);
    [C,ia,ic] = unique(edges_all, 'rows');
    countx = accumarray(ic, 1);
    edge_2_rm = C(countx == 2,:);
    for kk = 1:size(edge_2_rm,1)
        curve_vertex_to_compare3(edges_all(:,1) == edge_2_rm(kk,1) & edges_all(:,2) == edge_2_rm(kk,2),:) = [];
        edges_all(edges_all(:,1) == edge_2_rm(kk,1) & edges_all(:,2) == edge_2_rm(kk,2),:) = [];
    end
    nodeList = curve_vertex_to_compare3(:,1);

    
    writematrix(nodeList, fullfile(tamada_result_folder,[namex,'_cut.txt']));
end



% node_colorMap2 = zeros(size(Pts,1),1);
% node_colorMap2(nodeList(:)) = 1;
% 
% 
% 
%     figure; trisurf(Tri, Pts(:,1), Pts(:,2), Pts(:,3), 'FaceVertexCData',node_colorMap2,'Facecolor','interp') 