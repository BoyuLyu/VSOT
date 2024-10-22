clear
offFolder = '/work/boyu/EM_astrocyte/test_segmentation_samples/gt_300/surface_off_300';
annotationFolder = '/work/boyu/EM_astrocyte/test_segmentation_samples/gt_300/annotation_json_300';
annotationFolder_gt_curves = '/work/boyu/EM_astrocyte/test_segmentation_samples/gt_300/cut_cycle_ID_300';
listx = dir([offFolder, '/*.off']);
our_method_coordinate_folder = '/work/boyu/EM_astrocyte/test_segmentation_samples/our_segmentation_result_cut_coordinates';
our_method_cut_result_folder = '/work/boyu/EM_astrocyte/dendrite_segmentation_peer_methods/our_method/result';
for j = 1:length(listx)
    
    namex = listx(j).name;
    namex = namex(1:end-4);
    disp(namex)


    [Pts,Tri] = read_off(fullfile(offFolder, [namex , '.off']));
    Tri = Tri';
    Pts = Pts';

% 
% figure;trisurf(Tri,Pts(:,1),Pts(:,2),Pts(:,3),'Facecolor','red','FaceAlpha',0.1);
% 
% figure;trisurf(mF2,Vnieuwx(:,1),Vnieuwx(:,2),Vnieuwx(:,3),'FaceVertexCData',ccScore_fin,'Facecolor','interp');
% hold on;trisurf(mF2,Vnieuwx(:,1),Vnieuwx(:,2),Vnieuwx(:,3), 'Facecolor','blue','FaceAlpha',0.1)

%% find a cycle in the mesh for the annotation that is close to the cycle in the tetrahedra mesh
% for each vertex in the tetra-mesh ,find its closest vertex on the smooth
% mesh
    if(exist(fullfile(our_method_coordinate_folder,[namex,'.mrk.json.mat']), 'file'))
        aa = load(fullfile(our_method_coordinate_folder,[namex,'.mrk.json.mat']));
        coor_cut_tmp = aa.coor_cut_tmp;
        cut_vex_to_compare = zeros(size(coor_cut_tmp,1), 1);
        
        for i = 1:size(coor_cut_tmp,1)
            dist = sqrt( (Pts(:,1) - coor_cut_tmp(i,1)).^2 + (Pts(:,2) - coor_cut_tmp(i,2)).^2 + (Pts(:,3) - coor_cut_tmp(i,3)).^2);
            cut_vex_to_compare(i) = find(dist == min(dist),1);
        end
        cut_vex_to_compare2 = unique(cut_vex_to_compare,'stable');
    % cut_vex_coor = Pts(cut_vex_to_compare2(:),:);
    % graph_cut_vex = zeros(length(cut_vex_to_compare2),3);
    % for i = 1:size(cut_vex_to_compare2,1)
    %     dist_coor_diff = cut_vex_coor - cut_vex_coor(i,:); 
    %     dist = sqrt(dist_coor_diff(:,1).^2 + dist_coor_diff(:,2).^2 + dist_coor_diff(:,3).^2);
    %     [~,sorted_id] = sort(dist);
    %     graph_cut_vex(i,1) = sorted_id(1);
    %     graph_cut_vex(i,2) = sorted_id(2);
    %     graph_cut_vex(i,3) = sorted_id(3);
    % end
    % graph_cut_vex_G = graph([graph_cut_vex(:,1);graph_cut_vex(:,1)],[graph_cut_vex(:,2);graph_cut_vex(:,3)]); 
    % cut_vex_to_compare2
    
        v_graph = [Tri(:,1), Tri(:,2);Tri(:,1), Tri(:,3);Tri(:,2), Tri(:,3)];
        v_graph = sort(v_graph, 2);
        v_graph = unique(v_graph, 'rows');
        e_weight = sqrt((Pts(v_graph(:,1), 1) - Pts(v_graph(:,2), 1)).^2 + (Pts(v_graph(:,1), 2) - Pts(v_graph(:,2), 2)).^2 ...
            + (Pts(v_graph(:,1), 3) - Pts(v_graph(:,2), 3)).^2);
        G_vertex = graph(v_graph(:,1), v_graph(:,2), e_weight);
        
        index_ccPoints_loop = [cut_vex_to_compare2(:);cut_vex_to_compare2(1)];
        curve_vertex_to_compare3 = [];
        for i = 1:(length(index_ccPoints_loop) - 1) 
            s_path = shortestpath(G_vertex, index_ccPoints_loop(i), index_ccPoints_loop(i+1));
            s_path = s_path(:);
            curve_vertex_to_compare3 = [curve_vertex_to_compare3; s_path(1:end-1)];
        end
        writematrix(curve_vertex_to_compare3(:), fullfile(our_method_cut_result_folder,[namex,'_cut.txt']));
    end
    
    % 
%     node_colorMap2 = zeros(size(Pts,1),1);
%     node_colorMap2(cut_vex_to_compare2(:)) = 1;
% figure; trisurf(Tri, Pts(:,1), Pts(:,2), Pts(:,3), 'FaceVertexCData',node_colorMap2,'Facecolor','interp') 

end
