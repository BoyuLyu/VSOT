function num_control_points = gen_gt_from_controlPoints(offFolder, annotationFolder, annotationFolder_gt_curves)
% offFolder = '/work/boyu/EM_astrocyte/materials_for_paper_VSOT/VSOT_data_paired_with_codes/data/performance_testing_level_2_segmentation/large_dataset_3_annotator_L2345_w_stubby/gt_final_800/off_file';
% annotationFolder = '/work/boyu/EM_astrocyte/materials_for_paper_VSOT/VSOT_data_paired_with_codes/data/performance_testing_level_2_segmentation/large_dataset_3_annotator_L2345_w_stubby/gt_final_800/annotation_annotator_1';
% annotationFolder_gt_curves = '/work/boyu/EM_astrocyte/materials_for_paper_VSOT/VSOT_data_paired_with_codes/data/performance_testing_level_2_segmentation/large_dataset_3_annotator_L2345_w_stubby/gt_final_800/cutID_annotator_1';
% 

listx = dir([annotationFolder, '/*.json']);

num_control_points = zeros(length(listx),1);
for j = 1:length(listx)
    
    namex_all = listx(j).name;
    namex = strsplit(namex_all, '.');
    namex = namex{1};
    % disp(namex)

    if(exist(fullfile(offFolder, [namex , '.off']), 'file'))
        if(exist(fullfile(annotationFolder, [namex , '.json']), 'file') && exist(fullfile(offFolder, [namex , '.mrk.json']), 'file'))
            warning('duplicate files %s',namex);
        end
        [Pts,Tri] = read_off(fullfile(offFolder, [namex , '.off']));
        Tri = Tri';
        Pts = Pts';
        
        
        fname = fullfile(annotationFolder, namex_all); 
        fid = fopen(fname); 
        raw = fread(fid,inf); 
        str = char(raw'); 
        fclose(fid); 
        val = jsondecode(str);
        ccPoints_position = zeros(length(val.markups.controlPoints), 3);
        for i = 1:length(val.markups.controlPoints)
            ccPoints_position(i,:) = val.markups.controlPoints(i).position';
        end
        num_control_points(j) = length(val.markups.controlPoints);
        index_ccPoints = zeros(length(val.markups.controlPoints),1);
        for i = 1:length(index_ccPoints)
            distx = sqrt((Pts(:,1) - ccPoints_position(i,1)).^2 + (Pts(:,2) - ccPoints_position(i,2)).^2 + (Pts(:,3) - ccPoints_position(i,3)).^2);
            index_ccPoints(i) = find(distx == min(distx));
        end
        
        %% find the curve
        % build the graph of vertices with distance as the edge weight
        v_graph = [Tri(:,1), Tri(:,2);Tri(:,1), Tri(:,3);Tri(:,2), Tri(:,3)];
        v_graph = sort(v_graph, 2);
        v_graph = unique(v_graph, 'rows');
        e_weight = sqrt((Pts(v_graph(:,1), 1) - Pts(v_graph(:,2), 1)).^2 + (Pts(v_graph(:,1), 2) - Pts(v_graph(:,2), 2)).^2 ...
            + (Pts(v_graph(:,1), 3) - Pts(v_graph(:,2), 3)).^2);
        G_vertex = graph(v_graph(:,1), v_graph(:,2), e_weight);
        
        index_ccPoints_loop = [index_ccPoints(:);index_ccPoints(1)];
        curve_vertex = [];
        for i = 1:(length(index_ccPoints_loop) - 1) 
            s_path = shortestpath(G_vertex, index_ccPoints_loop(i), index_ccPoints_loop(i+1));
            s_path = s_path(:);
            curve_vertex = [curve_vertex; s_path(1:end-1)];
        end
        writematrix(curve_vertex, fullfile(annotationFolder_gt_curves, [namex,'.cut.txt']))
        %% split the surface 
        edge_graph = [Tri(:,1), Tri(:,2), [1:size(Tri,1)]';Tri(:,1), Tri(:,3),[1:size(Tri,1)]';Tri(:,2), Tri(:,3),[1:size(Tri,1)]'];
        edge_graph(:,1:2) = sort(edge_graph(:,1:2), 2);
        edge_graph = sortrows(edge_graph, [1,2]);
        cut_curve_edges = [curve_vertex(1:end),[curve_vertex(2:end);curve_vertex(1)]];
        cut_curve_edges = sort(cut_curve_edges, 2);
        edge_graph_2 = edge_graph;
        edge_graph_2(ismember(edge_graph_2(:,1:2), cut_curve_edges, 'rows'),:) = [];
        edge_G = graph(edge_graph_2(1:2:end,3), edge_graph_2(2:2:end,3));
        [bins, binsizes] = conncomp(edge_G);
        %% color label the cycle
    
        if(length(binsizes) == 1)
            warning('incomplete separation %s',namex);
        else
            colormapx = [255, 0, 0, 255;255, 215, 0,255;0, 255, 0, 255];
            node_colorMap = repmat(colormapx(2,:),length(Pts),1);
        
        
            node_colorMap(curve_vertex,:) = repmat(colormapx(1,:),length(curve_vertex),1) ;
            fid = fopen([annotationFolder_gt_curves,'/', namex,'.color.off'],'wt');
            if( fid==-1 )
                error('Can''t open the file.');
            end
            
            % header
            fprintf(fid, 'COFF\n');
            fprintf(fid, '%d %d 0\n', size(Pts,1), size(Tri,1));
            % write the points & faces
            fprintf(fid, '%f %f %f %d %d %d %d\n', [Pts, node_colorMap]');
            fprintf(fid, '3 %d %d %d\n', Tri'-1);
            fclose(fid);
        
        end
    end

end