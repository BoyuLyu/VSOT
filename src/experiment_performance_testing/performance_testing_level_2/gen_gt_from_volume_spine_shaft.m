function gen_gt_from_volume_spine_shaft(offFolder, tifFolder, output_cut_folder)
    % offFolder = '/work/boyu/EM_astrocyte/materials_for_paper_VSOT/test_segmentation_samples/dendrite_spine_head_neck/gt_batch_2/stubby_samples';
    % annotationFolder_gt_curves = '/work/boyu/EM_astrocyte/materials_for_paper_VSOT/test_segmentation_samples/dendrite_spine_head_neck/gt_batch_2/cut_cycle_ID';
    % directly output the spine-shaft intersetion as the boundary

    files = dir(fullfile(offFolder, '*.off'));
    filesname = {files.name};
    
    for j = 1:length(filesname)
        id_3 = strsplit(filesname{j}, {'_', '.'});
        n = str2double(id_3{1});
        m = str2double(id_3{2});
        k = str2double(id_3{3});
        
        name_file = [num2str(n), '_', num2str(m), '_', num2str(k)];
    
        labelx = tiffreadVolume(fullfile(tifFolder, [name_file, '.tif']));
        labelxidx = label2idx(labelx);
        [lenx, leny, lenz] = size(labelx);
        [SpineCoorx,SpineCoory, SpineCoorz]  = ind2sub([lenx, leny, lenz], labelxidx{2});
        [ShaftCoorx,ShaftCoory, ShaftCoorz]  = ind2sub([lenx, leny, lenz], labelxidx{1});
        test_region = labelx >0;
        [node,elem,face,regions]=vol2surf(double(test_region),1:size(test_region,1),1:size(test_region,2),1:size(test_region,3),1.01,1,'cgalsurf');
        Tri = elem(:,1:3);
        Pts = node;
        Pts(:,1) = Pts(:,1) .*2;
        Pts(:,2) = Pts(:,2) .*2;
        Pts(:,3) = Pts(:,3) .*5;
        [Pts] = taubinsmooth( Tri,[Pts(:,1),Pts(:,2),Pts(:,3)],10);
        Tri_center = [mean([Pts(Tri(:,1), 1),Pts(Tri(:,2), 1),Pts(Tri(:,3), 1)],2), ...
            mean([Pts(Tri(:,1), 2),Pts(Tri(:,2), 2),Pts(Tri(:,3), 2)],2),...
            mean([Pts(Tri(:,1), 3),Pts(Tri(:,2), 3),Pts(Tri(:,3), 3)],2)];
        head_neck_label = zeros(length(Tri_center),1);
    
        for i = 1:size(Tri_center,1)
            distHead = sqrt((Tri_center(i,1) - SpineCoorx*2).^2 + (Tri_center(i,2) - SpineCoory*2).^2 + (Tri_center(i,3) - SpineCoorz*5).^2);
            distNeck = sqrt((Tri_center(i,1) - ShaftCoorx*2).^2 + (Tri_center(i,2) - ShaftCoory*2).^2 + (Tri_center(i,3) - ShaftCoorz*5).^2);
            if(min(distHead) > min(distNeck))
                head_neck_label(i) = 1;
            else
                head_neck_label(i) = 2;
            end
        end
        edge_all = [[Tri(:,1), Tri(:,2), [1:size(Tri,1)]']; [Tri(:,1), Tri(:,3), [1:size(Tri,1)]'];[Tri(:,3), Tri(:,2), [1:size(Tri,1)]']];
        edge_all(:,1:2) = sort(edge_all(:,1:2), 2);
        edge_all = sortrows(edge_all,[1,2]);
        edge_all(:,3) = head_neck_label(edge_all(:,3));
        a1 = [1:2:size(edge_all,1)];
        a2 = [2:2:size(edge_all,1)];
        id0 = find(edge_all(a1,3) ~= edge_all(a2,3));
        edge_all = edge_all(a1(id0),:);
        edge_all = edge_all(:,1:2);
        edge_all = sort(edge_all, 2);
        edge_all = unique(edge_all, 'rows');
        % re-order the vertices to make it a loop
        graph_loop = graph(edge_all(:,1), edge_all(:,2));
        allCycles = allcycles(graph_loop);
        cycle_length = cellfun(@length, allCycles);
        [~, idx] = max(cycle_length);
        curve_vertex = allCycles{idx};

        % the surface file of the spine only 
        [Pts2,Tri2] = read_off(fullfile(offFolder, [name_file , '.off']));
        Tri2 = Tri2';
        Pts2 = Pts2';
        Tri_center = [mean([Pts2(Tri2(:,1), 1),Pts2(Tri2(:,2), 1),Pts2(Tri2(:,3), 1)],2), ...
        mean([Pts2(Tri2(:,1), 2),Pts2(Tri2(:,2), 2),Pts2(Tri2(:,3), 2)],2),...
        mean([Pts2(Tri2(:,1), 3),Pts2(Tri2(:,2), 3),Pts2(Tri2(:,3), 3)],2)];

        % find the curve close to the separatation curve between the dendritic spine and the shaft
        curve_vertex_2 = zeros(length(curve_vertex),1);
        for i = 1:length(curve_vertex)
            distx = sqrt((Pts2(:,1) - Pts(curve_vertex(i),1)).^2 + (Pts2(:,2) - Pts(curve_vertex(i),2)).^2 + (Pts2(:,3) - Pts(curve_vertex(i),3)).^2);
            curve_vertex_2(i) = find(distx == min(distx));
        end
        curve_vertex_2 = unique(curve_vertex_2,"stable");
        % connect these points into a complete loop
        v_graph = [Tri2(:,1), Tri2(:,2);Tri2(:,1), Tri2(:,3);Tri2(:,2), Tri2(:,3)];
        v_graph = sort(v_graph, 2);
        v_graph = unique(v_graph, 'rows');
        e_weight = sqrt((Pts2(v_graph(:,1), 1) - Pts2(v_graph(:,2), 1)).^2 + (Pts2(v_graph(:,1), 2) - Pts2(v_graph(:,2), 2)).^2 ...
            + (Pts2(v_graph(:,1), 3) - Pts2(v_graph(:,2), 3)).^2);
        G_vertex = graph(v_graph(:,1), v_graph(:,2), e_weight);
        
        index_ccPoints_loop = [curve_vertex_2(:); curve_vertex_2(1)]; % contains the index of the vertices on the cut curve e.g. [1,2,3,4,5,1]
        curve_vertex_output = [];
        for i = 1:(length(index_ccPoints_loop) - 1) 
            s_path = shortestpath(G_vertex, index_ccPoints_loop(i), index_ccPoints_loop(i+1));
            s_path = s_path(:);
            curve_vertex_output = [curve_vertex_output; s_path(1:end-1)];
        end


        edge_graph = [Tri2(:,1), Tri2(:,2), [1:size(Tri2,1)]';Tri2(:,1), Tri2(:,3),[1:size(Tri2,1)]';Tri2(:,2), Tri2(:,3),[1:size(Tri2,1)]'];
        edge_graph(:,1:2) = sort(edge_graph(:,1:2), 2);
        edge_graph = sortrows(edge_graph, [1,2]);
        cut_curve_edges = [curve_vertex_output(1:end),[curve_vertex_output(2:end);curve_vertex_output(1)]];
        cut_curve_edges = sort(cut_curve_edges, 2);
        edge_graph_2 = edge_graph;
        edge_graph_2(ismember(edge_graph_2(:,1:2), cut_curve_edges, 'rows'),:) = [];
        face_G = graph(edge_graph_2(1:2:end,3), edge_graph_2(2:2:end,3));
        [bins, binsizes] = conncomp(face_G);
        face_label_idx = label2idx(bins);
        face_label = bins;
        if(length(binsizes) == 1)
            warning('incomplete separation %s',name_file);
            cut_vertex = [];
            face_fin_label = [];
        else
            [~, sortedID] = sort(cellfun(@length, face_label_idx), 'descend');
            face_label_idx = face_label_idx(sortedID(1:2));
            part1 = face_label_idx{1};
            part2 = face_label_idx{2};
            face_center_label_1 = Tri_center(part1,:);
            face_center_label_2 = Tri_center(part2,:);
            ShaftCoorx = ShaftCoorx*2;
            ShaftCoory = ShaftCoory*2;
            ShaftCoorz = ShaftCoorz*5;
            dist_all_1 = sqrt((ShaftCoorx - face_center_label_1(:,1)').^2 + ...
                (ShaftCoory - face_center_label_1(:,2)').^2 + ...
                (ShaftCoorz - face_center_label_1(:,3)').^2);
            dist_all_2 = sqrt((ShaftCoorx - face_center_label_2(:,1)').^2 + ...
                (ShaftCoory - face_center_label_2(:,2)').^2 + ...
                (ShaftCoorz - face_center_label_2(:,3)').^2);
                face_label = face_label.*0;
            if(mean(min(dist_all_1, [], 1)) > mean(min(dist_all_2, [], 1)))
                face_label(part1) = 2;
                face_label(part2) = 1;
            else
                face_label(part2) = 2;
                face_label(part1) = 1;
            end
            head_neck_label = face_label;
            save(fullfile(output_cut_folder,[name_file,'.face_label.mat']), 'head_neck_label');
            writematrix(curve_vertex_output, fullfile(output_cut_folder, [name_file,'.cut.txt']))
        end




    end

end