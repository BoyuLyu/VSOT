% main for comparing the error between different methods
% clear
function out = main_segmentation_error_testing(offFolder,annotationFolder_gt_curves,tamada_result_folder,our_method_cut_result_folder,ofer_output_folder,dorkenwalk_result_folder,save_qauntification_folder)
    % offFolder = '../data/performance_testing_level_2_segmentation/gt_300/surface_off_300';
    % annotationFolder = '../data/performance_testing_level_2_segmentation/gt_300/annotation_json_300';
    % annotationFolder_gt_curves = '../data/performance_testing_level_2_segmentation/gt_300/cut_cycle_ID_300';
    listx = dir([offFolder, '/*.off']);
    % tamada_result_folder = '../data/performance_testing_level_2_segmentation/dendrite_segmentation_peer_methods/tamada/results';
    tamada_error_all = zeros(length(listx), 1) + nan;
    tamada_error_all_dist = tamada_error_all;
    % our_method_cut_result_folder = '../data/performance_testing_level_2_segmentation/dendrite_segmentation_peer_methods/our_method/result';
    our_method_error_all = zeros(length(listx), 1) + nan;
    our_method_error_all_dist = our_method_error_all;
    % ofer_output_folder = '../data/performance_testing_level_2_segmentation/ofer/result';
    ofer_method_error_all = zeros(length(listx), 1) + nan;
    ofer_method_error_all_dist = ofer_method_error_all;
    % dorkenwalk_result_folder = '../data/performance_testing_level_2_segmentation/dorkenwald/result';
    dorkenwalk_error_all = zeros(length(listx), 1) + nan;
    dorkenwalk_error_all_dist = dorkenwalk_error_all;

    for j = 1:length(listx)
        % j = 8 for testing
        namex = listx(j).name;
        namex = namex(1:end-4);
        disp(namex)


        [Pts,Tri] = read_off(fullfile(offFolder, [namex , '.off']));
        Tri = Tri';
        Pts = Pts';
        if(exist(fullfile(annotationFolder_gt_curves, [namex,'.cut.txt']), 'file'))
            gt_cut = readtable(fullfile(annotationFolder_gt_curves, [namex,'.cut.txt']));
            gt_cut = table2array(gt_cut);
            edge_graph = [Tri(:,1), Tri(:,2), [1:size(Tri,1)]';Tri(:,1), Tri(:,3),[1:size(Tri,1)]';Tri(:,2), Tri(:,3),[1:size(Tri,1)]'];
            edge_graph(:,1:2) = sort(edge_graph(:,1:2), 2);
            edge_graph = sortrows(edge_graph, [1,2]);
            % split into two parts
            g_edge = graph(edge_graph(1:2:end,1), edge_graph(1:2:end, 2), vecnorm(Pts(edge_graph(1:2:end,1),:) - Pts(edge_graph(1:2:end,2),:), 2,2));
            d = distances(g_edge);      
            [ss, tt] = find(d == max(d(:)));
            p1 = ss(1);
            p2 = tt(1);
            face1 = find(Tri(:,1) == p1|Tri(:,2) == p1|Tri(:,3) == p1, 1);
            face2 = find(Tri(:,1) == p2|Tri(:,2) == p2|Tri(:,3) == p2, 1); 
    
            cut_curve_edges = [gt_cut(1:end),[gt_cut(2:end);gt_cut(1)]];
            cut_curve_edges = sort(cut_curve_edges, 2);
            edge_graph_2 = edge_graph;
            edge_graph_2(ismember(edge_graph_2(:,1:2), cut_curve_edges, 'rows'),:) = [];
            face_G = graph(edge_graph_2(1:2:end,3), edge_graph_2(2:2:end,3));
            [bins, binsizes] = conncomp(face_G);
            [~, sortedID] = sort(binsizes, 'descend');
            if(length(binsizes) == 1)
                warning('incomplete separation %s',namex);
            else
                part1_can = find(bins == sortedID(1));
                if(any(part1_can == face1))
                    part1 = find(bins == sortedID(1));
                    part2 = find(bins == sortedID(2));
                else
                    part1 = find(bins == sortedID(2));
                    part2 = find(bins == sortedID(1));
                end
            end
            try 
            if(exist(fullfile(tamada_result_folder, [namex,'.cut.txt']), 'file'))
                tamada_cut = readtable(fullfile(tamada_result_folder, [namex,'.cut.txt']));
                tamada_cut = table2array(tamada_cut);
                [error_tamda_tmp, error_tamda_tmp_dist] = checkSegDiff(Tri, Pts, face1, part1, part2, edge_graph, gt_cut, tamada_cut);
                tamada_error_all(j) = error_tamda_tmp;
                tamada_error_all_dist(j) = error_tamda_tmp_dist;
            end
            
            catch ME
                continue;
            end
            try 
    
            if(exist(fullfile(our_method_cut_result_folder,[namex,'.cut.txt']), 'file'))
    
                our_method_cut = readtable(fullfile(our_method_cut_result_folder,[namex,'.cut.txt']));
                our_method_cut = table2array(our_method_cut);
                [error_our_tmp,error_our_tmp_dist]  = checkSegDiff(Tri, Pts, face1, part1, part2, edge_graph, gt_cut, our_method_cut);
                our_method_error_all(j) = error_our_tmp;
                our_method_error_all_dist(j) = error_our_tmp_dist;
            end
            catch ME
                continue;
            end
            try 
                if(exist(fullfile(ofer_output_folder,[namex,'.cut.txt']), 'file'))
        
                    ofer_method_cut = readtable(fullfile(ofer_output_folder,[namex,'.cut.txt']));
                    ofer_method_cut = table2array(ofer_method_cut);
                    [error_ofer_tmp,error_ofer_tmp_dist]= checkSegDiff(Tri, Pts, face1, part1, part2, edge_graph, gt_cut, ofer_method_cut);
                    ofer_method_error_all(j) = error_ofer_tmp;
                    ofer_method_error_all_dist(j) = error_ofer_tmp_dist;
                end
            catch ME
                continue;
            end
            try 
                if(exist(fullfile(dorkenwalk_result_folder,[namex,'.cut.txt']), 'file'))
        
                    dorkenwald_method_cut = readtable(fullfile(dorkenwalk_result_folder,[namex,'.cut.txt']));
                    dorkenwald_method_cut = table2array(dorkenwald_method_cut);
                    [error_dorkenwald_tmp,error_dorkenwald_tmp_dist] = checkSegDiff(Tri, Pts, face1, part1, part2, edge_graph, gt_cut, dorkenwald_method_cut);
                    dorkenwalk_error_all(j) = error_dorkenwald_tmp;
                    dorkenwalk_error_all_dist(j) = error_dorkenwald_tmp_dist;
                end
            catch ME
                continue;
            end
        end
    end
        out.tamada_error_all = tamada_error_all;
        out.tamada_error_all_dist = tamada_error_all_dist;
        out.our_method_error_all = our_method_error_all;
        out.our_method_error_all_dist = our_method_error_all_dist;
        out.ofer_method_error_all = ofer_method_error_all;
        out.ofer_method_error_all_dist = ofer_method_error_all_dist;
        out.dorkenwalk_error_all = dorkenwalk_error_all;
        out.dorkenwalk_error_all_dist = dorkenwalk_error_all_dist;
            % tamada_error_all(isnan(tamada_error_all)) = [];
            % tamada_error_all_dist(isnan(tamada_error_all_dist)) = [];
            % our_method_error_all(isnan(our_method_error_all)) = [];
            % our_method_error_all_dist(isnan(our_method_error_all_dist)) = [];
            % ofer_method_error_all(isnan(ofer_method_error_all)) = [];
            % ofer_method_error_all_dist(isnan(ofer_method_error_all_dist)) = [];
            % 
            % dorkenwalk_error_all(isnan(dorkenwalk_error_all)) = [];
            % dorkenwalk_error_all_dist(isnan(dorkenwalk_error_all_dist)) = [];
    
        % fprintf('area error: VSOT %f, O.N method %f, D.S method %f, T.H method %f \n', ...
        %     median(our_method_error_all), median(ofer_method_error_all), median(dorkenwalk_error_all), median(tamada_error_all))
        % 
        % fprintf('dist error: VSOT %f, O.N method %f, D.S method %f, T.H method %f \n', ...
        %     median(our_method_error_all_dist), median(ofer_method_error_all_dist), median(dorkenwalk_error_all_dist), median(tamada_error_all_dist))
        % 
        % % Save the print information to a txt file
        % fileID = fopen(fullfile(save_qauntification_folder, 'segmentation_error_summary.txt'), 'w');
        % fprintf(fileID, 'area error: VSOT %f, O.N method %f, D.S method %f, T.H method %f \n', ...
        %     median(our_method_error_all), median(ofer_method_error_all), median(dorkenwalk_error_all), median(tamada_error_all));
        % fprintf(fileID, 'dist error: VSOT %f, O.N method %f, D.S method %f, T.H method %f \n', ...
        %     median(our_method_error_all_dist), median(ofer_method_error_all_dist), median(dorkenwalk_error_all_dist), median(tamada_error_all_dist));
        % fclose(fileID);
end