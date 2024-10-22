function [error, error_dist,part1_test,part2_test] = checkSegDiff_v2(tri, pts, part1face_example, part1_face_gt, part2_face_gt, edge_graph, target_cycle, test_cycle)
    % error is defined as the ratio between the area between two cycle and
    % the area around each cycle
    error = nan;
    error_dist = nan;
    vectorMF3_1 = pts(tri(:,3),:) - pts(tri(:,1),:);
    vectorMF3_2 = pts(tri(:,3),:) - pts(tri(:,2),:);
    part1_test = 1;
    part2_test = 1;
    
    ss_tmp = cross(vectorMF3_1, vectorMF3_2,2);
    area_each_face = 1/2*(sqrt(ss_tmp(:,1).^2 + ss_tmp(:,2).^2 + ss_tmp(:,3).^2));
    cut_curve_edges = edge_graph(ismember(edge_graph(:,1),test_cycle) & ismember(edge_graph(:,2),test_cycle), 1:2);
    cut_curve_edges = sort(cut_curve_edges, 2);
    edge_graph_2 = edge_graph;
    edge_graph_2(ismember(edge_graph_2(:,1:2), cut_curve_edges, 'rows'),:) = [];
    face_G = graph(edge_graph_2(1:2:end,3), edge_graph_2(2:2:end,3));
    [bins, binsizes] = conncomp(face_G);
    [~, sortedID] = sort(binsizes, 'descend');
    if(length(binsizes) == 1)
        warning('incomplete separation');
    else
        part1_can = find(bins == sortedID(1));
        if(any(part1_can == part1face_example))
            part1_test = find(bins == sortedID(1));
            part2_test = find(bins == sortedID(2));
        else
            part1_test = find(bins == sortedID(2));
            part2_test = find(bins == sortedID(1));
        end
        % the middle area is calculated as whole area subtracted the
        % intersected area
        intersectID = [intersect(part1_face_gt(:), part1_test(:));intersect(part2_face_gt(:), part2_test(:))];
        middle_area = sum(area_each_face) - sum(area_each_face(intersectID));
        cycle_point_list = [target_cycle(:); test_cycle(:)];
        cycle_point_list = unique(cycle_point_list);
        cycleNeiSurface = edge_graph(ismember(edge_graph(:,1), cycle_point_list) & ismember(edge_graph(:,2), cycle_point_list), 3);
        cycleNeiSurface = unique(cycleNeiSurface);
%         cycleNeiSurface = union(cycleNeiSurface,intersectID);
        cycleNeiSurface_area = sum(area_each_face(cycleNeiSurface));
        error = middle_area/ cycleNeiSurface_area;

    end
    dist_average_edge = mean(vecnorm(pts(edge_graph(:,1),:) - pts(edge_graph(:,2),:),2,2));
    dist_min_all = zeros(size(test_cycle,1),1);
    for i = 1:length(dist_min_all)
        dist_min_all(i) = min(vecnorm(pts(test_cycle(i),:) - pts(target_cycle,:),2,2));
    end
    error_dist = mean(dist_min_all)/ dist_average_edge;


end