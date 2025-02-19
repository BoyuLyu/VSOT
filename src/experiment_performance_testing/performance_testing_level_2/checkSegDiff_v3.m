function [F1score, precision, recall, IOU] = checkSegDiff_v3(tri, pts, face_label_gt, face_label_test, edge_graph, target_cycle, test_cycle)
    % error is defined as the ratio between the area between two cycle and
    % the area around each cycle
    F1score = nan;
    precision = nan;
    recall = nan;
    IOU = nan;
    vectorMF3_1 = pts(tri(:,3),:) - pts(tri(:,1),:);
    vectorMF3_2 = pts(tri(:,3),:) - pts(tri(:,2),:);
    ss_tmp = cross(vectorMF3_1, vectorMF3_2,2);
    area_each_face = 1/2*(sqrt(ss_tmp(:,1).^2 + ss_tmp(:,2).^2 + ss_tmp(:,3).^2));
    intersected_ID = find(face_label_gt == face_label_test);
    intersected_ID = intersected_ID(:);

    TP_id = find((face_label_gt == 2) & (face_label_test == 2));
    FP_id = find((face_label_gt == 1) & (face_label_test == 2));
    FN_id = find((face_label_gt == 2) & (face_label_test == 1));
    F1score = (2*sum(area_each_face(TP_id)))/(sum(2*area_each_face(TP_id)) + sum(area_each_face(FP_id)) + sum(area_each_face(FN_id)));
    PP_id = find(face_label_test == 2);
    P_id = find(face_label_gt == 2);
    precision = sum(area_each_face(TP_id))/ sum(area_each_face(PP_id));
    recall = sum(area_each_face(TP_id))/ sum(area_each_face(P_id));
    IOU = sum(area_each_face(intersected_ID))/ sum(area_each_face);

%     middle_area = sum(area_each_face) - sum(area_each_face(intersected_ID));
%     cycle_point_list = [target_cycle(:); test_cycle(:)];
%     cycle_point_list = unique(cycle_point_list);
%     cycleNeiSurface = edge_graph(ismember(edge_graph(:,1), cycle_point_list) & ismember(edge_graph(:,2), cycle_point_list), 3);
%     cycleNeiSurface = unique(cycleNeiSurface);
% %         cycleNeiSurface = union(cycleNeiSurface,intersectID);
%     cycleNeiSurface_area = sum(area_each_face(cycleNeiSurface));
%     error = middle_area/ cycleNeiSurface_area;

%     dist_average_edge = mean(vecnorm(pts(edge_graph(:,1),:) - pts(edge_graph(:,2),:),2,2));
%     dist_min_all = zeros(size(test_cycle,1),1);
%     for i = 1:length(dist_min_all)
%         dist_min_all(i) = min(vecnorm(pts(test_cycle(i),:) - pts(target_cycle,:),2,2));
%     end
%     error_dist = mean(dist_min_all)/ dist_average_edge;



end