
current_branch = 'D5_Branch_1';
data_path = 'D:\Guoqiang_lab\VSOT_samples\VSOT_data_paired_with_codes\data';
rootFolder = [data_path, '/performance_testing_level_1_segmentation/dataset/',current_branch];
[Pts,Tri] = read_off(fullfile(rootFolder, ['dendrite.off']));
bin_img = tiffreadVolume(fullfile(rootFolder, [current_branch, '_dendrite_volume.tif.tif'])) > 0;

Tri = Tri';
Pts = Pts';
vectorMF3_1 = Pts(Tri(:,3),:) - Pts(Tri(:,1),:);
vectorMF3_2 = Pts(Tri(:,3),:) - Pts(Tri(:,2),:);
ss_tmp = cross(vectorMF3_1, vectorMF3_2,2);
ss = 1/2*(sqrt(ss_tmp(:,1).^2 + ss_tmp(:,2).^2 + ss_tmp(:,3).^2));

triangle_area = ss;
% colormapx = [255, 64, 64, 255;64, 255, 255,255;64, 255, 64, 255];
% node_colorMap = zeros(size(Pts,1),4);
face_centers = (Pts(Tri(:,1),:) + Pts(Tri(:,2),:)  + Pts(Tri(:,3),:) )/3;
dict_points = find(bin_img(:) == 1);
[lenx, leny, lenz] = size(bin_img);
[dict_pointsx, dict_pointsy, dict_pointsz] = ind2sub([lenx, leny, lenz], dict_points);

dict_pointsxyz = [dict_pointsx(:), dict_pointsy(:), dict_pointsz(:)];
idx = knnsearch(dict_pointsxyz, face_centers, 'K', 1);
%% ground truth
gt_shaft = tiffreadVolume(fullfile(rootFolder, [current_branch,'_shaft_volume.tif.tif'])) > 0;
gt_spine = bin_img - gt_shaft > 0;
gt_combined = gt_shaft*1 + gt_spine*2;
label_face_gt = gt_combined(dict_points(idx));
face_part1 = find(label_face_gt == 1);
face_part2 = find(label_face_gt == 2);
output_gt = zeros(size(face_centers,1),1);
output_gt(face_part1(:)) = 1;
output_gt(face_part2(:)) = 2; % spines are labeled as 2
% disp('gt finished')
%% check the result from neurd, 1st cgal score tuning
p2 = [5:10];
p1 = [0.01:0.01:0.09];
neurd_result_path = '/work/boyu/EM_astrocyte/materials_for_paper_VSOT/VSOT_data_paired_with_codes/data/performance_testing_level_1_segmentation/results/neurd_results/output_filtering_tune_parameter'
f1_score_map = zeros(length(p1), length(p2));
for i = 1:length(p1)
    for j = 1:length(p2)
        saved_result_path = fullfile(neurd_result_path, num2str(i - 1), num2str(j - 1), 'classification_result');
        if(~exist(saved_result_path, 'dir'))
            mkdir(saved_result_path);
        end
        cur_folder = fullfile(neurd_result_path, num2str(i - 1), num2str(j - 1));
        all_file_name = dir(fullfile(cur_folder, '*.off'));
        points_spine = [];
        for m = 1:length(all_file_name)
            [Pts2,Tri2] = read_off(fullfile(cur_folder, all_file_name(m).name));
            Tri2 = Tri2';
            Pts2 = Pts2';
            points_spine = [points_spine; Pts2];
        end
        vert_idx = knnsearch(Pts, points_spine, 'K', 1);
        vert_part1 = unique(vert_idx);
        spine_faces_id = find(ismember(Tri(:,1), vert_part1(:)) & ismember(Tri(:,2), vert_part1(:)) & ismember(Tri(:,3), vert_part1(:)));
        % face_idx = knnsearch(face_centers, points_spine, 'K',1);
        % spine_faces_id = unique(face_idx);
        output_neurd = ones(size(face_centers,1),1);
        output_neurd(spine_faces_id(:)) = 2;
        disp('output_neurd finished')
        writematrix(output_neurd, fullfile(saved_result_path, 'neurd_face_classification_result.txt'));
        [F1score_neurd, p_neurd, r_neurd, J_neurd] = level_1_performance_testing_helper(output_gt, output_neurd, Tri, triangle_area);
        f1_score_map(i,j) = F1score_neurd;
    end
end




writematrix(f1_score_map, fullfile(neurd_result_path, 'f1_score_map.csv'));





figure; 
for i = 1:size(f1_score_map,2)
    smoothness = p1;
    hold on; plot(f1_score_map(:,i), 'LineWidth',2);
end
legend('cluster number 5', 'cluster number 6', 'cluster number 7','cluster number 8','cluster number 9','cluster number 10')
% record all the folder IDs as well as the corresponding f1 score
%
xticks(1:length(p1)); xticklabels(p1)
%% check the result, part 2 tuning

p2 = [10:20:100]; %spine_sk_length_threshold_bare_min
p1 = [50:50:450];%filter_by_volume_threshold_bare_min
% neurd_root_folder = '/work/boyu/EM_astrocyte/materials_for_paper_VSOT/test_segmentation_samples/dendrite_spine_segmentation/neurd_results/Spine_Detection_On_Mesh_Branch/notebooks/notebooks'
neurd_result_path = 'D:\Guoqiang_lab\VSOT_samples\VSOT_data_paired_with_codes\data/performance_testing_level_1_segmentation/results/neurd_results/output_filtering_tune_parameter2'
f1_score_map = zeros(length(p1), length(p2));
for i = 1:length(p1)
    for j = 1:length(p2)
        saved_result_path = fullfile(neurd_result_path, num2str(i - 1), num2str(j - 1), 'classification_result');
        if(~exist(saved_result_path, 'dir'))
            mkdir(saved_result_path);
        end
        cur_folder = fullfile(neurd_result_path, num2str(i - 1), num2str(j - 1));
        all_file_name = dir(fullfile(cur_folder, '*.off'));
        points_spine = [];
        for m = 1:length(all_file_name)
            [Pts2,Tri2] = read_off(fullfile(cur_folder, all_file_name(m).name));
            Tri2 = Tri2';
            Pts2 = Pts2';
            points_spine = [points_spine; Pts2];
        end
        vert_idx = knnsearch(Pts, points_spine, 'K', 1);
        vert_part1 = unique(vert_idx);
        spine_faces_id = find(ismember(Tri(:,1), vert_part1(:)) & ismember(Tri(:,2), vert_part1(:)) & ismember(Tri(:,3), vert_part1(:)));
        % face_idx = knnsearch(face_centers, points_spine, 'K',1);
        % spine_faces_id = unique(face_idx);
        output_neurd = ones(size(face_centers,1),1);
        output_neurd(spine_faces_id(:)) = 2;
        disp('output_neurd finished')
        writematrix(output_neurd, fullfile(saved_result_path, 'neurd_face_classification_result.txt'));
        [F1score_neurd, p_neurd, r_neurd, J_neurd] = level_1_performance_testing_helper(output_gt, output_neurd, Tri, triangle_area);
        f1_score_map(i,j) = F1score_neurd;
    end
end


[m, n] = find(f1_score_map == max(f1_score_map(:)))

writematrix(f1_score_map, fullfile(neurd_result_path, 'f1_score_map.csv'));





figure; 
for i = 1:size(f1_score_map,2)
    smoothness = p1;
    hold on; plot(f1_score_map(:,i), 'LineWidth',2);
end
legend('spine length 10', 'spine length 30', 'spine length 50','spine length 70','spine length 90')
% record all the folder IDs as well as the corresponding f1 score
%
xticks(1:length(p1)); xticklabels(p1); ylabel('F1 score')
xlabel('spine area threshold')