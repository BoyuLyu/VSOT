current_branch = 'D5_Branch_1';
data_path = '/work/boyu/EM_astrocyte/materials_for_paper_VSOT/VSOT_data_paired_with_codes/data';
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
%% our result
results_folder = '/work/boyu/EM_astrocyte/materials_for_paper_VSOT/test_segmentation_samples/dendrite_spine_segmentation/our_method_result/D5_Branch_1';
lambda_all = [0.25:0.25:2];
f1_score_all = zeros(length(lambda_all),1);
for i = 1:length(lambda_all)
    curFolder = fullfile(results_folder, ['lambda_index_', num2str(i)]);
    us_shaft = tiffreadVolume(fullfile(curFolder, 'shaft_segmentation.tif')) > 0;
    se1 = strel('sphere',1);
    us_shaft = imdilate(us_shaft, se1);
    us_spine = bin_img - us_shaft > 0;
    us_combined = us_shaft*1 + us_spine*2;
    label_face_us = us_combined(dict_points(idx));
    face_part1 = find(label_face_us == 1);
    face_part2 = find(label_face_us == 2);
    output_us = zeros(size(face_centers,1),1);
    output_us(face_part1(:)) = 1;
    output_us(face_part2(:)) = 2;
    disp('VSOT finished')
    [F1score_vsot, p_vsot, r_vsot, J_vsot] = level_1_performance_testing_helper(output_gt, output_us, Tri, triangle_area);
    f1_score_all(i) = F1score_vsot;
end
