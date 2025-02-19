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


%% check spinetool result, which is best
results_folder = '/work/boyu/EM_astrocyte/materials_for_paper_VSOT/test_segmentation_samples/dendrite_spine_segmentation/spinetool_result/D5_Branch_1_dendrite_volume.tif';

p1 = [-3:0.3:-0.3];
p2 = [-3:1:2];
f1_score_all = zeros(length(p1), length(p2));
for m = 1:length(p1)
    for n = 1:length(p2)
    points_spine = [];
    curfolder = fullfile(results_folder,  num2str((m-1)), num2str((n-1)));
    all_file_name = dir(fullfile(curfolder, '*.off'));
        for j = 1:length(all_file_name)
            [Pts2,Tri2] = read_off(fullfile(curfolder, all_file_name(j).name));
            Tri2 = Tri2';
            Pts2 = Pts2';
            points_spine = [points_spine; Pts2];
        end
    vert_idx = knnsearch(Pts, points_spine, 'K', 1);
    vert_part1 = unique(vert_idx);
    spine_faces_id = find(ismember(Tri(:,1), vert_part1(:)) & ismember(Tri(:,2), vert_part1(:)) & ismember(Tri(:,3), vert_part1(:)));
    output_spinetools = ones(size(face_centers,1),1);
    output_spinetools(spine_faces_id(:)) = 2;
    disp('spinetool finished')
    [F1score_spinetool, p_spinetool, r_spinetool, J_spinetool] = level_1_performance_testing_helper(output_gt, output_spinetools, Tri, triangle_area);
    f1_score_all(m, n) = F1score_spinetool;
    end
end

writematrix(f1_score_all, fullfile(neurd_result_path, 'f1_score_map.csv'));



[m,n] = find(f1_score_all == max(f1_score_all(:)));

figure; 
for i = 1:size(f1_score_all,2)
    smoothness = p1;
    hold on; plot(f1_score_all(:,i), 'LineWidth',2);
end
legend('correction -3', 'correction -2','correction -1','correction 0','correction 1','correction 2');
xticks(1:length(p1)); xticklabels(p1);
xlabel('log(sensitivity)'); ylabel('F1-score');title('SpineTool Parameter Tuning')