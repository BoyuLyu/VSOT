%% biotechnology paper that uses erosion and then dilation to obtain the dendrite spines
% implementation in MATLAB\

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
%% check different kernel size for erosion and dilation
kernel1 = [1,2,3,4];
kernel2 = [1,2,3,4];
f1_score_all = zeros(length(kernel1), length(kernel2));
for m = 1:length(kernel1)
    for n = 1:length(kernel2)
        se0 = strel('sphere', kernel1(m));
        bin_img2 = imopen(bin_img, se0);
        bin_img2_roi = bwlabeln(bin_img2);
        bin_img2_roi_idx = label2idx(bin_img2_roi);
        bin_img2 = bin_img2.*0;
        len_cell = cellfun(@length, bin_img2_roi_idx);
        bin_img2(bin_img2_roi_idx{(len_cell == max(len_cell))}) = 1;
        se1 = strel('sphere', kernel2(n));
        bin_img3 = imdilate(bin_img2, se1);
        % figure; volshow(bin_img3)
        spine_region = (bin_img - bin_img3) > 0;
        spine_region_roi = bwlabeln(spine_region);
        spine_region_roi_idx = label2idx(spine_region_roi);
        spine_region_roi_idx(cellfun(@length, spine_region_roi_idx) < 3) = [];
        spine_region_roi_idx = spine_region_roi_idx(:);

        new_combined = double(bin_img);
        new_combined(cell2mat(spine_region_roi_idx)) = 2;
        node_colorMap = zeros(size(Pts,1),4);
        face_centers = (Pts(Tri(:,1),:) + Pts(Tri(:,2),:)  + Pts(Tri(:,3),:) )/3;
        %find the nearest faces to t
        dict_points = find(bin_img(:) == 1);
        [lenx, leny, lenz] = size(bin_img);
        [dict_pointsx, dict_pointsy, dict_pointsz] = ind2sub([lenx, leny, lenz], dict_points);

        dict_pointsxyz = [dict_pointsx(:), dict_pointsy(:), dict_pointsz(:)];
        idx = knnsearch(dict_pointsxyz, face_centers, 'K', 1);
        label_face = new_combined(dict_points(idx));
        face_part1 = find(label_face == 1);
        face_part2 = find(label_face == 2);
        output_morph = ones(size(face_centers,1),1);
        output_morph(face_part1(:)) = 1;
        output_morph(face_part2(:)) = 2;
        disp('output_morph finished')
        [F1score_morph, p_morph, r_morph, J_morph] = level_1_performance_testing_helper(output_gt, output_morph, Tri, triangle_area);
        f1_score_all(m, n) = F1score_morph;
    end
end
[m,n] = find(f1_score_all == max(f1_score_all(:)));

figure; 
for i = 2:size(f1_score_all,1)
    smoothness = p1;
    hold on; plot(f1_score_all(i,:), 'LineWidth',2);
end
legend('Erosion kernel size 2', 'Erosion kernel size 3','Erosion kernel size 4');
xticks(1:length(kernel2)); xticklabels(kernel2);
xlabel('Dilation kenel size'); ylabel('F1-score');title('Morph Parameter Tuning')