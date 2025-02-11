%% test the methods for level-1 segmentation 
function gen_final_results_all_methods(data_path)
    list_of_folders = ["D5_Apical_Spines","D5_Branch_1","D5_Branch_2","D5_Branch_3","D5_Branch_4","D5_Branch_5"];
    for i = 1:length(list_of_folders)
        current_branch = convertStringsToChars(list_of_folders(i));
        % data_path = '/data';
        rootFolder = [data_path, '/performance_testing_level_1_segmentation/dataset/',current_branch];
        results_folder = [data_path,'/performance_testing_level_1_segmentation/results/'];
        output_face_classification_folder = [data_path,'/performance_testing_level_1_segmentation/results/summary_of_results/', current_branch];
        se0 = strel('sphere', 5);
        % rootFolder = '/work/boyu/EM_astrocyte/test_segmentation_samples/dendrite_spine_segmentation/spinetools_and_our_and_method_3_biotech_result/D5_Apical_Spines';
        [Pts,Tri] = read_off(fullfile(rootFolder, ['dendrite.off']));
        bin_img = tiffreadVolume(fullfile(rootFolder, [current_branch, '_dendrite_volume.tif.tif'])) > 0;

        Tri = Tri';
        Pts = Pts';
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
        %write gt_surface color


        % %% generate our own result of surface 
        us_shaft = tiffreadVolume(fullfile(results_folder, 'VSOT_result' ,current_branch, 'shaft_segmentation.tif')) > 0;
        se1 = strel('sphere',2);
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
        %% organize the result from spinetool
        % read in all the points that are detected as spines

        % sptools_result_path = fullfile(results_folder, 'spineTool_results', current_branch,'spinetools_result');
        % all_file_name = dir(fullfile(sptools_result_path, '*.off'));
        % points_spine = [];
        % for j = 1:length(all_file_name)
        %     [Pts2,Tri2] = read_off(fullfile(sptools_result_path, all_file_name(j).name));
        %     Tri2 = Tri2';
        %     Pts2 = Pts2';
        %     points_spine = [points_spine; Pts2];
        % end
        % vert_idx = knnsearch(Pts, points_spine, 'K', 1);
        % vert_part1 = unique(vert_idx);
        % spine_faces_id = find(ismember(Tri(:,1), vert_part1(:)) & ismember(Tri(:,2), vert_part1(:)) & ismember(Tri(:,3), vert_part1(:)));
        % output_spinetools = ones(size(face_centers,1),1);
        % output_spinetools(spine_faces_id(:)) = 2;
        % disp('spinetool finished')
        % % %% neurd 
        % neurd_root_folder = fullfile(results_folder, 'neurd_results');
        % neurd_result_path = fullfile(neurd_root_folder, current_branch);
        % all_file_name = dir(fullfile(neurd_result_path, '*.off'));
        % points_spine = [];
        % for j = 1:length(all_file_name)
        %     [Pts2,Tri2] = read_off(fullfile(neurd_result_path, all_file_name(j).name));
        %     Tri2 = Tri2';
        %     Pts2 = Pts2';
        %     points_spine = [points_spine; Pts2];
        % end
        % vert_idx = knnsearch(Pts, points_spine, 'K', 1);
        % vert_part1 = unique(vert_idx);
        % spine_faces_id = find(ismember(Tri(:,1), vert_part1(:)) & ismember(Tri(:,2), vert_part1(:)) & ismember(Tri(:,3), vert_part1(:)));
        % output_neurd = ones(size(face_centers,1),1);
        % output_neurd(spine_faces_id(:)) = 2;
        % disp('output_neurd finished')
        % %% rootFolder = '/work/boyu/EM_astrocyte/test_segmentation_samples/dendrite_spine_segmentation/bio_tech_paper_result/D5_Branch_4';
        % %% morphological method
        % se0 = strel('sphere', 4);
        % bin_img2 = imopen(bin_img, se0);
        % bin_img2_roi = bwlabeln(bin_img2);
        % bin_img2_roi_idx = label2idx(bin_img2_roi);
        % bin_img2 = bin_img2.*0;
        % len_cell = cellfun(@length, bin_img2_roi_idx);
        % bin_img2(bin_img2_roi_idx{(len_cell == max(len_cell))}) = 1;
        % se1 = strel('sphere', 2);
        % bin_img3 = imdilate(bin_img2, se1);
        % % figure; volshow(bin_img3)
        % spine_region = (bin_img - bin_img3) > 0;
        % spine_region_roi = bwlabeln(spine_region);
        % spine_region_roi_idx = label2idx(spine_region_roi);
        % spine_region_roi_idx(cellfun(@length, spine_region_roi_idx) < 3) = [];
        % spine_region_roi_idx = spine_region_roi_idx(:);

        % new_combined = double(bin_img);
        % new_combined(cell2mat(spine_region_roi_idx)) = 2;

        % face_centers = (Pts(Tri(:,1),:) + Pts(Tri(:,2),:)  + Pts(Tri(:,3),:) )/3;
        % %find the nearest faces to t
        % dict_points = find(bin_img(:) == 1);
        % [lenx, leny, lenz] = size(bin_img);
        % [dict_pointsx, dict_pointsy, dict_pointsz] = ind2sub([lenx, leny, lenz], dict_points);

        % dict_pointsxyz = [dict_pointsx(:), dict_pointsy(:), dict_pointsz(:)];
        % idx = knnsearch(dict_pointsxyz, face_centers, 'K', 1);
        % label_face = new_combined(dict_points(idx));
        % face_part1 = find(label_face == 1);
        % face_part2 = find(label_face == 2);
        % output_morph = ones(size(face_centers,1),1);
        % output_morph(face_part1(:)) = 1;
        % output_morph(face_part2(:)) = 2;
        % disp('output_morph finished')
        % bio_tech_root_folder = fullfile(results_folder, 'bio_tech_paper_result', current_branch);
        % % bin_img = tiffreadVolume(fullfile(rootFolder, [current_branch, '_dendrite_volume.tif.tif'])) > 0;

        % bin_img2 = imopen(bin_img, se0);
        % bin_img2_roi = bwlabeln(bin_img2);
        % bin_img2_roi_idx = label2idx(bin_img2_roi);
        % bin_img2 = bin_img2.*0;
        % len_cell = cellfun(@length, bin_img2_roi_idx);
        % bin_img2(bin_img2_roi_idx{(len_cell == max(len_cell))}) = 1;
        % se1 = strel('sphere', 1);
        % bin_img3 = imdilate(bin_img2, se1);
        % % figure; volshow(bin_img3)
        % spine_region = (bin_img - bin_img3) > 0;
        % spine_region_roi = bwlabeln(spine_region);
        % spine_region_roi_idx = label2idx(spine_region_roi);
        % % remove the small components that are caused by the process of erosion
        % spine_region_roi_idx(cellfun(@length, spine_region_roi_idx) < 3) = []; 
        % spine_region_roi_idx = spine_region_roi_idx(:);
        % new_combined = zeros(size(spine_region));
        % new_combined = double(bin_img);
        % new_combined(cell2mat(spine_region_roi_idx)) = 2;
        % % figure; volshow(new_combined)
        % [Pts,Tri] = read_off(fullfile(bio_tech_root_folder, ['dendrite.off']));
        % Tri = Tri';
        % Pts = Pts';
        % colormapx = [255, 64, 64, 255;64, 255, 255,255;64, 255, 64, 255];
        % node_colorMap = zeros(size(Pts,1),4);
        % face_centers = (Pts(Tri(:,1),:) + Pts(Tri(:,2),:)  + Pts(Tri(:,3),:) )/3;
        % %find the nearest faces to t
        % dict_points = find(bin_img(:) == 1);
        % [lenx, leny, lenz] = size(bin_img);
        % [dict_pointsx, dict_pointsy, dict_pointsz] = ind2sub([lenx, leny, lenz], dict_points);

        % dict_pointsxyz = [dict_pointsx(:), dict_pointsy(:), dict_pointsz(:)];
        % idx = knnsearch(dict_pointsxyz, face_centers, 'K', 1);
        % label_face = new_combined(dict_points(idx));
        % face_part1 = find(label_face == 1);
        % face_part2 = find(label_face == 2);
        % output_morph = ones(size(face_centers,1),1);
        % output_morph(face_part1(:)) = 1;
        % output_morph(face_part2(:)) = 2;
        % disp('output_morph finished')

        % writematrix(output_gt, fullfile(output_face_classification_folder, 'gt_face_classification_result.txt'));
        writematrix(output_us, fullfile(output_face_classification_folder, 'our_face_classification_result.txt'));
        % writematrix(output_neurd, fullfile(output_face_classification_folder, 'neurd_face_classification_result.txt'));
        % writematrix(output_spinetools, fullfile(output_face_classification_folder, 'spinetool_face_classification_result.txt'));
        % writematrix(output_morph, fullfile(output_face_classification_folder, 'morph_face_classification_result.txt'));
    end
end