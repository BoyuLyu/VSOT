function save_performance_testing_results(data_path, quantification_save_path)
    %% test the methods for level-1 segmentation 
    us_score = zeros(6,4);
    neurd_score = zeros(6,4);
    spinetool_score = zeros(6,4);
    morph_score = zeros(6,4);
    branchlist = [{'D5_Branch_1'}, {'D5_Branch_2'}, {'D5_Branch_3'}, {'D5_Branch_4'},{'D5_Branch_5'}, {'D5_Apical_Spines'}];
    for j = 1:length(branchlist)
        % disp(j)
        current_branch = branchlist{j};
        rootFolder = [data_path,'/performance_testing_level_1_segmentation/results/summary_of_results/', current_branch];
        rootFolder_0 = [data_path, '/performance_testing_level_1_segmentation/dataset/', current_branch];
        % rootFolder = ['/work/boyu/EM_astrocyte/test_segmentation_samples/dendrite_spine_segmentation/spinetools_and_our_and_method_3_biotech_result/',current_branch];
        % rootFolder = '/work/boyu/EM_astrocyte/test_segmentation_samples/dendrite_spine_segmentation/spinetools_and_our_and_method_3_biotech_result/D5_Apical_Spines';
        [Pts,Tri] = read_off(fullfile(rootFolder_0, ['dendrite.off']));
        % bin_img = tiffreadVolume(fullfile(rootFolder, [current_branch, '_dendrite_volume.tif.tif'])) > 0;
        
        Tri = Tri';
        Pts = Pts';
        vectorMF3_1 = Pts(Tri(:,3),:) - Pts(Tri(:,1),:);
        vectorMF3_2 = Pts(Tri(:,3),:) - Pts(Tri(:,2),:);
        ss_tmp = cross(vectorMF3_1, vectorMF3_2,2);
        ss = 1/2*(sqrt(ss_tmp(:,1).^2 + ss_tmp(:,2).^2 + ss_tmp(:,3).^2));
        
        triangle_area = ss;
        
        output_gt = table2array(readtable(fullfile(rootFolder, 'gt_face_classification_result.txt')));
        output_us = table2array(readtable(fullfile(rootFolder, 'our_face_classification_result.txt')));
        output_neurd = table2array(readtable(fullfile(rootFolder, 'neurd_face_classification_result.txt')));
        output_spinetools = table2array(readtable(fullfile(rootFolder, 'spinetool_face_classification_result.txt')));
        output_morph = table2array(readtable(fullfile(rootFolder, 'morph_face_classification_result.txt')));
        
        %% consider F1 score
        [F1score_us, p_us, r_us, J_us] = level_1_performance_testing_helper(output_gt, output_us, Tri, triangle_area);
        [F1score_neurd, p_neurd, r_neurd, J_neurd] = level_1_performance_testing_helper(output_gt, output_neurd, Tri, triangle_area);
        [F1score_spinetool, p_spinetool, r_spinetool, J_spinetool] = level_1_performance_testing_helper(output_gt, output_spinetools, Tri, triangle_area);
        [F1score_morph, p_morph, r_morph, J_morph] = level_1_performance_testing_helper(output_gt, output_morph, Tri, triangle_area);
        us_score(j,:) = [F1score_us, p_us, r_us, J_us];
        neurd_score(j,:) = [F1score_neurd, p_neurd, r_neurd, J_neurd];
        spinetool_score(j,:) = [F1score_spinetool, p_spinetool, r_spinetool, J_spinetool];
        morph_score(j,:) = [F1score_morph, p_morph, r_morph, J_morph];
    end

    fileID = fopen(fullfile(quantification_save_path, 'level_1_segmentation_accuracy.txt'), 'w');
    fprintf(fileID, 'VSOT: F1 %f, precision %f,recall %f, Jeccard %f \n', ...
        mean(us_score(:,1)), mean(us_score(:,2)), mean(us_score(:,3)),mean(us_score(:,4)));
    fprintf(fileID, 'Neurd: F1 %f, precision %f,recall %f, Jeccard %f \n', ...
        mean(neurd_score(:,1)), mean(neurd_score(:,2)), mean(neurd_score(:,3)),mean(neurd_score(:,4)));
    fprintf(fileID, 'Spinetool: F1 %f, precision %f,recall %f, Jeccard %f \n', ...
        mean(spinetool_score(:,1)), mean(spinetool_score(:,2)), mean(spinetool_score(:,3)),mean(spinetool_score(:,4)));
    fprintf(fileID, 'Morph: F1 %f, precision %f,recall %f, Jeccard %f \n', ...
        mean(morph_score(:,1)), mean(morph_score(:,2)), mean(morph_score(:,3)),mean(morph_score(:,4)));
    fclose(fileID);


end