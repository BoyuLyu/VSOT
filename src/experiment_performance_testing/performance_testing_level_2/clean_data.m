offFolder = '/work/boyu/EM_astrocyte/materials_for_paper_VSOT/VSOT_data_paired_with_codes/data/performance_testing_level_2_segmentation/large_dataset_3_annotator_L2345_w_stubby/gt_final_800/off_file';
annotationFolder = '/work/boyu/EM_astrocyte/materials_for_paper_VSOT/VSOT_data_paired_with_codes/data/performance_testing_level_2_segmentation/large_dataset_3_annotator_L2345_w_stubby/gt_final_800/annotation_annotator_1';
annotationFolder_gt_curves = '/work/boyu/EM_astrocyte/materials_for_paper_VSOT/VSOT_data_paired_with_codes/data/performance_testing_level_2_segmentation/large_dataset_3_annotator_L2345_w_stubby/gt_final_800/cutID_annotator_1';
if(~exist(annotationFolder_gt_curves, 'dir'))
    mkdir(annotationFolder_gt_curves)

end
num_control_points1 = gen_gt_from_controlPoints(offFolder, annotationFolder, annotationFolder_gt_curves);
% figure; trisurf(Tri, Pts(:,1), Pts(:,2), Pts(:,3),'Facecolor','red','FaceAlpha',0.1); hold on;plot3(Pts(curve_vertex, 1), Pts(curve_vertex, 2), Pts(curve_vertex, 3), 'LineWidth', 3)
annotationFolder = '/work/boyu/EM_astrocyte/materials_for_paper_VSOT/VSOT_data_paired_with_codes/data/performance_testing_level_2_segmentation/large_dataset_3_annotator_L2345_w_stubby/gt_final_800/annotation_annotator_2';
annotationFolder_gt_curves = '/work/boyu/EM_astrocyte/materials_for_paper_VSOT/VSOT_data_paired_with_codes/data/performance_testing_level_2_segmentation/large_dataset_3_annotator_L2345_w_stubby/gt_final_800/cutID_annotator_2';
if(~exist(annotationFolder_gt_curves, 'dir'))
    mkdir(annotationFolder_gt_curves)

end
num_control_points2 = gen_gt_from_controlPoints(offFolder, annotationFolder, annotationFolder_gt_curves);

annotationFolder = '/work/boyu/EM_astrocyte/materials_for_paper_VSOT/VSOT_data_paired_with_codes/data/performance_testing_level_2_segmentation/large_dataset_3_annotator_L2345_w_stubby/gt_final_800/annotation_annotator_3';
annotationFolder_gt_curves = '/work/boyu/EM_astrocyte/materials_for_paper_VSOT/VSOT_data_paired_with_codes/data/performance_testing_level_2_segmentation/large_dataset_3_annotator_L2345_w_stubby/gt_final_800/cutID_annotator_3';
if(~exist(annotationFolder_gt_curves, 'dir'))
    mkdir(annotationFolder_gt_curves)

end
num_control_points3 = gen_gt_from_controlPoints(offFolder, annotationFolder, annotationFolder_gt_curves);

gt_folder = '/work/boyu/EM_astrocyte/materials_for_paper_VSOT/VSOT_data_paired_with_codes/data/performance_testing_level_2_segmentation/large_dataset_3_annotator_L2345_w_stubby/gt_final_800';
an1 = 'cutID_annotator_1';
an2 = 'cutID_annotator_2';
an3 = 'cutID_annotator_3';
volume_w_shaft_folder = 'gt_volume_w_shaft';
gen_ground_truth(gt_folder, an1, an2,an3,volume_w_shaft_folder)


rootfolder = '/work/boyu/EM_astrocyte/materials_for_paper_VSOT/VSOT_data_paired_with_codes/data/performance_testing_level_2_segmentation/large_dataset_3_annotator_L2345_w_stubby/gt_final_800/stubby_samples';
offFolder = fullfile(rootfolder, 'off_file/');
tifFolder = fullfile(rootfolder, 'tif_file/');
output_cut_folder = fullfile(rootfolder, 'cut_index/');
gen_gt_from_volume_spine_shaft(offFolder, tifFolder, output_cut_folder)
figure; trisurf(Tri2, Pts2(:,1), Pts2(:,2), Pts2(:,3),'Facecolor','red','FaceAlpha',0.1); hold on;plot3(Pts2(curve_vertex_output, 1), Pts2(curve_vertex_output, 2), Pts2(curve_vertex_output, 3), 'LineWidth', 3)
figure; trisurf(Tri, Pts(:,1), Pts(:,2), Pts(:,3),'Facecolor','red','FaceAlpha',0.1); hold on;plot3(Pts(curve_vertex, 1), Pts(curve_vertex, 2), Pts(curve_vertex, 3), 'LineWidth', 3)
figure; trisurf(Tri2(part1,:), Pts2(:,1), Pts2(:,2), Pts2(:,3),'Facecolor','red'); hold on; 
trisurf(Tri2(part2,:), Pts2(:,1), Pts2(:,2), Pts2(:,3),'Facecolor','blue');
%assign the label to head and neck part of the surface 
gt_folder = '/work/boyu/EM_astrocyte/materials_for_paper_VSOT/VSOT_data_paired_with_codes/data/performance_testing_level_2_segmentation/large_dataset_3_annotator_L2345_w_stubby/gt_final_800/uniformly_selected_samples';
offFolder = fullfile(gt_folder, 'off_file/');
an1 = fullfile(gt_folder,'cutID_annotator_1');
an2 = fullfile(gt_folder, 'cutID_annotator_2');
an3 = fullfile(gt_folder, 'cutID_annotator_3');
tifFolder = fullfile(gt_folder, 'gt_volume_w_shaft/');

gen_face_head_neck_label(offFolder, an1, an2, an3, tifFolder)

out_record = check_gt_consistency(offFolder, an1, an2, an3);
writetable(out_record, fullfile(gt_folder, 'annotator_between_iou2.csv'));
iou_score = table2array(out_record(:,2:end));
selectedID = []; % the 2:4 column contains the selection of annotator
threshold1 = 0.9;
threshold2 = 0.95;
for i = 1:size(out_record,1)
    if(all(iou_score(i,:) > threshold1))
        selectedID = [selectedID; [i, 1,1,1]];
    elseif(any(iou_score(i,:) > threshold2))
        idx = find(iou_score(i,:) == max(iou_score(i,:))); % should output only one
        % 1 -> 1 2
        % 2 -> 1 3
        % 3 -> 2 3
        if(idx == 1)
            selectedID = [selectedID; [i, 1,1,0]];
        elseif(idx == 2)
            selectedID = [selectedID; [i, 1,0,1]];
        else
            selectedID = [selectedID; [i, 0,1,1]];
        end

    end
end
mean_curve_folder = fullfile(gt_folder, 'mean_annotation_curve');
find_mean_curve(selectedID, out_record, offFolder, an1, an2, an3, gt_folder)



out_map = nan(size(Vnieuw2, 1),1);
out_map(idSelected) = ccScore;
% out_map(neckpoint) = 1;
% out_map(id_point) = 2;
% out_map(ccNeck) = 1;
figure;trisurf(mF2,Vnieuw2(:,1),Vnieuw2(:,2),Vnieuw2(:,3),'FaceVertexCData',out_map)
hold on; plot3(Vnieuw2(cutVexcellUnique2{id_point},1),Vnieuw2(cutVexcellUnique2{id_point},2),Vnieuw2(cutVexcellUnique2{id_point},3), 'LineWidth',3)


out_map = nan(size(Vnieuw2, 1),1);
out_map(ccNeck) = 1;
figure;trisurf(mF2,Vnieuw2(:,1),Vnieuw2(:,2),Vnieuw2(:,3),'FaceVertexCData',out_map)
hold on; plot3(pathidxyz(:,1), pathidxyz(:,2), pathidxyz(:,3), 'LineWidth',3)

hold on; plot3(Vnieuw2(cc_vertexCell{id_point},1),Vnieuw2(cc_vertexCell{id_point},2),Vnieuw2(cc_vertexCell{id_point},3), 'LineWidth',3)


 figure; trisurf(Tri(face_label == 2,:), Pts(:,1), Pts(:,2), Pts(:,3), 'Facecolor','red');
 hold on; trisurf(Tri(face_label == 1,:), Pts(:,1), Pts(:,2), Pts(:,3), 'Facecolor','blue');



  figure; trisurf(Tri(face_gt == 2,:), Pts(:,1), Pts(:,2), Pts(:,3), 'Facecolor','red');
 hold on; trisurf(Tri(face_gt == 1,:), Pts(:,1), Pts(:,2), Pts(:,3), 'Facecolor','blue');