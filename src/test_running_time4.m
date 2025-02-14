curvature_scale = 30;
resx_s = 2;
resy_s = 2;
resz_s = 5;
data_path = 'D:\Guoqiang_lab\VSOT_samples\VSOT_data_paired_with_codes\data\performance_testing_level_2_segmentation\large_dataset_3_annotator_L2345_w_stubby\gt_final_800\uniformly_selected_samples';
result_path = 'D:\Guoqiang_lab\VSOT_samples\VSOT_data_paired_with_codes\data\performance_testing_level_2_segmentation\large_dataset_3_annotator_L2345_w_stubby\segmentation_results';
spine_save_folder = fullfile(data_path, "tif_file\");
offFolder = fullfile(data_path, 'all/');
spine_head_neck_save_folder = fullfile(result_path, "our_result_all_tmp_test\volume_segmentation200/");
if(~exist(spine_head_neck_save_folder, 'dir'))
    mkdir(spine_head_neck_save_folder)
end
coordinate_output_folder = fullfile(result_path, "our_result_all_tmp_test\head_neck_coor_result200/");
if(~exist(coordinate_output_folder, 'dir'))
    mkdir(coordinate_output_folder)
end
our_method_cut_result_folder = fullfile(result_path, "our_result_all_tmp_test\cut_result200/");
if(~exist(our_method_cut_result_folder, 'dir'))
    mkdir(our_method_cut_result_folder);
end
disp(200);
tic;
comSeg.level_2_segmentation_main_speedup_output_cutVex_v2(spine_save_folder, spine_head_neck_save_folder, coordinate_output_folder, curvature_scale, resx_s, resy_s, resz_s,200,1);
toc;
segmentation_our_own_method_from_volume(offFolder, spine_head_neck_save_folder, our_method_cut_result_folder, spine_save_folder);