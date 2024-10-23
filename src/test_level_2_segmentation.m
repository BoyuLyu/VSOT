function test_level_2_segmentation(data_path, save_path)
offFolder = fullfile(data_path,'/performance_testing_level_2_segmentation/gt_300/surface_off_300');
annotationFolder = fullfile(data_path,'/performance_testing_level_2_segmentation/gt_300/annotation_json_300');
annotationFolder_gt_curves = fullfile(data_path,'/performance_testing_level_2_segmentation/gt_300/cut_cycle_ID_300');
listx = dir([offFolder, '/*.off']);
tamada_result_folder = fullfile(data_path,'/performance_testing_level_2_segmentation/dendrite_segmentation_peer_methods/tamada/results');
our_method_cut_result_folder = fullfile(save_path,'/performance_testing_level_2_segmentation/dendrite_segmentation_peer_methods/our_method/result');
ofer_output_folder = fullfile(data_path,'/performance_testing_level_2_segmentation/dendrite_segmentation_peer_methods/ofer/result');
dorkenwalk_result_folder = fullfile(data_path,'/performance_testing_level_2_segmentation/dendrite_segmentation_peer_methods/dorkenwald/result');
spine_save_folder = fullfile(data_path,'/performance_testing_level_2_segmentation/gt_300/volume_300_w_shaft');
spine_head_neck_save_folder = fullfile(save_path,'/performance_testing_level_2_segmentation/dendrite_segmentation_peer_methods/our_method/volume_result');
coordinate_output_folder = fullfile(save_path,'/performance_testing_level_2_segmentation/dendrite_segmentation_peer_methods/our_method/coor_result');

if(~exist(spine_head_neck_save_folder, 'dir'))
    mkdir(spine_head_neck_save_folder);
end

if(~exist(coordinate_output_folder, 'dir'))
    mkdir(coordinate_output_folder);
end

curvature_scale = 30;
resx_s = 2;
resy_s = 2;
resz_s = 5;
if(~exist(spine_head_neck_save_folder,"dir"))
    mkdir(spine_head_neck_save_folder);
else
    rmdir(spine_head_neck_save_folder,'s');
    mkdir(spine_head_neck_save_folder);
end

if(~exist(coordinate_output_folder,"dir"))
    mkdir(coordinate_output_folder);
else
    rmdir(coordinate_output_folder,'s');
    mkdir(coordinate_output_folder);
end

if(~exist(our_method_cut_result_folder,"dir"))
    mkdir(our_method_cut_result_folder);
else
    rmdir(our_method_cut_result_folder,'s');
    mkdir(our_method_cut_result_folder);
end

comSeg.level_2_segmentation_main_speedup_output_cutVex(spine_save_folder, spine_head_neck_save_folder, coordinate_output_folder, curvature_scale, resx_s, resy_s, resz_s,1200,1);
% match the segmentation results to the surface label to fairly compare with peer methods
head_coordinates_match_to_gt_surface(offFolder, annotationFolder, annotationFolder_gt_curves,coordinate_output_folder,our_method_cut_result_folder);

main_segmentation_error_testing(offFolder,annotationFolder, annotationFolder_gt_curves,tamada_result_folder,our_method_cut_result_folder,ofer_output_folder,dorkenwalk_result_folder,save_path);