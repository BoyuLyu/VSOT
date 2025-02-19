function test_level_2_segmentation_time(data_path, result_path)
offFolder = fullfile(data_path,'off_file');
tifFolder = fullfile(data_path, 'tif_file');
annotationFolder = fullfile(data_path,'/performance_testing_level_2_segmentation/gt_300/annotation_json_300');
annotationFolder_gt_curves = fullfile(data_path,'refined_cutID');
listx = dir([offFolder, '/*.off']);
tamada_result_folder = fullfile(result_path,'tamada_result');
our_method_cut_result_folder = fullfile(result_path,'our_result');
ofer_output_folder = fullfile(result_path,'ofer_result');
dorkenwald_result_folder = fullfile(result_path,'dorkenwald_result');
if(~exist(tamada_result_folder, 'dir'))
    mkdir(tamada_result_folder);
end
if(~exist(our_method_cut_result_folder, 'dir'))
    mkdir(our_method_cut_result_folder);
end
if(~exist(ofer_output_folder, 'dir'))
    mkdir(ofer_output_folder);
end
if(~exist(dorkenwald_result_folder, 'dir'))
    mkdir(dorkenwald_result_folder);
end
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

profile on 
data_path = '/work/boyu/EM_astrocyte/materials_for_paper_VSOT/VSOT_data_paired_with_codes/data/performance_testing_level_2_segmentation/large_dataset_3_annotator_L2345_w_stubby/gt_final_800/subtypes';
result_path = '/work/boyu/EM_astrocyte/materials_for_paper_VSOT/VSOT_data_paired_with_codes/data/performance_testing_level_2_segmentation/large_dataset_3_annotator_L2345_w_stubby/segmentation_results';
spine_save_folder = fullfile(data_path, "all_tif/");
spine_head_neck_save_folder = fullfile(result_path, "our_result_all/volume_segmentation/");
if(~exist(spine_head_neck_save_folder, 'dir'))
    mkdir(spine_head_neck_save_folder)
end
coordinate_output_folder = fullfile(result_path, "our_result_all/head_neck_coor_result/");
if(~exist(coordinate_output_folder, 'dir'))
    mkdir(coordinate_output_folder)
end
our_method_cut_result_folder = fullfile(result_path, "our_result_all/cut_result/");
if(~exist(our_method_cut_result_folder, 'dir'))
    mkdir(our_method_cut_result_folder);
end
comSeg.level_2_segmentation_main_speedup_output_cutVex_v2(spine_save_folder, spine_head_neck_save_folder, coordinate_output_folder, curvature_scale, resx_s, resy_s, resz_s,1200,1);
segmentation_our_own_method(offFolder, coordinate_output_folder, our_method_cut_result_folder, spine_save_folder);



data_path = 'D:\Guoqiang_lab\VSOT_samples\VSOT_data_paired_with_codes\data\performance_testing_level_2_segmentation\large_dataset_3_annotator_L2345_w_stubby\gt_final_800\uniformly_selected_samples';
result_path = 'D:\Guoqiang_lab\VSOT_samples\VSOT_data_paired_with_codes\data\performance_testing_level_2_segmentation\large_dataset_3_annotator_L2345_w_stubby\segmentation_results';
spine_save_folder = fullfile(data_path, "tif_file\");
offFolder = fullfile(data_path, 'all/');
spine_head_neck_save_folder = fullfile(result_path, "our_result_all_tmp_test\volume_segmentation1200/");
if(~exist(spine_head_neck_save_folder, 'dir'))
    mkdir(spine_head_neck_save_folder)
end
coordinate_output_folder = fullfile(result_path, "our_result_all_tmp_test\head_neck_coor_result1200/");
if(~exist(coordinate_output_folder, 'dir'))
    mkdir(coordinate_output_folder)
end
our_method_cut_result_folder = fullfile(result_path, "our_result_all_tmp_test\cut_result1200/");
if(~exist(our_method_cut_result_folder, 'dir'))
    mkdir(our_method_cut_result_folder);
end
disp(1200);
tic;
comSeg.level_2_segmentation_main_speedup_output_cutVex_v2(spine_save_folder, spine_head_neck_save_folder, coordinate_output_folder, curvature_scale, resx_s, resy_s, resz_s,1200,1);
toc;
segmentation_our_own_method_from_volume(offFolder, spine_head_neck_save_folder, our_method_cut_result_folder, spine_save_folder);




data_path = 'D:\Guoqiang_lab\VSOT_samples\VSOT_data_paired_with_codes\data\performance_testing_level_2_segmentation\large_dataset_3_annotator_L2345_w_stubby\gt_final_800\uniformly_selected_samples';
result_path = 'D:\Guoqiang_lab\VSOT_samples\VSOT_data_paired_with_codes\data\performance_testing_level_2_segmentation\large_dataset_3_annotator_L2345_w_stubby\segmentation_results';
spine_save_folder = fullfile(data_path, "tif_file\");
offFolder = fullfile(data_path, 'all/');
spine_head_neck_save_folder = fullfile(result_path, "our_result_all_tmp_test\volume_segmentation800/");
if(~exist(spine_head_neck_save_folder, 'dir'))
    mkdir(spine_head_neck_save_folder)
end
coordinate_output_folder = fullfile(result_path, "our_result_all_tmp_test\head_neck_coor_result800/");
if(~exist(coordinate_output_folder, 'dir'))
    mkdir(coordinate_output_folder)
end
our_method_cut_result_folder = fullfile(result_path, "our_result_all_tmp_test\cut_result800/");
if(~exist(our_method_cut_result_folder, 'dir'))
    mkdir(our_method_cut_result_folder);
end
disp(800);
tic;
comSeg.level_2_segmentation_main_speedup_output_cutVex_v2(spine_save_folder, spine_head_neck_save_folder, coordinate_output_folder, curvature_scale, resx_s, resy_s, resz_s,800,1);
toc;
vsot_save_face_label(offFolder, spine_head_neck_save_folder, our_method_cut_result_folder, spine_save_folder);




data_path = 'D:\Guoqiang_lab\VSOT_samples\VSOT_data_paired_with_codes\data\performance_testing_level_2_segmentation\large_dataset_3_annotator_L2345_w_stubby\gt_final_800\uniformly_selected_samples';
result_path = 'D:\Guoqiang_lab\VSOT_samples\VSOT_data_paired_with_codes\data\performance_testing_level_2_segmentation\large_dataset_3_annotator_L2345_w_stubby\segmentation_results';
spine_save_folder = fullfile(data_path, "tif_file\");
offFolder = fullfile(data_path, 'all/');
spine_head_neck_save_folder = fullfile(result_path, "our_result_all_tmp_test\volume_segmentation400/");
if(~exist(spine_head_neck_save_folder, 'dir'))
    mkdir(spine_head_neck_save_folder)
end
coordinate_output_folder = fullfile(result_path, "our_result_all_tmp_test\head_neck_coor_result400/");
if(~exist(coordinate_output_folder, 'dir'))
    mkdir(coordinate_output_folder)
end
our_method_cut_result_folder = fullfile(result_path, "our_result_all_tmp_test\cut_result400/");
if(~exist(our_method_cut_result_folder, 'dir'))
    mkdir(our_method_cut_result_folder);
end
disp(400);
tic;
comSeg.level_2_segmentation_main_speedup_output_cutVex_v2(spine_save_folder, spine_head_neck_save_folder, coordinate_output_folder, curvature_scale, resx_s, resy_s, resz_s,400,1);
toc;
segmentation_our_own_method_from_volume(offFolder, spine_head_neck_save_folder, our_method_cut_result_folder, spine_save_folder);




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


data_path = 'D:\Guoqiang_lab\VSOT_samples\VSOT_data_paired_with_codes\data\performance_testing_level_2_segmentation\large_dataset_3_annotator_L2345_w_stubby\gt_final_800\uniformly_selected_samples';
result_path = 'D:\Guoqiang_lab\VSOT_samples\VSOT_data_paired_with_codes\data\performance_testing_level_2_segmentation\large_dataset_3_annotator_L2345_w_stubby\segmentation_results';
spine_save_folder = fullfile(data_path, "tif_file\");
offFolder = fullfile(data_path, 'all/');
spine_head_neck_save_folder = fullfile(result_path, "our_result_all_tmp_test\volume_segmentation100/");
if(~exist(spine_head_neck_save_folder, 'dir'))
    mkdir(spine_head_neck_save_folder)
end
coordinate_output_folder = fullfile(result_path, "our_result_all_tmp_test\head_neck_coor_result100/");
if(~exist(coordinate_output_folder, 'dir'))
    mkdir(coordinate_output_folder)
end
our_method_cut_result_folder = fullfile(result_path, "our_result_all_tmp_test\cut_result100/");
if(~exist(our_method_cut_result_folder, 'dir'))
    mkdir(our_method_cut_result_folder);
end
disp(100);
tic;
comSeg.level_2_segmentation_main_speedup_output_cutVex_v2(spine_save_folder, spine_head_neck_save_folder, coordinate_output_folder, curvature_scale, resx_s, resy_s, resz_s,100,1);
toc;
segmentation_our_own_method_from_volume(offFolder, spine_head_neck_save_folder, our_method_cut_result_folder, spine_save_folder);





profile viewer
profile off
% match the segmentation results to the surface label to fairly compare with peer methods
% head_coordinates_match_to_gt_surface(offFolder, annotationFolder, annotationFolder_gt_curves,coordinate_output_folder,our_method_cut_result_folder);
% preprocessing on the ground truth data
%% run peer methods


main_segmentation_error_testing(offFolder,annotationFolder, annotationFolder_gt_curves,tamada_result_folder,our_method_cut_result_folder,ofer_output_folder,dorkenwalk_result_folder,save_path);



%ofer result
segmentation_peer_method_ofer(offFolder, ofer_output_folder, tifFolder)

data_path = '/work/boyu/EM_astrocyte/materials_for_paper_VSOT/VSOT_data_paired_with_codes/data/performance_testing_level_2_segmentation/large_dataset_3_annotator_L2345_w_stubby/gt_final_800/uniformly_selected_samples';
root_folder = data_path;
root_folder = '/work/boyu/EM_astrocyte/materials_for_paper_VSOT/VSOT_data_paired_with_codes/data/performance_testing_level_2_segmentation/large_dataset_3_annotator_L2345_w_stubby/gt_final_800/stubby_samples';
result_folder = '/work/boyu/EM_astrocyte/materials_for_paper_VSOT/VSOT_data_paired_with_codes/data/performance_testing_level_2_segmentation/large_dataset_3_annotator_L2345_w_stubby/segmentation_results';
surfaceFolder = fullfile(root_folder, 'off_file');
volumeGTOutfolder = fullfile(root_folder, 'tif_file');
skelfolder = fullfile(root_folder, 'spine_skel');
dorkenwalk_result_folder = fullfile(result_folder, 'dorkenwald_result/');
segmentation_peer_method_dorkenwald(surfaceFolder, volumeGTOutfolder,skelfolder,dorkenwalk_result_folder)

root_folder = '/work/boyu/EM_astrocyte/materials_for_paper_VSOT/VSOT_data_paired_with_codes/data/performance_testing_level_2_segmentation/large_dataset_3_annotator_L2345_w_stubby/gt_final_800/uniformly_selected_samples';
result_folder = '/work/boyu/EM_astrocyte/materials_for_paper_VSOT/VSOT_data_paired_with_codes/data/performance_testing_level_2_segmentation/large_dataset_3_annotator_L2345_w_stubby/segmentation_results';
offFolder = fullfile(root_folder, 'off_file');
% listx = dir(fullfile(offFolder, '*.off'));
% source_coordinate_folder = '/work/boyu/EM_astrocyte/materials_for_paper_VSOT/VSOT_data_paired_with_codes/data/performance_testing_level_2_segmentation/large_dataset_3_annotator_L2345_w_stubby/segmentation_results/our_result/head_neck_coor_result';
% listx_name = {listx.name};
our_method_coordinate_folder = fullfile(result_folder, 'our_result/head_neck_coor_result');
our_method_cut_result_folder = fullfile(result_folder, 'vsot_result_total');
% if(~exist(our_method_coordinate_folder, 'dir'))
%     mkdir(our_method_coordinate_folder);
% end
% if(~exist(our_method_cut_result_folder, 'dir'))
%     mkdir(our_method_cut_result_folder);
% end
% for i = 1:length(listx_name)
%     tmp_name = listx_name{i};
%     tmp_name = strsplit(tmp_name, '.');
%     namex = tmp_name{1};
%     if(exist(fullfile(source_coordinate_folder,[namex,'_head.mat']), 'file'))
%         copyfile(fullfile(source_coordinate_folder, [namex,'_head.mat']), fullfile(our_method_coordinate_folder,[namex,'_head.mat']));
%         copyfile(fullfile(source_coordinate_folder, [namex,'_neck.mat']), fullfile(our_method_coordinate_folder,[namex,'_neck.mat']));
%     end
% end
segmentation_our_own_method(offFolder, our_method_coordinate_folder, our_method_cut_result_folder)

listx = dir(fullfile(surfaceFolder, '*.off'));
listx_name = {listx.name};
for i = 1:length(listx_name)
    namex =listx_name{i};
    namex = strsplit(namex, '.');
    namex = namex{1};
    if(~exist(fullfile(volumeGTOutfolder, [namex, '.tif']), 'file'))
        delete(fullfile(surfaceFolder,[namex, '.off']));
    end
end


tif_folder = fullfile(data_path, 'tif_file/');
tamada_result_folder = fullfile(result_folder, 'tamada_result_stubby/');
segmentation_peer_method_tamada(offFolder, tif_folder,tamada_result_folder)



data_path = '/work/boyu/EM_astrocyte/materials_for_paper_VSOT/VSOT_data_paired_with_codes/data/performance_testing_level_2_segmentation/large_dataset_3_annotator_L2345_w_stubby/gt_final_800/subtypes';
    
result_path = '/work/boyu/EM_astrocyte/materials_for_paper_VSOT/VSOT_data_paired_with_codes/data/performance_testing_level_2_segmentation/large_dataset_3_annotator_L2345_w_stubby/segmentation_results';

offFolder = fullfile(data_path, 'all/');
annotationFolder_gt_curves = fullfile(data_path, 'cut_index');
tamada_result_folder = fullfile(result_path, 'tamada_result_total/');
% our_method_cut_result_folder = fullfile(result_path, 'vsot_result_total');
our_method_cut_result_folder = fullfile(result_path, 'our_result_all/cut_result/');
ofer_output_folder = fullfile(result_path, 'ofer_total/');
dorkenwalk_result_folder = fullfile(result_path,'dorkenwald_total/' );
out1 = main_segmentation_error_testing_v3(offFolder,annotationFolder_gt_curves,tamada_result_folder,our_method_cut_result_folder,ofer_output_folder,dorkenwalk_result_folder,[]);

tamada_error_f1 = out1.tamada_error_f1;
tamada_error_precision = out1.tamada_error_precision;
tamada_error_recall = out1.tamada_error_recall;
tamada_error_iou = out1.tamada_error_iou;
our_method_error_f1 = out1.our_method_error_f1;
our_method_error_precision = out1.our_method_error_precision;
our_method_error_recall = out1.our_method_error_recall;
our_method_error_iou = out1.our_method_error_iou;
ofer_method_error_f1 = out1.ofer_method_error_f1;
ofer_method_error_precision = out1.ofer_method_error_precision;
ofer_method_error_recall = out1.ofer_method_error_recall;
ofer_method_error_iou = out1.ofer_method_error_iou;
dorkenwalk_error_f1 = out1.dorkenwalk_error_f1;
dorkenwalk_error_precision = out1.dorkenwalk_error_precision;
dorkenwalk_error_recall = out1.dorkenwalk_error_recall;
dorkenwalk_error_iou = out1.dorkenwalk_error_iou;

nan_rm1 = find(isnan(tamada_error_f1)| isnan(our_method_error_f1) |isnan(ofer_method_error_f1)|isnan(dorkenwalk_error_f1));
nan_rm2 = find(isnan(tamada_error_precision)| isnan(our_method_error_precision) |isnan(ofer_method_error_precision)|isnan(dorkenwalk_error_precision));
nan_rm3 = find(isnan(tamada_error_recall)| isnan(our_method_error_recall) |isnan(ofer_method_error_recall)|isnan(dorkenwalk_error_recall));
nan_rm4 = find(isnan(tamada_error_iou)| isnan(our_method_error_iou) |isnan(ofer_method_error_iou)|isnan(dorkenwalk_error_iou));

tamada_error_f1(nan_rm1) = [];
our_method_error_f1(nan_rm1) = [];
ofer_method_error_f1(nan_rm1) = [];
dorkenwalk_error_f1(nan_rm1) = [];
tamada_error_precision(nan_rm2) = [];
our_method_error_precision(nan_rm2) = [];
ofer_method_error_precision(nan_rm2) = [];
dorkenwalk_error_precision(nan_rm2) = [];
tamada_error_recall(nan_rm3) = [];
our_method_error_recall(nan_rm3) = [];
ofer_method_error_recall(nan_rm3) = [];
dorkenwalk_error_recall(nan_rm3) = [];
tamada_error_iou(nan_rm4) = [];
our_method_error_iou(nan_rm4) = [];
ofer_method_error_iou(nan_rm4) = [];
dorkenwalk_error_iou(nan_rm4) = [];


fprintf('f1: VSOT %.2f, O.N method %.2f,  T.H method %.2f ,D.S method %.2f \n' ...
    , mean(our_method_error_f1), mean(ofer_method_error_f1), mean(tamada_error_f1), mean(dorkenwalk_error_f1));
fprintf('precision: VSOT %.2f, O.N method %.2f,  T.H method %.2f ,D.S method %.2f \n' ...
    , mean(our_method_error_precision), mean(ofer_method_error_precision), mean(tamada_error_precision), mean(dorkenwalk_error_precision));
fprintf('recall: VSOT %.2f, O.N method %.2f,  T.H method %.2f ,D.S method %.2f \n' ...
    , mean(our_method_error_recall), mean(ofer_method_error_recall), mean(tamada_error_recall), mean(dorkenwalk_error_recall));
fprintf('iou: VSOT %.2f, O.N method %.2f,  T.H method %.2f ,D.S method %.2f \n' ...
    , mean(our_method_error_iou), mean(ofer_method_error_iou), mean(tamada_error_iou), mean(dorkenwalk_error_iou));



fprintf('f1 std: VSOT %.2f, O.N method %.2f,  T.H method %.2f ,D.S method %.2f \n' ...
    , std(our_method_error_f1), std(ofer_method_error_f1), std(tamada_error_f1), std(dorkenwalk_error_f1));
fprintf('precision std: VSOT %.2f, O.N method %.2f,  T.H method %.2f ,D.S method %.2f \n' ...
    , std(our_method_error_precision), std(ofer_method_error_precision), std(tamada_error_precision), std(dorkenwalk_error_precision));
fprintf('recall std: VSOT %.2f, O.N method %.2f,  T.H method %.2f ,D.S method %.2f \n' ...
    , std(our_method_error_recall), std(ofer_method_error_recall), std(tamada_error_recall), std(dorkenwalk_error_recall));
fprintf('iou std: VSOT %.2f, O.N method %.2f,  T.H method %.2f ,D.S method %.2f \n' ...
    , std(our_method_error_iou), std(ofer_method_error_iou), std(tamada_error_iou), std(dorkenwalk_error_iou));


% fprintf('area error: VSOT %f, O.N method %f, D.S method %f, T.H method %f \n', ...
%     median(our_method_error_all1), median(ofer_method_error_all1), median(dorkenwalk_error_all1), median(tamada_error_all1))
% fprintf('area error dispersion: VSOT %f, O.N method %f, D.S method %f, T.H method %f \n', ...
%     std(our_method_error_all1), std(ofer_method_error_all1), std(dorkenwalk_error_all1), std(tamada_error_all1))

% fprintf('dist error: VSOT %f, O.N method %f, D.S method %f, T.H method %f \n', ...
%     median(our_method_error_all_dist1), median(ofer_method_error_all_dist1), median(dorkenwalk_error_all_dist1), median(tamada_error_all_dist1))

% fprintf('dist error dispersion: VSOT %f, O.N method %f, D.S method %f, T.H method %f \n', ...
%     std(our_method_error_all_dist1), std(ofer_method_error_all_dist1), std(dorkenwalk_error_all_dist1), std(tamada_error_all_dist1))

% fprintf('iou : VSOT %f, O.N method %f, D.S method %f, T.H method %f \n', ...
%     median(1 - our_method_error_iou1), median(1 - ofer_method_error_iou1), median(1 - dorkenwalk_error_iou1), median(1 - tamada_error_all_iou1))

% fprintf('iou std: VSOT %f, O.N method %f, D.S method %f, T.H method %f \n', ...
%     std(1 - our_method_error_iou1), std(1 - ofer_method_error_iou1), std(1 - dorkenwalk_error_iou1), std(1 - tamada_error_all_iou1))

% figure; boxplot([our_method_error_all1;ofer_method_error_all1;dorkenwalk_error_all1;tamada_error_all1], [ones(size(our_method_error_all1));...
%     ones(size(our_method_error_all1))*2 ;ones(size(our_method_error_all1))*3;ones(size(our_method_error_all1))*4])



data_path = '/work/boyu/EM_astrocyte/materials_for_paper_VSOT/VSOT_data_paired_with_codes/data/performance_testing_level_2_segmentation/large_dataset_3_annotator_L2345_w_stubby/gt_final_800/stubby_samples';

result_path = '/work/boyu/EM_astrocyte/materials_for_paper_VSOT/VSOT_data_paired_with_codes/data/performance_testing_level_2_segmentation/large_dataset_3_annotator_L2345_w_stubby/segmentation_results';

offFolder = fullfile(data_path, 'off_file/');
annotationFolder_gt_curves = fullfile(data_path, 'cut_index');
tamada_result_folder = fullfile(result_path, 'tamada_result_stubby/');
our_method_cut_result_folder = fullfile(result_path, 'our_result_stubby/cut_result/');
ofer_output_folder = fullfile(result_path, 'ofer_result_stubby/');
dorkenwalk_result_folder = fullfile(result_path,'dorkenwald_result_stubby/' );
out2 = main_segmentation_error_testing_v2(offFolder,annotationFolder_gt_curves,tamada_result_folder,our_method_cut_result_folder,ofer_output_folder,dorkenwalk_result_folder,[]);

tamada_error_all2 = out2.tamada_error_all;
tamada_error_all_dist2 = out2.tamada_error_all_dist;
tamada_error_iou2 = out2.tamada_error_iou;
our_method_error_all2 = out2.our_method_error_all;
our_method_error_all_dist2 = out2.our_method_error_all_dist ;
our_method_error_iou2 = out2.our_method_error_iou;
ofer_method_error_all2 = out2.ofer_method_error_all;
ofer_method_error_all_dist2  =    out2.ofer_method_error_all_dist;
ofer_method_error_iou2 = out2.ofer_method_error_iou;
dorkenwalk_error_all2 = out2.dorkenwalk_error_all;
dorkenwalk_error_all_dist2=  out2.dorkenwalk_error_all_dist;
dorkenwalk_error_iou2 = out2.dorkenwalk_error_iou;

nan_rm1 = find(isnan(tamada_error_all2)| isnan(our_method_error_all2) |isnan(ofer_method_error_all2)|isnan(dorkenwalk_error_all2));
nan_rm2 = find(isnan(tamada_error_all_dist2)| isnan(our_method_error_all_dist2) |isnan(ofer_method_error_all_dist2)|isnan(dorkenwalk_error_all_dist2));
nan_rm3 = find(isnan(tamada_error_iou2)| isnan(our_method_error_iou2) |isnan(ofer_method_error_iou2)|isnan(dorkenwalk_error_iou2));
tamada_error_all2(nan_rm1) = [];
our_method_error_all2(nan_rm1) = [];
ofer_method_error_all2(nan_rm1) = [];
dorkenwalk_error_all2(nan_rm1) = [];
tamada_error_all_dist2(nan_rm2) = [];
our_method_error_all_dist2(nan_rm2) = [];
ofer_method_error_all_dist2(nan_rm2) = [];
dorkenwalk_error_all_dist2(nan_rm2) = [];
tamada_error_iou2(nan_rm3) = [];
our_method_error_iou2(nan_rm3) = [];
ofer_method_error_iou2(nan_rm3) = [];
dorkenwalk_error_iou2(nan_rm3) = [];
fprintf('area error: VSOT %f, O.N method %f, D.S method %f, T.H method %f \n', ...
    mean(our_method_error_all2), mean(ofer_method_error_all2), mean(dorkenwalk_error_all2), mean(tamada_error_all2))

fprintf('dist error: VSOT %f, O.N method %f, D.S method %f, T.H method %f \n', ...
    mean(our_method_error_all_dist2), mean(ofer_method_error_all_dist2), mean(dorkenwalk_error_all_dist2), mean(tamada_error_all_dist2))

fprintf('IOU: VSOT %f, O.N method %f, D.S method %f, T.H method %f \n', ...
    mean(our_method_error_iou2), mean(ofer_method_error_iou2), mean(dorkenwalk_error_iou2), mean(tamada_error_iou2))





fprintf('area error: VSOT %f, O.N method %f, D.S method %f, T.H method %f \n', ...
    median([our_method_error_all1;our_method_error_all2]), median([ofer_method_error_all1;ofer_method_error_all2]), median([dorkenwalk_error_all1;dorkenwalk_error_all2 ]), median([tamada_error_all1;tamada_error_all2]))

fprintf('dist error: VSOT %f, O.N method %f, D.S method %f, T.H method %f \n', ...
    median([our_method_error_all_dist1;our_method_error_all_dist2]), median([ofer_method_error_all_dist1;ofer_method_error_all_dist2]), median([dorkenwalk_error_all_dist1;dorkenwalk_error_all_dist2 ]), median([tamada_error_all_dist1;tamada_error_all_dist2]))


data_path = '/work/boyu/EM_astrocyte/materials_for_paper_VSOT/VSOT_data_paired_with_codes/data/performance_testing_level_2_segmentation/large_dataset_3_annotator_L2345_w_stubby/gt_final_800/uniformly_selected_samples';

result_path = '/work/boyu/EM_astrocyte/materials_for_paper_VSOT/VSOT_data_paired_with_codes/data/performance_testing_level_2_segmentation/large_dataset_3_annotator_L2345_w_stubby/segmentation_results';

offFolder = fullfile(data_path, 'off_file/');
tifFolder = fullfile(data_path, 'tif_file/');
cutFolder = fullfile(data_path, 'refined_cutID/');
gen_surface_cut_gt(offFolder, tifFolder, cutFolder);



data_path = '/work/boyu/EM_astrocyte/materials_for_paper_VSOT/VSOT_data_paired_with_codes/data/performance_testing_level_2_segmentation/large_dataset_3_annotator_L2345_w_stubby/gt_final_800/stubby_samples';
offFolder = fullfile(data_path, 'off_file/');
tifFolder = fullfile(data_path, 'tif_file/');
cutFolder = fullfile(data_path, 'cut_index/');
gen_surface_cut_gt(offFolder, tifFolder, cutFolder);
