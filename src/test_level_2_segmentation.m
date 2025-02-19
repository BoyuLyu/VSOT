function test_level_2_segmentation(data_path, result_path)
    curvature_scale = 30;
    resx_s = 2;
    resy_s = 2;
    resz_s = 5;
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
    offFolder = fullfile(data_path, 'all/');
    comSeg.level_2_segmentation_main_speedup_output_cutVex_v2(spine_save_folder, spine_head_neck_save_folder, coordinate_output_folder, curvature_scale, resx_s, resy_s, resz_s,1200,1);
    vsot_save_face_label(offFolder, coordinate_output_folder, our_method_cut_result_folder, spine_save_folder);

    % compare the results with the peer method
    offFolder = fullfile(data_path, 'all/');
    annotationFolder_gt_curves = fullfile(data_path, 'cut_index');
    tamada_result_folder = fullfile(result_path, 'tamada_result_total/');
    % our_method_cut_result_folder = fullfile(result_path, 'vsot_result_total');
    our_method_cut_result_folder = fullfile(result_path, 'our_result_all/cut_result/');
    ofer_output_folder = fullfile(result_path, 'ofer_total/');
    dorkenwalk_result_folder = fullfile(result_path,'dorkenwald_total/' );
    out1 = main_segmentation_error_testing_v3(offFolder,annotationFolder_gt_curves,tamada_result_folder,our_method_cut_result_folder,ofer_output_folder,dorkenwalk_result_folder,[]);
    [tamada_error_f1, our_method_error_f1, ofer_method_error_f1, dorkenwalk_error_f1,...
    tamada_error_precision, our_method_error_precision,ofer_method_error_precision,dorkenwalk_error_precision,...
        tamada_error_recall,our_method_error_recall,ofer_method_error_recall,dorkenwalk_error_recall,...
            tamada_error_iou,our_method_error_iou,ofer_method_error_iou, dorkenwalk_error_iou] = measure_performance(out1);
    
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


end
