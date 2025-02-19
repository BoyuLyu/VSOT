%% main functions of VSOT
%%==========================================================================
%Build the mex functions
addpath(genpath('./src/'));
build_mex('./resources/');
addpath('./resources/curvatures/')
addpath('./resources/ImprovedSurfaceSmooth/')
addpath('./resources/CurvatureEstimation')
addpath('./resources/edt_mex/edt_mex')
addpath('./resources/TAUBIN/TAUBIN/')
addpath('./resources/iso2mesh/')
addpath('./resources/data_to_off/')
addpath('./resources/src_mex/mex_EM_analysis/mex_EM_analysis')
addpath('./resources/tinymesh_mex/tiny_mesh_mex')
addpath('./resources/graph_related/graph_mex/')
addpath('./resources/sigstar-master/')
addpath('./resources/Violinplot-Matlab/')
addpath(genpath('./src/'))



addpath('./resources/curvatures/')
addpath('./resources/ImprovedSurfaceSmooth/')
addpath('./resources/CurvatureEstimation')
addpath('/home/boyu/Documents/edt_mex/edt_mex/edt_mex')
addpath('./resources/TAUBIN/TAUBIN/')
addpath('/home/boyu/Documents/tinymesh_mex/tiny_mesh_mex')
addpath('/home/boyu/Documents/iso2mesh/')
addpath('./resources/data_to_off/')
addpath('./resources/data_to_off/')
addpath('/home/boyu/Documents/graph_related/graph_mex/')

disp('Build mex functions finished');
data_path = '/work/boyu/EM_astrocyte/materials_for_paper_VSOT/VSOT_data_paired_with_codes/data/';
%%==========================================================================
% Run the whole pipeline of VSOT for one neuron
save_main_path = '/work/boyu/EM_astrocyte/materials_for_paper_VSOT/VSOT_data_paired_with_codes/results/';
main_whole_pipeline_one_neuron_code_ocean(data_path, save_main_path);
disp('whole pipeline finished running');
%%==========================================================================
% Run VSOT for level-1 segmentation and test the performance
save_path_1 = '/work/boyu/EM_astrocyte/materials_for_paper_VSOT/VSOT_data_paired_with_codes/results/test_level_1_segmentation/';
if(~exist(save_path_1, 'dir'))
    mkdir(save_path_1);
end
tic;
test_level_1_segmentation(data_path, save_path_1);
disp('Level-1 segmentation performance testing finished');
toc;
%%==========================================================================
% Run VSOT for level-2 segmentation and test the performance
save_path_2 = '/work/boyu/EM_astrocyte/materials_for_paper_VSOT/VSOT_data_paired_with_codes/results/test_level_2_segmentation/';
if(~exist(save_path_2, 'dir'))
    mkdir(save_path_2);
end
test_level_2_segmentation(data_path,save_path_2);
disp('Level-2 segmentation performance testing finished');

%%==========================================================================
% plot the quantification results of structural analysis
save_path_3 = '/work/boyu/EM_astrocyte/materials_for_paper_VSOT/VSOT_data_paired_with_codes/results/quantification_plot/';
if(~exist(save_path_3, 'dir'))
    mkdir(save_path_3);
end
main_plot_quantification(data_path,save_path_3);
disp('example quantification plots finished');

%%==========================================================================