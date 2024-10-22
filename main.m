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
addpath(genpath('./src/'))
disp('Build mex functions finished');
data_path = '/data/';
%%==========================================================================
% Run the whole pipeline of VSOT for one neuron
save_main_path = '/results/';
main_whole_pipeline_one_neuron_code_ocean(data_path, save_main_path);
disp('whole pipeline finished running');
%%==========================================================================
% Run VSOT for level-1 segmentation and test the performance
save_path_1 = '/results/test_level_1_segmentation/';
if(~exist(save_path_1, 'dir'))
    mkdir(save_path_1);
end
test_level_1_segmentation(data_path, save_path_1);
disp('Level-1 segmentation performance testing finished');
%%==========================================================================
% Run VSOT for level-2 segmentation and test the performance
save_path_2 = '/results/test_level_2_segmentation/';
if(~exist(save_path_2, 'dir'))
    mkdir(save_path_2);
end
test_level_2_segmentation(data_path,save_path_2);
disp('Level-2 segmentation performance testing finished');

%%==========================================================================
% plot the quantification results of structural analysis
save_path_3 = '/results/quantification_plot/';
main_plot_quantification(data_path);
disp('example quantification plots finished');

%%==========================================================================