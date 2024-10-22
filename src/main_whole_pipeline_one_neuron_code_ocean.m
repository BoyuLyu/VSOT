function main_whole_pipeline_one_neuron_code_ocean(data_path, save_main_path)
%% VSOT: 
% showcase the whole pipeline 
% input: A binary segmentation of a neuron, the resolution of the structure
% output: 
% (1) Segmentation of the dendrite into dendrite spine and dendrite shaft
%    Stored in a folder containing all the spines and a part of dendrite shaft that is attached to the spine (spine: 2, shaft: 1)
% (2) segmentation of dendrite spine into dendrite spine head and neck
%    Stored in a folder containing all the spine head and neck regions (spinne head: 3, spine neck: 2, dendrite shaft: 1)
% (3) Quantification of the morphological features of dendrite spine
% and astrocyte related quantification


%%==========================================================================
% Setup the paths & environments
% Generate skeleton and then remove soma to obtain individual dendritic
% trunks
% addpath('./resources/curvatures/')
% addpath('./resources/ImprovedSurfaceSmooth/')
% addpath('./resources/CurvatureEstimation')
% addpath('./resources/edt_mex/edt_mex')
% addpath('./resources/TAUBIN/TAUBIN/')
% addpath('./resources/iso2mesh/')
% addpath('./resources/data_to_off/')
% addpath('./resources/src_mex/mex_EM_analysis/mex_EM_analysis')
% addpath('./resources/tinymesh_mex/tiny_mesh_mex')
% addpath('./resources/graph_related/graph_mex/')
% addpath(genpath('./src/'))

resx = 16;
resy = 16;
resz = 40;
% to reduce the computation time use the greatest devisor of the original resolution in level-2 segmentation
resx_s = 2; 
resy_s = 2;
resz_s = 5;

rootFolder = [data_path, '/segmentation_example_data/EM_segmentation_MICrONS/'];
input_folder = [data_path,'/segmentation_example_data/EM_segmentation_MICrONS/astrocyte_28_whole_dataset'];
skel_folder = [data_path,'/segmentation_example_data/EM_segmentation_MICrONS/864691135403747310_neuron_ds'];
% Here we divide the whole dataset into 5x5x5 chunks
curpsID = '864691135403747310';
neuro_ds_folder = [save_main_path,'/864691135403747310_neuron_ds'];
if(~exist(neuro_ds_folder, 'dir'))
    mkdir(neuro_ds_folder);
end
spine_save_folder = [save_main_path,'/864691135403747310_spine_shaft'];
if(~exist(spine_save_folder, 'dir'))
    mkdir(spine_save_folder);
end
branched_spine_save_folder = [save_main_path,'/864691135403747310_branchedspine'];
if(~exist(branched_spine_save_folder, 'dir'))
    mkdir(branched_spine_save_folder);
end
spine_head_neck_save_folder = [save_main_path,'/864691135403747310_spine_head_neck'];
if(~exist(spine_head_neck_save_folder, 'dir'))
    mkdir(spine_head_neck_save_folder);
end 
%%==========================================================================
% parameters:
cfg.ds_ratio = [5,5,5]; % the downsampled resolution will be 80x80x200 nm
cfg.ratio_radius = 3;
cfg.offset = 1200;
cfg.ds_ratio_computation = [2.5,2.5,5]; % the downsampled resolution for computation will be 32x32x40 nm which is closer to isotropic resolution
cfg.minSpineSize = 200;
cfg.second_path_length = 20;
cfg.curvature_scale = 30;
cfg.numP_level2 = 100;
cfg.smooth_level2 = 4;
%%==========================================================================

% 1. Preprocessing
ds_ratio = cfg.ds_ratio;
% Downsample the dendrite
[lenx, leny, lenz, neuron_x] = preprocessing.downsample_dendrite(input_folder, curpsID, ds_ratio, neuro_ds_folder);
% Generate skeleton
% run using Kimimaro python package at a lower resolution
% saved in the neuron_ds folder as '*_skel.tif'
% obtain skel_x
copyfile(fullfile(skel_folder, [curpsID,'_skel.tif']), fullfile(neuro_ds_folder, [curpsID,'_skel.tif']));
skel_x = tiffreadVolume(fullfile(neuro_ds_folder, [curpsID,'_skel.tif']));
% Remove soma
radius_threshold = 1000/ds_ratio(1); % remove soma from the downsampled neuron 
ratio_radius = cfg.ratio_radius; % the 
offset = cfg.offset;
mask_dendrite = preprocessing.remove_soma(neuro_ds_folder, curpsID, neuron_x, skel_x, resx, resy,resz,radius_threshold, ratio_radius, offset);

%%==========================================================================

% 2. Level-1 compartment segmentation
% mask_dendrite_dist = comSeg.cal_edt_dendrite(mask_dendrite, resx, resy, resz);
ds_ratio_computation = cfg.ds_ratio_computation;
minSpineSize = cfg.minSpineSize;
comSeg.level_1_segmentation_extract_spine(input_folder,neuro_ds_folder,curpsID,resx, resy, resz,ds_ratio, ds_ratio_computation, minSpineSize);

%%==========================================================================

% 3. Level-2 compartment segmentation
% extract non-branch dendrite spines and save each individual spine into a separate file
% each file contains the segmentation of both the individual spine and also the dendrite shaft region attached to it
second_path_length = cfg.second_path_length; % this decides which are the branched spines
preprocessing.level_2_preprocessing_remove_branch(input_folder, spine_save_folder, branched_spine_save_folder, curpsID, lenx, leny, lenz, resx, resy, resz, second_path_length);
% run volume-surface optimized segmentation of dendrite spine head and neck
curvature_scale = cfg.curvature_scale;
run_in_parallel = 1;
minpoints = cfg.numP_level2;
smooth_kernel = cfg.smooth_level2;
comSeg.level_2_segmentation_main_speed_up(spine_save_folder, spine_head_neck_save_folder, [], curvature_scale, resx_s, resy_s, resz_s,minpoints,smooth_kernel);

%%==========================================================================

% 4. quantification of the morphological features of dendrite spine and dendrite branches
% split the dendrite into branches and quantify the morphological features of each branch
dendrite_branches_Quant_save_folder = fullfile(save_main_path, [curpsID, '_dendrite_branches_quantification']);
if(~exist(dendrite_branches_Quant_save_folder, 'dir'))
    mkdir(dendrite_branches_Quant_save_folder);
end
structQuant.quantify_dendrite_branch(curpsID, lenx, leny, lenz,resx, resy, resz,input_folder,spine_save_folder,neuro_ds_folder, dendrite_branches_Quant_save_folder);
% quantify the morphological features of dendrite spines as well as the wrapping ratio of astrocyte
dendrite_spine_Quant_save_folder = fullfile(save_main_path, [curpsID, '_dendrite_spine_quantification']);
if(~exist(dendrite_spine_Quant_save_folder, 'dir'))
    mkdir(dendrite_spine_Quant_save_folder);
end

structQuant.quantify_dendrite_spine_score(curpsID, lenx, leny, lenz, input_folder, spine_save_folder,spine_head_neck_save_folder, dendrite_spine_Quant_save_folder, resx, resy, resz);

% classify the dendrite spines into different categories using the features obtained from the spine quantification
spine_classification_folder = fullfile(save_main_path, [curpsID, '_spine_classification']);
if(~exist(spine_classification_folder, 'dir'))
    mkdir(spine_classification_folder);
end
outlabel = structQuant.classification_dendrite_spine(spine_save_folder,dendrite_spine_Quant_save_folder, spine_classification_folder, resx, resy, resz);


%%==========================================================================

% 5. Summarize the results (plots...)

% % example plot the quantifications 
% multilayer_summary_path = ''
% plot_save_folder = fullfile(rootFolder, [curpsID, '_quantification_plots']);
% if(~exist(plot_save_folder, 'dir'))
%     mkdir(plot_save_folder);
% end
% % read in the summary from the quantification for each layer
% % both the quantification of the dendrite branches as well the quantification for the dendrite spines.

% [layer23_dendrite_summary, layer23_spine_summary, layer4_dendrite_summary, layer4_spine_summary, layer5_dendrite_summary, layer5_spine_summary] = structQuant.read_in_summary(multilayer_summary_path);
% structQuant.plot_summary_plot_example(layer23_dendrite_summary, layer23_spine_summary, layer4_dendrite_summary, layer4_spine_summary, layer5_dendrite_summary, layer5_spine_summary);

end