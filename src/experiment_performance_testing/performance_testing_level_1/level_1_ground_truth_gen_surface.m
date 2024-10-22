clear
rootFolder = '/work/boyu/EM_astrocyte/astro_11_28_16_16_40/';
neuronListFolder = '/work/boyu/EM_astrocyte/astro_11_28_64_64_80/';
outputFolder = '/work/boyu/EM_astrocyte/test_segmentation_samples/dendrite_spines';
addpath('../resources/curvatures/')
addpath('../resources/ImprovedSurfaceSmooth/')
addpath('../resources/CurvatureEstimation')
addpath('/home/boyu/Documents/edt_mex/edt_mex/edt_mex')
addpath('../resources/TAUBIN/TAUBIN/')
addpath('/home/boyu/Documents/tinymesh_mex/tiny_mesh_mex')
addpath('/home/boyu/Documents/iso2mesh/')
addpath('../resources/data_to_off/')
addpath('/home/boyu/Documents/graph_related/graph_mex/')
segment_name = 'D5_Apical_Spines'; %list_of_segments = ["D5_Apical_Spines","D5_Branch_1","D5_Branch_2","D5_Branch_3","D5_Branch_4","D5_Branch_5"];
rootfolder= ['../data/performance_testing_level_1_segmentation/', segment_name];
test_region = tiffreadVolume(fullfile(rootfolder, [segment_name,'_dendrite_volume.tif.tif'])) > 0;
%filter the test region so that there is only one segmentation
test_region_roi = bwlabeln(test_region);
test_region_roi_idx = label2idx(test_region_roi);
len_all= cellfun(@length, test_region_roi_idx);
tmp_id = test_region_roi_idx(len_all == max(len_all));
test_region = test_region.*0;
test_region(tmp_id{1}) = 1;
test_region2 = test_region == 1;
[node,elem,face,regions]=vol2surf((test_region2),1:size(test_region2,1),1:size(test_region2,2),1:size(test_region2,3),1.1,1,'cgalsurf');
Tri = elem(:,1:3);
Pts = node;
figure;trisurf(Tri,Pts(:,1),Pts(:,2),Pts(:,3),'Facecolor','red','FaceAlpha',0.1);
% write outputs as the format of the input for NeuRD, txt
writematrix(Tri, fullfile(rootfolder, 'tri.txt'));
writematrix(Pts, fullfile(rootfolder, 'vert.txt'));
% write .off file for later visualization
data_to_off(Tri',Pts',fullfile(rootfolder, 'dendrite.off'));


