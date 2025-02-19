function test_level_1_segmentation(data_path, save_path)

    % contains the main code for running VSOT for the level-1 segmentation experiment

    % experiment 1
    % The data is from the dendrite #5 which was used in Kasthuri 2015 paper.
    % The workflow:
    % 1. since there is no soma in this test data, extract the apical dendrite
    % and its branches. Basde on the names of the structures in the annotation
    % file, as well as the hierachical relationshiop between different objects,
    % obtain each individual branches.
    % 2. Run our level-1 segmentation method
    % 3. Generate surface reconstruction for the volume which might be used in
    % other methods

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


    %% step 1 preprocessing the ground truth dataset
    % obtain the segmentation volume as well as the skeleton
    %
    % Run kimimaro for the skeleton generation 
    %% step 2 run our method
    list_of_folders = ["D5_Apical_Spines","D5_Branch_1","D5_Branch_2","D5_Branch_3","D5_Branch_4","D5_Branch_5"];

    resx = 24;
    resy = 24;
    resz = 30;

    rootfolder = [data_path,'/performance_testing_level_1_segmentation/dataset'];
    outputfolder = [data_path,'/performance_testing_level_1_segmentation/results/VSOT_result'];

    if(~exist(outputfolder, 'dir'))
        mkdir(outputfolder);
    end

    for i = 1:length(list_of_folders)
        cur_folder = convertStringsToChars(list_of_folders(i));
        mask_dendrite = tiffreadVolume(fullfile(rootfolder, cur_folder, [cur_folder,'_dendrite_volume.tif.tif'])) > 0;
        mask_skel = tiffreadVolume(fullfile(rootfolder, cur_folder, 'ds_skeleton.tif')) > 0;
        dist_dendrite = cal_edt_dendrite(mask_dendrite, resx, resy, resz);
        [lenx, leny ,lenz] = size(mask_dendrite);
        % when the size of the dendrite shaft is too large, like the length
        % is larger than 500, then split the volume into several chunks and
        % do the segmentation in each separate chunk
        % running directly on a very large volume will likely cost many
        % memory and at the same time it will make the surface biased to the
        % global constraint. In stead of local smooth surface.
        if(lenx > 600 || leny > 600 || lenz > 600)
            
            maskRemoveAll = false(lenx, leny ,lenz);
            lmx = min(lenx, 400);
            lmy = min(leny, 400);
            lmz = min(lenz, 400);
            num_x = floor(lenx/ lmx);
            num_y = floor(leny/ lmy);
            num_z = floor(lenz/ lmz);
            for ix = 0:num_x
                for iy = 0:num_y
                    for iz = 0:num_z
                        disp([ix, iy, iz])
                        winx = [max((1 + ix*lmx - round(lmx/2)), 1), min((ix + 1)*lmx + round(lmx/2), lenx);max((1 + iy*lmy - round(lmy/2)), 1), min((iy + 1)*lmy + round(lmy/2), leny);max((1 + iz*lmz - round(lmz/2)), 1), min((iz + 1)*lmz + round(lmz/2), lenz)];
                        tmp_den = mask_dendrite(winx(1,1):winx(1,2), winx(2,1):winx(2,2), winx(3,1):winx(3,2));
                        tmp_skel = mask_skel(winx(1,1):winx(1,2), winx(2,1):winx(2,2), winx(3,1):winx(3,2));
                        if(sum(tmp_skel(:)) > 100) 
                            tmp_dist = dist_dendrite(winx(1,1):winx(1,2), winx(2,1):winx(2,2), winx(3,1):winx(3,2));
                            % run segmentation based on the surface_volume optimization
                            maskRemove_tmp = comSeg.level_1_segmentation_function_fin(tmp_den, tmp_skel, tmp_dist,resx, resy, resz);
                            maskRemoveAll(winx(1,1):winx(1,2), winx(2,1):winx(2,2), winx(3,1):winx(3,2)) = maskRemove_tmp | maskRemoveAll(winx(1,1):winx(1,2), winx(2,1):winx(2,2), winx(3,1):winx(3,2));
                        end
                    end
                end
            end
        else
            maskRemoveAll = comSeg.level_1_segmentation_function_fin(mask_dendrite, mask_skel, dist_dendrite,resx, resy, resz);
        end
    %     out_dendrite = level_1_segmentation_function_fin(mask_dendrite,mask_skel,dist_dendrite, resx, resy, resz);
        mask_spine = xor(mask_dendrite,maskRemoveAll);
        if(~exist(fullfile(outputfolder, cur_folder), 'dir'))
            mkdir(fullfile(outputfolder, cur_folder));
        end
        tifwrite(uint8(255*mask_spine), fullfile(outputfolder, cur_folder, 'spine_segmentation'));
        tifwrite(uint8(255*maskRemoveAll), fullfile(outputfolder, cur_folder, 'shaft_segmentation')); 
    end
    mask_spine = [];
    maskRemoveAll = [];
    dist_dendrite = [];
    mask_dendrite = [];
    mask_skel = [];
    gen_final_results_all_methods(data_path);
    save_performance_testing_results(data_path, save_path);
end