function export_dendrite_image_w_spinelabel(typelabels, curpsID, lenx, leny, lenz,fullSegfolder_root)

    mask_spine = false(5*lenx, 5*leny,5*lenz);
    mask_dendrite = false(5*lenx, 5*leny, 5*lenz);
    spine_save_folder = fullfile(fullSegfolder_root, num2str(curpsID));
    for ix = 0:4
        for iy = 0:4
            for iz = 0:4
                fullSegfolder = fullfile(fullSegfolder_root, num2str(ix), num2str(iy),num2str(iz));
                if(exist(fullfile(fullSegfolder,['spine_new_no_soma_revised_',num2str(curpsID),'.tif']),'file') )
                    tmp = tiffreadVolume(fullfile(fullSegfolder,['spine_new_no_soma_revised_',num2str(curpsID),'.tif']));
                    mask_spine((1 + ix*lenx):(ix+1)*lenx, (1 + iy*leny):(iy+1)*leny, (1 + iz*lenz):(iz+1)*lenz) = tmp == 2;
                    mask_dendrite((1 + ix*lenx):(ix+1)*lenx, (1 + iy*leny):(iy+1)*leny, (1 + iz*lenz):(iz+1)*lenz) = tmp == 1;
    %                 tifwrite(uint8(tmp),fullfile(fullSegfolder,['spine_new_no_soma_revised_',num2str(curpsID)]))                        
                elseif(exist(fullfile(fullSegfolder,['spine_new_no_soma_',num2str(curpsID),'.tif']),'file'))
                    tmp = tiffreadVolume(fullfile(fullSegfolder,['spine_new_no_soma_',num2str(curpsID),'.tif']));
                    mask_spine((1 + ix*lenx):(ix+1)*lenx, (1 + iy*leny):(iy+1)*leny, (1 + iz*lenz):(iz+1)*lenz) = tmp == 2;
                    mask_dendrite((1 + ix*lenx):(ix+1)*lenx, (1 + iy*leny):(iy+1)*leny, (1 + iz*lenz):(iz+1)*lenz) = tmp == 1;
    %               
                end
            end
        end
    end 
    spine_coordinates = table2array(readtable(fullfile(spine_save_folder, 'spine_coordinate.csv')));
    % store the saved spine ID, saved dendrite ID, the length of the corresponding dendrite, the thickness of the corresponding dendrite
    se = strel("sphere", 1);
    mask_spine_erode = imerode(mask_spine, se);
    mask_spine_roi = bwlabeln(mask_spine);
    mask_spine_roi_idx = label2idx(mask_spine_roi);
    output_image_original_scale = uint8(mask_dendrite);
    for i = 1:size(spine_coordinates,1)
        spine_small = tiffreadVolume(fullfile(spine_save_folder, [num2str(i), '.tif']));
        [lenx_small, leny_small, lenz_small] = size(spine_small);
        spine_labels = mask_spine_roi(spine_coordinates(i,1):(spine_coordinates(i,1) + lenx_small - 1)...
            , spine_coordinates(i, 2): (spine_coordinates(i, 2) + leny_small - 1), ...
            spine_coordinates(i, 3):(spine_coordinates(i, 3)+lenz_small - 1));
        spine_labels(spine_labels(:) == 0) = [];
        spine_label = mode(spine_labels);
        if(spine_label ~= 0 && ~isnan(spine_label))
            output_image_original_scale(mask_spine_roi_idx{spine_label}) = typelabels(i) + 1;
        end
    end
    % output_image_ds = imresize3(output_image_original_scale, [lenx, leny, 5*lenz], 'Method', 'nearest');
    level_1_type_visualization_folder = '/work/boyu/EM_astrocyte/test_level_1_segmentation/tripartite_roi_with_type_label';
    tifwrite(uint8(output_image_original_scale), fullfile(level_1_type_visualization_folder, [num2str(curpsID),'dendrite_with_spine_label.tif']));
end