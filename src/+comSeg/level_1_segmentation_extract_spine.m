function level_1_segmentation_extract_spine(fullSegfolder_root, neuro_ds_folder, curpsID, resx, resy, resz, ds_ratio, ds_ratio_computation, minSpineSize)

xxshift = zeros(3,3);
yyshift = zeros(3,3);
for i = -1:1
    for j = -1:1
        xxshift((i+2), (j+2)) = i;
        yyshift((i+2), (j+2)) = j;
    end
end

xxshift3D = zeros(3,3,3);
yyshift3D = zeros(3,3,3);
zzshift3D = zeros(3,3,3);
for i = -1:1
    for j =-1:1
        for k = -1:1
            xxshift3D((i+2), (j+2), (k+2)) = i;
            yyshift3D((i+2), (j+2), (k+2)) = j;
            zzshift3D((i+2), (j+2), (k+2)) = k;
        end
    end
end

tmpImg = tiffreadVolume(fullfile(fullSegfolder_root, '0','0', '0', 'new_astrocyte_seg.tif'));
[lenx, leny, lenz] = size(tmpImg);
tmpImg = [];

neuron_x_removal = tiffreadVolume(fullfile(neuro_ds_folder, [num2str(curpsID),'_dendrite_soma.tif'])) > 0;
neuron_x_removal_roi = bwlabeln(neuron_x_removal);
neuron_x_removal_roi_idx = label2idx(neuron_x_removal_roi);
neuron_x_removal_roi_idx(cellfun(@length, neuron_x_removal_roi_idx) < 1000) = [];
neuron_x_removal_roi_idx = neuron_x_removal_roi_idx(:);
neuron_x_removal = false(size(neuron_x_removal));
neuron_x_removal(cell2mat(neuron_x_removal_roi_idx)) = 1;
skel_x = tiffreadVolume(fullfile(neuro_ds_folder, [curpsID,'_skel.tif']));
skel_x = skel_x > 0;
se = strel('sphere', 1);
neuron_x_removal = imdilate(neuron_x_removal, se);
skel_x = skel_x.*(1 - neuron_x_removal) > 0;

combined_region = false(ds_ratio(1)*lenx, ds_ratio(2)*leny, ds_ratio(3)*lenz);
neuron_x_removal_ori_scale = imresize3(neuron_x_removal,[ds_ratio(1)*lenx, ds_ratio(2)*leny, ds_ratio(3)*lenz], "Method","nearest");
neuron_x_removal = [];
try
    for ix = 0:4
        for iy = 0:4
            for iz = 0:4
                fullSegfolder = fullfile(fullSegfolder_root, num2str(ix), num2str(iy),num2str(iz));
                if(exist(fullfile(fullSegfolder,['dendrite_',curpsID,'.tif']),'file') )
                    mask_dendrite = tiffreadVolume(fullfile(fullSegfolder,['dendrite_',curpsID,'.tif'])) > 0;
                    combined_region((1 + ix*lenx):(ix+1)*lenx, (1 + iy*leny):(iy+1)*leny, (1 + iz*lenz):(iz+1)*lenz) = mask_dendrite;
                end
            end
        end
    end
    combined_region2 = combined_region.*(1 - neuron_x_removal_ori_scale) > 0;

    combined_region2 = imresize3(combined_region2, [ds_ratio_computation(1)*lenx, ds_ratio_computation(2)*leny, ds_ratio_computation(3)*lenz], "method","nearest");% 32nmx32nmx40nm
    neuron_x_removal_ori_scale = [];
    clear neuron_x_removal neuron_x_removal_ori_scale

    % To avoid the redundant memory usage, we use bounding box to extract each dendrite branches by splitting the skeleton using the graph split

    [out_skel_small_scale, dendrite_cell,out_bounding_box] = comSeg.get_bbx_from_split_skeleton(skel_x,combined_region2,xxshift3D, yyshift3D, zzshift3D);

    resx_new = resx/ds_ratio_computation(1)*ds_ratio(1);
    resy_new = resy/ds_ratio_computation(2)*ds_ratio(2);
    resz_new = resz/ds_ratio_computation(3)*ds_ratio(3);
    tic;
    mask_spine_cell = comSeg.level_1_segmentation_function_v4(dendrite_cell, xxshift, yyshift,resx_new, resy_new, resz_new, minSpineSize);
    toc;
        
    output_spine = false(size(combined_region2));
    for i = 1:length(mask_spine_cell)
        if(~isempty(mask_spine_cell{i}))
        output_spine(out_bounding_box(i,1):out_bounding_box(i,2), out_bounding_box(i,3):out_bounding_box(i,4), out_bounding_box(i,5):out_bounding_box(i,6)) =...
            output_spine(out_bounding_box(i,1):out_bounding_box(i,2), out_bounding_box(i,3):out_bounding_box(i,4), out_bounding_box(i,5):out_bounding_box(i,6))...
            + mask_spine_cell{i};
        end
    end
    mask_spine_cell = [];
    clear mask_spine_cell
    mask_spine = imresize3(output_spine, [ds_ratio(1)*lenx, ds_ratio(2)*leny, ds_ratio(3)*lenz],  "method","nearest");
    % tifwrite(uint8(output_spine + imresize3(combined_region, [2.5*lenx, 2.5*leny, 5*lenz], "method","nearest")),['/work/boyu/EM_astrocyte/test_level_1_segmentation/tripartite_ROI/',num2str(curpsID),'_32_32_40'])                                         
    se = strel('sphere',1);
    for ix = 0:4
        for iy = 0:4
            for iz = 0:4
                fullSegfolder = fullfile(fullSegfolder_root, num2str(ix), num2str(iy),num2str(iz));
                if(exist(fullfile(fullSegfolder,['dendrite_',num2str(curpsID),'.tif']),'file') )
%                             mask_dendrite = logical(tiffreadVolume(fullfile(fullSegfolder,['dendrite_',num2str(curpsID),'.tif'])));
                    mask_spine_ds = mask_spine((1 + ix*lenx):(ix+1)*lenx, (1 + iy*leny):(iy+1)*leny, (1 + iz*lenz):(iz+1)*lenz);
                    combined_region_local = combined_region((1 + ix*lenx):(ix+1)*lenx, (1 + iy*leny):(iy+1)*leny, (1 + iz*lenz):(iz+1)*lenz);
                    tmp = double(mask_spine_ds) + double(combined_region_local);
                    tifwrite(uint8(tmp),fullfile(fullSegfolder,['spine_new_no_soma_revised_',curpsID]))                        
                end
            end
        end
    end 
    skel_x = [];
    mask_spine = [];
    output_spine = [];
    combined_region2 = [];
    combined_region = [];
    clear skel_x mask_spine output_spine combined_region2 combined_region
catch ME
    warning('Problem occured');
end

end







