function [lenx, leny, lenz, combined_region]  = downsample_dendrite(fullSegfolder_root,curpsID, ds_ratio, output_folder)
    % DOWNSAMPLE_DENDRITE
    % fullSegfolder_root: the root folder of the full segmentation. This folder should contain all the chunks inside the target region.
    % curpsID: the ID of the current post-synaptic neuron
    % ds_ratio: the downsample ratio of the full segmentation (for 3D volume this is a 3x1 vector)
    % output_folder: the folder to save the downsampled dendrite (neuron_ds folder)

% addpath('../resources/curvatures/')
% addpath('../resources/ImprovedSurfaceSmooth/')
% addpath('../resources/CurvatureEstimation')
% addpath('/home/boyu/Documents/edt_mex/edt_mex/edt_mex')
% addpath('../resources/TAUBIN/TAUBIN/')
% addpath('/home/boyu/Documents/iso2mesh/')
% addpath('../resources/data_to_off/')
% addpath('/home/boyu/Documents/src_mex/mex_EM_analysis/mex_EM_analysis')
xxshift = zeros(3,3);
yyshift = zeros(3,3);
for i = -1:1
    for j = -1:1
        xxshift((i+2), (j+2)) = i;
        yyshift((i+2), (j+2)) = j;
    end
end

% read in a random chunk to get the size of the full segmentation
tmpImg = tiffreadVolume(fullfile(fullSegfolder_root, '0','0', '0', 'new_astrocyte_seg.tif'));
[lenx, leny, lenz] = size(tmpImg);
tmpImg = [];
% we define the size of the downsampled for the whole region to be the size of each chunk
combined_region = false(lenx, leny, lenz); % normally set the size
for ix = 0:4
    for iy = 0:4
        for iz = 0:4
            fullSegfolder = fullfile(fullSegfolder_root, num2str(ix), num2str(iy),num2str(iz));
            if(exist(fullfile(fullSegfolder,['dendrite_',num2str(curpsID),'.tif']),'file') )
                mask_dendrite = logical(tiffreadVolume(fullfile(fullSegfolder,['dendrite_',num2str(curpsID),'.tif'])));
                mask_dendrite_ds = imresize3(mask_dendrite, [lenx/ds_ratio(1), leny/ds_ratio(2), lenz/ds_ratio(3)],'nearest');
                combined_region((1 + ix/ds_ratio(1)*lenx):(ix+1)/ds_ratio(1)*lenx, (1 + iy/ds_ratio(2)*leny):(iy+1)/ds_ratio(2)*leny, (1 + iz/ds_ratio(3)*lenz):(iz+1)/ds_ratio(3)*lenz) = mask_dendrite_ds;
            end
        end
    end
    
end
if(max(combined_region(:)) ==0)
    disp([fullSegfolder_root, num2str(curpsID)])
else
    tifwrite(uint8(double(combined_region)), fullfile(output_folder, [num2str(curpsID)]));
end



end