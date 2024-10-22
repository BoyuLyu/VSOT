% se3 = strel('sphere', 5);
% hvolume = volshow(imdilate(skel_x_parts_roi, se3), 'Colormap',cmp);
neuron_x_removal_big = imresize3(neuron_x_removal, [5*lenx, 5*leny, 5*lenz], 'Method', 'nearest');



final_mask_original_scale2 = final_mask_original_scale2;
final_mask_original_scale2_idx = label2idx(final_mask_original_scale2);
final_mask_original_scale2_idx = [final_mask_original_scale2_idx(:)];

final_mask_original_scale3 = zeros(size(final_mask_original_scale2));
final_mask_original_scale3_idx = final_mask_original_scale2_idx(randperm(length(final_mask_original_scale2_idx)));
for i = 1:length(final_mask_original_scale3_idx)
    final_mask_original_scale3(final_mask_original_scale3_idx{i}) = i;
end

final_mask_original_scale3 = final_mask_original_scale3 + (final_mask_original_scale3 > 0);
final_mask_original_scale3 = final_mask_original_scale3 + neuron_x_removal_big;
tifwrite(imresize3(uint8(final_mask_original_scale3), [lenx, leny, 2.5*lenz],'Method','nearest'), '/work/boyu/EM_astrocyte/materials_for_paper_VSOT/application_3/dendrite_segmentation/roi/dendrite_seg')

hvolume = volshow(imresize3(uint8(final_mask_original_scale3), [lenx, leny, lenz], 'Method','nearest') + uint8(neuron_x_removal),'Colormap',cmp);
viewer = hvolume.Parent;
viewer.CameraZoom = 1.5;
viewer.BackgroundColor="white";
viewer.BackgroundGradient="off";

hvolume = volshow(mask_dendrite_ds_all);
viewer = hvolume.Parent;
viewer.CameraZoom = 1.5;
viewer.BackgroundColor="white";
viewer.BackgroundGradient="off";

se2 = strel('sphere', 3);
hvolume = volshow(imdilate(skel_x_parts_roi >0 , se2));
viewer = hvolume.Parent;
viewer.CameraZoom = 1.5;
viewer.BackgroundColor="white";
viewer.BackgroundGradient="off";


cmp = jet;
% [~,cmp] = cmpermute([], jet);
se2 = strel('sphere', 3);
hvolume = volshow(imdilate(final_mask_ds_scale3.*uint8(skel_x_parts_roi >0) , se2), 'Colormap',cmp);
viewer = hvolume.Parent;
viewer.CameraZoom = 1.5;
viewer.BackgroundColor="white";
viewer.BackgroundGradient="off";


final_mask_ds_scale3 = imresize3(uint8(final_mask_original_scale3), [lenx, leny, lenz], 'Method','nearest');
hvolume = volshow( final_mask_ds_scale3,'Colormap',cmp);
viewer = hvolume.Parent;
viewer.CameraZoom = 1.5;
viewer.BackgroundColor="white";
viewer.BackgroundGradient="off";






cmp = jet;
% [~,cmp] = cmpermute([], jet);
se2 = strel('sphere', 3);
hvolume = volshow(mask_dendrite_ds_dist, 'Colormap',cmp);
viewer = hvolume.Parent;
viewer.CameraZoom = 1.5;
viewer.BackgroundColor="white";
viewer.BackgroundGradient="off";




