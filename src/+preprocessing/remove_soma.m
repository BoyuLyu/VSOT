function neuron_x = remove_soma(folder, filename, neuron_x, skel_x, resx, resy,resz, radius_threshold, ratio_radius, offset)

 % input: neuron segmentation and skeleton
 % output: dendrite segmentation without soma
 % save the results to the folder
 % ratio_radius: the ratio of the radius of the soma to the radius of the dendrite, defaulth 3
 % offset: the offset of the radius of the soma to the radius of the dendrite, default 1200
    neuron_x = neuron_x > 0;
    neuron_x_before = neuron_x;
    skel_x = skel_x > 0;
    % se2 = strel('sphere', 10);
    % skel_x2 = imdilate(skel_x, se2);
    [lenx, leny, lenz] = size(neuron_x);
    mask_dendrite_1D = neuron_x(:);
    dist_dendrite_1D = edt_mex(mask_dendrite_1D, lenx, leny, lenz, resx,resy,resz);
    mask_dendrite_1D = [];
    dist_dendrite = reshape(dist_dendrite_1D, lenx, leny, lenz);
    dist_skel = dist_dendrite.*double(skel_x);
    while max(dist_skel(:)) > radius_threshold
        neuron_x_roi = bwlabeln(neuron_x);
        neuron_x_roi_idx = label2idx(neuron_x_roi);
        for i = 1:length(neuron_x_roi_idx)
            tmpid = neuron_x_roi_idx{i};
            tmpid_skel = tmpid(skel_x(tmpid)>0);
            if(max(dist_dendrite(tmpid_skel)) == max(dist_skel(:)))
                centerID = tmpid(find(dist_dendrite(tmpid) == max(dist_skel(:)), 1));
                mask0 = false(lenx, leny, lenz);
                mask0(centerID) = 1;
                inputx = false(lenx, leny, lenz);
                inputx(tmpid) = 1;
                mask01D = mask0(:);
                inputx1D = inputx(:);
                dist0 = imchaferDist3D(inputx1D,mask01D,lenx, leny, lenz, resx, resy, resz);
                dist0(dist0 == 1e10) = 0;
                dist_2_center = reshape(dist0, lenx, leny, lenz);
                outx = dist_2_center > (ratio_radius*max(dist_skel(:)) + offset); % 3 & 500 are the coefficients like teaser algorithm
                removed_ID = tmpid(outx(tmpid) == 0);
                dist_skel(removed_ID) = 0;
                skel_x(removed_ID) = 0;
                neuron_x(removed_ID) = 0;
            end
        end
    end
    
    dist_dendrite = [];
    dist_dendrite_1D = [];
    tifwrite(double(neuron_x), fullfile(folder, [filename,'_dendrite_without_soma.tif']));
    removed_x = neuron_x_before - neuron_x;
    se2 = strel('disk',10);
    removed_x = imdilate(removed_x , se2).*neuron_x_before;
    tifwrite(double(removed_x), fullfile(folder, [filename,'_dendrite_soma']));
end