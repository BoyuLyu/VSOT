function mask_spine_cell = level_1_segmentation_function_v4(mask_dendrite_cell, xxshift, yyshift,resx, resy, resz, minSpineSize)
%32xresyxresz
mask_spine_cell = cell(length(mask_dendrite_cell),1);
disp('start extracting shafts')

for i = 1:length(mask_dendrite_cell)
%     disp(i)
    mask_dendrite = mask_dendrite_cell{i} > 0;
    mask_skel = mask_dendrite_cell{i} == 2;
    [lenx, leny, lenz] = size(mask_dendrite);
    if(lenx > 20 && leny > 20 && lenz > 20)
        % certain chunks contain only a small fratgment of dendrite remove it
        if(lenx > 500 && leny > 500 && lenz > 500)
            % edt_mex can get stuck when the input region is really large
            mask_dendrite_dist = zeros(lenx, leny ,lenz);
            lmx = ceil(lenx/5);
            lmy = ceil(leny/5);
            lmz = ceil(lenz/5);
            for ix = 0:4
                for iy = 0:4
                    for iz = 0:4
    %                      disp([ix, iy, iz])
                        winx = [max((1 + ix*lmx - round(lmx/2)), 1), min((ix + 1)*lmx + round(lmx/2), lenx);max((1 + iy*lmy - round(lmy/2)), 1), min((iy + 1)*lmy + round(lmy/2), leny);max((1 + iz*lmz - round(lmz/2)), 1), min((iz + 1)*lmz + round(lmz/2), lenz)];
                        tmp_den = mask_dendrite(winx(1,1):winx(1,2), winx(2,1):winx(2,2), winx(3,1):winx(3,2));
                        tmp_den_1d = tmp_den(:);
                        if(sum(tmp_den_1d) > 10000)
                            disp([ix, iy, iz])
        %                     tic;
                            tmp_dist_1d = edt_mex(tmp_den_1d, size(tmp_den,1), size(tmp_den,2), size(tmp_den,3), resx,resy,resz);
                            tmp_dist = reshape(tmp_dist_1d, size(tmp_den));
                            mask_dendrite_dist(winx(1,1):winx(1,2), winx(2,1):winx(2,2), winx(3,1):winx(3,2)) = max(mask_dendrite_dist(winx(1,1):winx(1,2), winx(2,1):winx(2,2), winx(3,1):winx(3,2)), tmp_dist);
                        end
                    end
                end
            end
        else
            mask_dendrite(:,:,lenz) = 0;
            mask_dendrite(:,:,1) = 0;
            mask_dendrite(:,1,:) = 0;
            mask_dendrite(:,leny,:) = 0;
            mask_dendrite(1,:,:) = 0;
            mask_dendrite(lenx,:,:) = 0;
            mask_dendrite_1d = mask_dendrite(:);
            mask_dendrite_dist_1d = edt_mex(mask_dendrite_1d, size(mask_dendrite,1), size(mask_dendrite,2), size(mask_dendrite,3), resx,resy,resz);
            mask_dendrite_dist = reshape(mask_dendrite_dist_1d, size(mask_dendrite));
        
        end
    disp('start extracting shaft for first cluster')
        maskRemoveAll = comSeg.extract_shaft_v5(mask_dendrite, mask_skel, mask_dendrite_dist,xxshift,yyshift,resx, resy, resz);
        mask_spine_tmp = xor(mask_dendrite,maskRemoveAll);
        mask_spine_tmp_roi = bwlabeln(mask_spine_tmp);
        mask_spine_tmp_roi_idx = label2idx(mask_spine_tmp_roi);
        len_spine = cellfun(@length, mask_spine_tmp_roi_idx);
        mask_spine_tmp_roi_idx(len_spine < minSpineSize) = []; % remove really small noisy structures, not really necessary
        mask_spine_tmp_roi_idx = mask_spine_tmp_roi_idx(:);
        if(~isempty(mask_spine_tmp_roi_idx))
            mask_spine_tmp = false(lenx, leny, lenz);
            mask_spine_tmp(cell2mat(mask_spine_tmp_roi_idx)) = 1;
            mask_spine_cell{i} = mask_spine_tmp;
        end
    end




end
end