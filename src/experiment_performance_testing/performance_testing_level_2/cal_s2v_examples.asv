function cal_s2v_examples(offFolder, tifFolder, s2vSavePath)
    files = dir(fullfile(offFolder, '*.off'));
    filesname = {files.name};
    s2v_all = zeros(length(filesname), 1) + nan;
    for j = 1:length(filesname)
        id_3 = strsplit(filesname{j}, { '.'});
        
        name_file = id_3{1};
        if(exist(fullfile(tifFolder, [name_file, '.tif_head_neck.tif']), 'file'))
            labelx = tiffreadVolume(fullfile(tifFolder, [name_file, '.tif_head_neck.tif']));
            [lenx, leny, lenz] = size(labelx);
            spineHead = labelx == 3;
            spineHeadROI = bwlabeln(spineHead);
            spineHeadROIIdx = label2idx(spineHeadROI);
            spineHeadROIIdx = spineHeadROIIdx(:);
            len_region = cellfun(@length, spineHeadROIIdx);
            [~, idx] = max(len_region);
            test_region = spineHeadROI == idx;
            [node,elem,face,regions]=vol2surf(double(test_region),1:size(test_region,1),1:size(test_region,2),1:size(test_region,3),1.01,1,'cgalsurf');
            Tri = elem(:,1:3);
            Pts = node;
            Pts(:,1) = Pts(:,1) .*2;
            Pts(:,2) = Pts(:,2) .*2;
            Pts(:,3) = Pts(:,3) .*5;
            [Pts] = taubinsmooth( Tri,[Pts(:,1),Pts(:,2),Pts(:,3)],10);
            vectorMF3_1 = Pts(Tri(:,3),:) - Pts(Tri(:,1),:);
            vectorMF3_2 = Pts(Tri(:,3),:) - Pts(Tri(:,2),:);
            ss_tmp = cross(vectorMF3_1, vectorMF3_2,2);
            ss = 1/2*sum(sqrt(ss_tmp(:,1).^2 + ss_tmp(:,2).^2 + ss_tmp(:,3).^2));
            volume = abs(comSeg.calVolfromMesh(Tri, Pts));                                               
            r1 = (volume/(4*pi/3))^(1/3);
            r2 = sqrt(ss/(4*pi));
            s2v_all(j) = abs(r1 - r2)/r2;
        end
    end
    T = table(filesname(:), 1 - s2v_all, 'VariableNames', {'Filename', 'S2V'});
    writetable(T, fullfile(s2vSavePath, 's2v_result.csv'));
end
