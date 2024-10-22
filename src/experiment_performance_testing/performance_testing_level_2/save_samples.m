% Save the surface file for the dendrite spine volumes which will be used for testing the performance.'
% the resolution is changed to [2x2x5]

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
addpath('../resources/data_to_off/')
addpath('/home/boyu/Documents/graph_related/graph_mex/')
xxshift = zeros(3,3);
yyshift = zeros(3,3);
for i = -1:1
    for j = -1:1
        xxshift((i+2), (j+2)) = i;
        yyshift((i+2), (j+2)) = j;
    end
end
for nn = 23:28
        % nn=22; m =3;ix = 1;iy = 2;iz= 2 ; i = 3
    disp(nn)
    fullSegfolder_root = [rootFolder,'astro_', num2str(nn), '_minnie65'];
    neuronList = [neuronListFolder, 'astro_', num2str(nn), '_minnie65/', 'top20_neuron_id_no_soma.txt'];
    opts = delimitedTextImportOptions("NumVariables", 1);

    % Specify range and delimiter
    opts.DataLines = [1, Inf];
    opts.Delimiter = ",";
    
    % Specify column names and types
    opts.VariableNames = "VarName1";
    opts.VariableTypes = "uint64";
    
    % Specify file level properties
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";    
    neuronList_str = table2array(readtable(neuronList, opts));
    [xxshift, yyshift, zzshift] = genShiftMatrixTetra;
    for m = 1:length(neuronList_str)
%         disp([nn, m])
         for ix = 0:4
            for iy = 0:4
                for iz = 0:4
                    try
                        fullSegfolder = fullfile(fullSegfolder_root, num2str(ix), num2str(iy),num2str(iz));
                        curpsID = neuronList_str(m);
                        if(exist(fullfile(fullSegfolder,['spine_new_no_soma_',num2str(curpsID),'.tif']),'file') && exist(fullfile(fullSegfolder,['spine_branched_ROI_',num2str(curpsID),'.tif']),'file'))
    %                             disp([ix, iy, iz])
                                targetFolder = fullSegfolder;
                                spineROI = tiffreadVolume(fullfile(targetFolder,['spine_branched_ROI_',num2str(curpsID),'.tif']));
                                spineMask0 = tiffreadVolume(fullfile(targetFolder,['spine_new_no_soma_',num2str(curpsID),'.tif']));
                                
                                shaft_mask = spineMask0 == 1;
    %                             spineMask0 = [];
                                spineROIidx = label2idx(spineROI);
                                [lenx, leny, lenz] = size(spineROI);
                                spine_summary = table2array(readtable(fullfile(targetFolder,['spine_branched_',num2str(curpsID),'_branch_summary.csv'])));
                                if(~isempty(spine_summary))
                                    spine_summary = spine_summary(:,1);
                                    spineHead = double(shaft_mask);
                                    for i = 1:length(spineROIidx)
                                        if(spine_summary(i) == 0)
                                            %(i = 4, 11, 10, 12)
    
                                            [lenx, leny, lenz] = size(spineROI);
                                            test_region = spineROI == i;
                        %                     disp(i)
                                            curID0 = spineROIidx{i};
                                            [idtmpx, idtmpy, idtmpz] = ind2sub([lenx, leny, lenz], curID0);
                                            bbx = [max(min(idtmpx) -5, 1), min(max(idtmpx) + 5, lenx);max(min(idtmpy) - 5, 1), min(max(idtmpy)+5, leny);max(min(idtmpz) - 5, 1), min(max(idtmpz)+5, lenz)];
                                            test_region = test_region(bbx(1,1): bbx(1,2),bbx(2,1):bbx(2,2),bbx(3,1): bbx(3,2));
                                            tifwrite(double(test_region),fullfile(outputFolder,[num2str(nn), '_', num2str(ix), num2str(iy),num2str(iz), '_', num2str(m), '_', num2str(i)]));
    
                                            [node,elem,face,regions]=vol2surf(double(test_region),1:size(test_region,1),1:size(test_region,2),1:size(test_region,3),2,1,'cgalsurf');
                                            Tri = elem(:,1:3);
                                            Pts = node;
                                            Pts(:,1) = Pts(:,1) .*2;
                                            Pts(:,2) = Pts(:,2) .*2;
                                            Pts(:,3) = Pts(:,3) .*5;
                                            [Pts] = taubinsmooth( Tri,[Pts(:,1),Pts(:,2),Pts(:,3)],10);
    %                                         figure; trisurf(Tri, Pts(:,1), Pts(:,2), Pts(:,3));
                                            data_to_off(Tri',Pts',fullfile(outputFolder, [num2str(nn), '_', num2str(ix), num2str(iy),num2str(iz), '_', num2str(m), '_', num2str(i),'.off']));
    
    
                                        end
                                    end
                                end
                        end
                    catch ME
                        continue;
                    end
                end
            end
         end
    end
end

                                 