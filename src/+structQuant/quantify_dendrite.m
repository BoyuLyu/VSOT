clear

rootFolder = '/work/boyu/EM_astrocyte/astro_11_33_16_16_40/';
neuronListFolder = '/work/boyu/EM_astrocyte/astro_11_33_64_64_80/';

addpath('../resources/curvatures/')
addpath('../resources/ImprovedSurfaceSmooth/')
addpath('../resources/CurvatureEstimation')
addpath('/home/boyu/Documents/edt_mex/edt_mex/edt_mex')
addpath('../resources/TAUBIN/TAUBIN/')
addpath('/home/boyu/Documents/iso2mesh/')
addpath('../resources/data_to_off/')
addpath('/home/boyu/Documents/src_mex/mex_EM_analysis/mex_EM_analysis')
xxshift = zeros(3,3);
yyshift = zeros(3,3);
for i = -1:1
    for j = -1:1
        xxshift((i+2), (j+2)) = i;
        yyshift((i+2), (j+2)) = j;
    end
end
group_info = []; % contains the [id of astrocyte, neuronList_str, lenx, leny, lenz]
for i = 11:33
            
%     disp(i)
    fullSegfolder_root = [rootFolder,'astro_', num2str(i), '_minnie65'];
    neuronList = [neuronListFolder, 'astro_', num2str(i), '_minnie65/', 'top20_neuron_id_no_soma_filtered.txt'];
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
    tmpImg = tiffreadVolume(fullfile(fullSegfolder_root, '1','1', '1', 'astro_segMask.tif'));
    [lenx, leny, lenz] = size(tmpImg);
    tmpImg = [];
    group_info_j = [repmat(i,length(neuronList_str),1), neuronList_str(:), repmat([lenx, leny, lenz],length(neuronList_str),1)];
    group_info = [group_info; group_info_j];            


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
% 
for m = 1:size(group_info,1)
% for m = [16,120,180,240,300,330]
    tic;
    try
        curpsID = group_info(m,2);
        kk = group_info(m,1);
        fullSegfolder_root = [rootFolder,'astro_', num2str(kk), '_minnie65'];
        disp(curpsID)
        lenx = group_info(m, 3);
        leny = group_info(m, 4);
        lenz = group_info(m, 5);
        % dendrite_branches_save_folder = fullfile(fullSegfolder_root, [num2str(curpsID),'_dendrite_branches']);
        % if(~exist(dendrite_branches_save_folder, 'dir'))
        %     dendrite_feature  = quantify_dendrite_sub_func(curpsID, lenx, leny, lenz, group_info, rootFolder,fullSegfolder_root);
        %     if(~exist(dendrite_branches_save_folder, 'dir'))
        %         mkdir(dendrite_branches_save_folder);
        %     end
        %     writematrix(dendrite_feature, fullfile(dendrite_branches_save_folder,'dendrite_feature_w_mapping.csv'));
        % end
        outlabel = classification_dendrite_spine(fullSegfolder_root, curpsID);
        % export_dendrite_image_w_spinelabel(outlabel, curpsID, lenx, leny, lenz,fullSegfolder_root);
    catch ME
        continue;
    end


    toc;
end