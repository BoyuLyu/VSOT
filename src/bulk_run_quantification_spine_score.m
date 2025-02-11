
clear

rootFolder = '/work/boyu/EM_astrocyte/astro_11_33_16_16_40/';
neuronListFolder = '/work/boyu/EM_astrocyte/astro_11_33_64_64_80/';
xxshift = zeros(3,3);
yyshift = zeros(3,3);
for i = -1:1
    for j = -1:1
        xxshift((i+2), (j+2)) = i;
        yyshift((i+2), (j+2)) = j;
    end
end
group_info = []; % contains the [id of astrocyte, neuronList_str, lenx, leny, lenz]
total_count_spine = zeros(23,1);
for i = 11:33
            
%     disp(i)
    fullSegfolder_root = [rootFolder,'astro_', num2str(i), '_minnie65'];
    neuronList = [neuronListFolder, 'astro_', num2str(i), '_minnie65/', 'top20_neuron_id_no_soma_filtered.txt'];
    if(~exist(neuronList, "file"))
        neuronList = [neuronListFolder, 'astro_', num2str(nn), '_minnie65/', 'top20_neuronID_spine_centered.txt'];
        if(~exist(neuronList, "file"))
            warning('No neuronList file found')
            continue;
        end
    end
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
    for k = 1:length(neuronList_str)
        curpsID = neuronList_str(k);
        % spine_save_folder = fullfile(fullSegfolder_root, num2str(curpsID));
        spine_save_folder = fullfile(fullSegfolder_root, [num2str(curpsID)]);
        files_x = dir(fullfile(spine_save_folder, '*.tif'));
    end
end



parfor m = 1:size(group_info,1)
    try
    tic;
    curpsID = num2str(group_info(m,2));
    lenx = double(group_info(m,3));
    leny = double(group_info(m,4));
    lenz = double(group_info(m,5));
    kk = group_info(m,1);
    fullSegfolder_root = [rootFolder,'astro_', num2str(kk), '_minnie65'];
    matSaveFolder = fullfile(fullSegfolder_root, [num2str(curpsID),'_quantify_score_spine_only']);
    spine_save_folder = fullfile(fullSegfolder_root, num2str(curpsID));
    spine_head_neck_save_folder = fullfile(fullSegfolder_root, [num2str(curpsID),'_spine_head_neck']);
    dendrite_spine_Quant_save_folder = fullfile(fullSegfolder_root, [num2str(curpsID),'_quantify_score_v4']);
    resx = 16;
    resy = 16;
    resz = 40;
    input_folder = fullSegfolder_root;
    structQuant.quantify_dendrite_spine_score_single_psd(curpsID, lenx, leny, lenz, input_folder, spine_save_folder,spine_head_neck_save_folder, dendrite_spine_Quant_save_folder, resx, resy, resz);
    toc;
    catch ME
        warning('error happend at %s', num2str(m));
        continue;
    end
end