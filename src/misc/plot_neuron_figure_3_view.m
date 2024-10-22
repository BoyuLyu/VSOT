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


parfor m = 1:size(group_info, 1)
    disp(m)
    tic;
    curpsID = group_info(m,2);
    lenx = group_info(m, 3);
    leny = group_info(m, 4);
    lenz = group_info(m, 5);
    mask_spine = false(5*lenx, 5*leny,5*lenz);
    mask_dendrite = false(5*lenx, 5*leny, 5*lenz);
    kk = group_info(m,1);
    disp(curpsID)

    fullSegfolder_root = [rootFolder,'astro_', num2str(kk), '_minnie65'];
    branched_spine_save_folder = fullfile(fullSegfolder_root, [num2str(curpsID), '_branchedspine']);
    branched_spine_file = dir(fullfile(branched_spine_save_folder,'*.tif'));

    for i = 1:length(branched_spine_file)
        try
            tmp_mask = tiffreadVolume(fullfile(branched_spine_save_folder, branched_spine_file(i).name));
            fig = figure('Visible', 'off');
            
            % Create the first subplot
            subplot(1, 3, 1);
            imagesc(max(tmp_mask, [],3));
            title('front view  (Z)');
            subplot(1, 3, 2);
            % Create the second subplot
            imagesc(squeeze(max(tmp_mask, [],2)));
            title('side view (Y)');
             subplot(1, 3, 3);           
            % Create the third subplot
            imagesc(squeeze(max(tmp_mask, [],1)));
            title('top view (X)');
            
            % Save the figure to a file
            saveas(fig, fullfile(branched_spine_save_folder, [branched_spine_file(i).name, '.jpg']));
            
            % Close the figure
            close(fig);
        catch ME
             continue;
        end
    end
    toc;
end


% 
% 
% 
% for m = 1:size(group_info,1)
%     try 
%     i = group_info(m,1);
%     disp(m)
%     fullSegfolder_root = [rootFolder,'astro_', num2str(i), '_minnie65'];
%        curpsID = group_info(m,2);
%        lenx = group_info(m, 3);
%        leny = group_info(m, 4);
%        lenz = group_info(m, 5);
%        combined_region = false(2.5*lenx, 2.5*leny, 2.5*lenz);
%        for ix = 0:4
%             for iy = 0:4
%                 for iz = 0:4
%                         fullSegfolder = fullfile(fullSegfolder_root, num2str(ix), num2str(iy),num2str(iz));
%                         if(exist(fullfile(fullSegfolder,['dendrite_',num2str(curpsID),'.tif']),'file') )
%                             mask_dendrite = logical(tiffreadVolume(fullfile(fullSegfolder,['dendrite_',num2str(curpsID),'.tif'])));
%                             mask_dendrite_ds = imresize3(mask_dendrite, [lenx/2, leny/2, lenz/2],'nearest');
%                             combined_region((1 + ix*0.5*lenx):(ix+1)*0.5*lenx, (1 + iy*0.5*leny):(iy+1)*0.5*leny, (1 + iz*0.5*lenz):(iz+1)*0.5*lenz) = mask_dendrite_ds;
%                         end
%                 end
%             end
%        end
%        close all
%         combined_region = imresize3(combined_region, [2*lenx, 2*leny, 2*lenz]);
%         id0 = find(combined_region(:) == 1);
%         [id0x, id0y, id0z] = ind2sub([2*lenx, 2*leny, 2*lenz], id0);
%         figure1 = volshow(combined_region(min(id0x): max(id0x), min(id0y): max(id0y), min(id0z): max(id0z)),'BackgroundColor','w');
%         h1 = figure1.Parent;
%         hFig = h1.Parent;
%         drawnow
%         I = getframe(hFig);
%         imwrite(I.cdata,fullfile('/work/boyu/EM_astrocyte/astro_11_28_16_16_40/dendrite_shape_figure', [num2str(i), '_', num2str(curpsID),'.jpg']))
% 
%     catch ME
%         warning('Problem occured at %s',[num2str(m)]);
%         continue;
%     end
% end
% 
% 
% 
% %         for ix = 0:4
% %             for iy = 0:4
% %                 for iz = 0:4
% %                     try
% %                     fullSegfolder = fullfile(fullSegfolder_root, num2str(ix), num2str(iy),num2str(iz));
% % 
% % 
% %                         if(exist(fullfile(fullSegfolder,['dendrite_',num2str(curpsID),'.tif']),'file') && ~exist(fullfile(fullSegfolder,['spine_new_no_soma_',num2str(curpsID),'.tif']), 'file'))
% %                             tic;
% %                             mask_dendrite = tiffreadVolume(fullfile(fullSegfolder,['dendrite_',num2str(curpsID),'.tif']));
% %                             mask_spine = level_1_segmentation_function_v3(mask_dendrite, xxshift, yyshift);
% %                             tifwrite(uint8(double(mask_spine > 0) + double(mask_dendrite > 0)),fullfile(fullSegfolder,['spine_new_no_soma_',num2str(curpsID)]))
% %                             toc;
% %     
% %                         end
% %                     catch ME
% %                         warning('Problem occured at %s',[num2str(ix),'/', num2str(iy),'/',num2str(iz), ' ', num2str(curpsID)]);
% %                         continue;
% %                     end
% % 
% %                 end
% %             end
% %         end
% % 
% %     end
% % 
% % end
% % 
% 
% 
% 
% 
% 
