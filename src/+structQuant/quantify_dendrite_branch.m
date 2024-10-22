function quantify_dendrite_branch(curpsID, lenx, leny, lenz,resx, resy, resz,fullSegfolder_root, spine_save_folder,neuron_ds_folder, dendrite_branches_save_folder)
    dendrite_feature  = structQuant.quantify_dendrite_sub_func(curpsID, lenx, leny, lenz, resx, resy, resz, fullSegfolder_root,spine_save_folder, neuron_ds_folder);
    if(~exist(dendrite_branches_save_folder, 'dir'))
        mkdir(dendrite_branches_save_folder);
    end
    column_names = {'Spine ID', 'Dendrite ID', 'Length', 'Radius'}; % Replace with actual column names
    dendrite_feature_table = array2table(dendrite_feature, 'VariableNames', column_names);
    writetable(dendrite_feature_table, fullfile(dendrite_branches_save_folder, 'dendrite_feature_w_mapping.csv'));
    % writematrix(dendrite_feature, fullfile(dendrite_branches_save_folder,'dendrite_feature_w_mapping.csv'));
end