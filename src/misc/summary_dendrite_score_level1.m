function layer2_summary_level1_score = summary_dendrite_score_level1(nf_list,rootFolder,neuronListFolder)
    %% record the 

    spine_dendrite_type_label = []; % [1: filopodia, 2: mushroom, 3: long thin, 4: thin, 5: stubby, dendrite_label] Ns x 2
    dendrite_length  = []; % Nd x 1
    dendrite_radius = []; % Nd x 1 
    dendrite_each_spine_type_number = []; % specifically add the 6th type as the branched spine, simply to check between different layers Nd x 6
    neuron_total_type_number = []; % Nn x 6
    max_dendrite_label = 0;
    ttc = 0;
    for nf = nf_list
        disp(nf)
        fullSegfolder_root = [rootFolder,'astro_', num2str(nf), '_minnie65'];
        neuronList = [neuronListFolder, 'astro_', num2str(nf), '_minnie65/', 'top20_neuron_id_no_soma_filtered.txt'];
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
        for m = 1:length(neuronList_str)
            curpsID = neuronList_str(m);
            matSaveFolder = fullfile(fullSegfolder_root, [num2str(curpsID),'_quantify_score_v4']);
            if(exist(fullfile(matSaveFolder, 'summary.mat'), 'file'))
                load(fullfile(matSaveFolder, 'summary.mat')); %load as the struct with name summaryStructure
            else
                continue;
            end
            if(~isempty(summaryStructure.ss1))
                ttc  = ttc + length(summaryStructure.ss1);
                spineSaveFolder = fullfile(fullSegfolder_root, [num2str(curpsID)]);
                dendrite_score_folder = fullfile(fullSegfolder_root, [num2str(curpsID), '_dendrite_branches']);
                dendrite_score = readtable(fullfile(dendrite_score_folder, 'dendrite_feature_w_mapping.csv'));
                dendrite_score = table2array(dendrite_score);
                single_spine_coordinates = table2array(readtable(fullfile(spineSaveFolder, 'spine_coordinate.csv')));
    
                cur_max_dendrite_label = max(dendrite_score(:,2));
                dendrite_each_spine_type_number_tmp = zeros(cur_max_dendrite_label, 6); % stores the total number of each type of spine for each dendrite branch
                filoFolder = fullfile(spineSaveFolder, 'filopodia');
                mushFolder = fullfile(spineSaveFolder, 'mushroom');
                longthinFolder = fullfile(spineSaveFolder, 'longthin');
                thinFolder = fullfile(spineSaveFolder, 'thin');
                stubbyFolder = fullfile(spineSaveFolder, 'stubby');
                branchedFolder = fullfile(fullSegfolder_root, [num2str(curpsID),'_branchedspine']);
                branched_spine_coordinates = table2array(readtable(fullfile(branchedFolder, 'spine_coordinate.csv')));
                % other than the branched spine folder, we will stack all the labels for different types of spines from the rest of folders
                single_spine_length = length(summaryStructure.ss1);
                spine_dendrite_type_label_tmp = zeros(single_spine_length, 2); %store both the type of each dendrite spine as well as the label of the corresponding dendrite branch
                % go over all the spine type folders and record their types
                filo_count = 0;
                mush_count = 0;
                longthin_count = 0;
                thin_count = 0;
                stubby_count = 0;
                branched_count = 0;
                if(exist(filoFolder, 'dir'))
                    filoFiles = dir(fullfile(filoFolder, '*.tif'));
                    for i = 1:length(filoFiles)
                        filename = filoFiles(i).name;
                        % do something with the filename
                        [~, filename_no_ext, ~] = fileparts(filename);
                        spine_dendrite_type_label_tmp(str2double(filename_no_ext),1) = 1;
                        ccdendrite_label = dendrite_score(dendrite_score(:,1) == str2double(filename_no_ext), 2);
                        if(ccdendrite_label ~= 0)
                            dendrite_each_spine_type_number_tmp(ccdendrite_label, 1) = dendrite_each_spine_type_number_tmp(ccdendrite_label, 1) + 1;
                            ccdendrite_label = ccdendrite_label + max_dendrite_label;
    
                        end
                        spine_dendrite_type_label_tmp(str2double(filename_no_ext),2) = ccdendrite_label;
                        filo_count = filo_count + 1;
                    end
                end
                if(exist(mushFolder, 'dir'))
                    mushFiles = dir(fullfile(mushFolder, '*.tif'));
                    for i = 1:length(mushFiles)
                        filename = mushFiles(i).name;
                        % do something with the filename
                        [~, filename_no_ext, ~] = fileparts(filename);
                        spine_dendrite_type_label_tmp(str2double(filename_no_ext),1) = 2;
                        ccdendrite_label = dendrite_score(dendrite_score(:,1) == str2double(filename_no_ext), 2);
                        if(ccdendrite_label ~= 0)
                            dendrite_each_spine_type_number_tmp(ccdendrite_label, 2) = dendrite_each_spine_type_number_tmp(ccdendrite_label, 2) + 1;
                            ccdendrite_label = ccdendrite_label + max_dendrite_label;
                        end
                        spine_dendrite_type_label_tmp(str2double(filename_no_ext),2) = ccdendrite_label;
                        mush_count = mush_count + 1;
                    end
                end
                if(exist(longthinFolder, 'dir'))
                    longthinFiles = dir(fullfile(longthinFolder, '*.tif'));
                    for i = 1:length(longthinFiles)
                        filename = longthinFiles(i).name;
                        % do something with the filename
                        [~, filename_no_ext, ~] = fileparts(filename);
                        spine_dendrite_type_label_tmp(str2double(filename_no_ext),1) = 3;
                        ccdendrite_label = dendrite_score(dendrite_score(:,1) == str2double(filename_no_ext), 2);
                        if(ccdendrite_label ~= 0)
                            dendrite_each_spine_type_number_tmp(ccdendrite_label, 3) = dendrite_each_spine_type_number_tmp(ccdendrite_label, 3) + 1;
                            ccdendrite_label = ccdendrite_label + max_dendrite_label;
                        end
                        spine_dendrite_type_label_tmp(str2double(filename_no_ext),2) = ccdendrite_label;
                        longthin_count = longthin_count + 1;
                    end
                end
                if(exist(thinFolder, 'dir'))
                    thinFiles = dir(fullfile(thinFolder, '*.tif'));
                    for i = 1:length(thinFiles)
                        filename = thinFiles(i).name;
                        % do something with the filename
                        [~, filename_no_ext, ~] = fileparts(filename);
                        spine_dendrite_type_label_tmp(str2double(filename_no_ext),1) = 4;
                        ccdendrite_label = dendrite_score(dendrite_score(:,1) == str2double(filename_no_ext), 2);
                        if(ccdendrite_label ~= 0)
                            dendrite_each_spine_type_number_tmp(ccdendrite_label, 4) = dendrite_each_spine_type_number_tmp(ccdendrite_label, 4) + 1;
                            ccdendrite_label = ccdendrite_label + max_dendrite_label;
                        end
                        spine_dendrite_type_label_tmp(str2double(filename_no_ext),2) = ccdendrite_label;
                        thin_count = thin_count + 1;
                    end
                end
                if(exist(stubbyFolder, 'dir'))
                    stubbyFiles = dir(fullfile(stubbyFolder, '*.tif'));
                    for i = 1:length(stubbyFiles)
                        filename = stubbyFiles(i).name;
                        % do something with the filename
                        [~, filename_no_ext, ~] = fileparts(filename);
                        spine_dendrite_type_label_tmp(str2double(filename_no_ext),1) = 5;
                        ccdendrite_label = dendrite_score(dendrite_score(:,1) == str2double(filename_no_ext), 2);
                        if(ccdendrite_label ~= 0)
                            dendrite_each_spine_type_number_tmp(ccdendrite_label, 5) = dendrite_each_spine_type_number_tmp(ccdendrite_label, 5) + 1;
                            ccdendrite_label = ccdendrite_label + max_dendrite_label;
                        end
                        spine_dendrite_type_label_tmp(str2double(filename_no_ext),2) = ccdendrite_label;
                        stubby_count = stubby_count + 1;
                    end
                end
                % go over the branched spine folder and record the type as 6
                if(exist(branchedFolder, 'dir'))
                    branchedFiles = dir(fullfile(branchedFolder, '*.tif.jpg'));
                    for i = 1:length(branchedFiles)
                        filename = branchedFiles(i).name;
                        % do something with the filename
                        [~, filename_with_tif, ~] = fileparts(filename);
                        [~, filename_no_ext, ~] = fileparts(filename_with_tif);
                        % find the closest spine in the single_spine_coordinates
                        cur_spine_coor = branched_spine_coordinates(str2double(filename_no_ext),:);
                        selected_single_spines = knnsearch(single_spine_coordinates,cur_spine_coor(:)');
                        ccdendrite_label = dendrite_score(dendrite_score(:,1) == selected_single_spines, 2);
                        if(ccdendrite_label ~= 0)
                            dendrite_each_spine_type_number_tmp(ccdendrite_label, 6) = dendrite_each_spine_type_number_tmp(ccdendrite_label, 6) + 1;
                            ccdendrite_label = ccdendrite_label + max_dendrite_label;
                        end
                        branched_count = branched_count + 1;
                    end
                end
                max_dendrite_label = max_dendrite_label + cur_max_dendrite_label;
                dendrite_each_spine_type_number = [dendrite_each_spine_type_number; dendrite_each_spine_type_number_tmp];
                spine_dendrite_type_label = [spine_dendrite_type_label; spine_dendrite_type_label_tmp];
                dendrite_label_idx = label2idx(dendrite_score(:,2));
                dendrite_length_tmp = zeros(length(dendrite_label_idx), 1);
                dendrite_radius_tmp = zeros(length(dendrite_label_idx), 1);
                for n = 1:length(dendrite_label_idx)
                    curID = dendrite_label_idx{n};
                    if(~isempty(curID))
                        dendrite_length_tmp(n) = dendrite_score(curID(1), 3);
                        dendrite_radius_tmp(n) = dendrite_score(curID(1), 4);
                    end
                end
                
                dendrite_length = [dendrite_length; dendrite_length_tmp]; % nm
                dendrite_radius = [dendrite_radius; dendrite_radius_tmp]; % nm
                neuron_total_type_number = [neuron_total_type_number; [filo_count, mush_count, longthin_count, thin_count, stubby_count, branched_count]];
            end
        end
    end
    layer2_summary_level1_score.dendrite_each_spine_type_number = dendrite_each_spine_type_number;
    layer2_summary_level1_score.spine_dendrite_type_label = spine_dendrite_type_label;
    layer2_summary_level1_score.dendrite_length = dendrite_length;
    layer2_summary_level1_score.dendrite_radius = dendrite_radius;
    layer2_summary_level1_score.neuron_total_type_number = neuron_total_type_number;
    disp(ttc)
end