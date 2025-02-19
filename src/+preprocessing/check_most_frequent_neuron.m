% select the top neuron IDs from the dataset at 128nmx128nm x80nm

nueronId = load('/work/boyu/EM_astrocyte/astro_11_33_64_64_80/neuron_id_micron_dataset.mat');
nueronId_string = string(nueronId.neuron_id_microns);
for m = [11:33]
    targetFolder = ['/work/boyu/EM_astrocyte/astro_11_33_64_64_80/astro_', num2str(m),'_minnie65'];
    seg_mask_full = tiffreadVolume(fullfile(targetFolder,'segMaskFull.tif'));
    seg_mask_full_idx = label2idx(seg_mask_full);
    len_cell = cellfun(@length, seg_mask_full_idx);
    [sorted_len, sorted_id] = sort(len_cell, 'descend');
    fid = fopen(fullfile(targetFolder, 'seg_mapping.txt')); % Opening the file
    raw = fread(fid,inf); % Reading the contents
    str = char(raw'); % Transformation
    fclose(fid); % Closing the file
    datax = jsondecode(str); % Using the jsondecode function to parse JSON from string
    segRootIDfields = fieldnames(datax);
    segNewID = struct2array(datax);
    neuron_id_high_rank = cell(20, 1);
    count = 1;
    for i = 1:30
        newID_top = sorted_id(i);
        rootID_uint64 = segRootIDfields{segNewID == newID_top};
        tmp = rootID_uint64(2:end);
        if(ismember(tmp, nueronId_string))
                neuron_id_high_rank{count} = tmp;
                count = count + 1;
        end
    
        if  count == 21
            break
        end
    
    end
    writecell(neuron_id_high_rank, fullfile(targetFolder, 'top20_neuronID_spine_centered.txt'));
end