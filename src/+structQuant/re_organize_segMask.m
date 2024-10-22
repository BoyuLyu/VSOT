function [seg_mask_small, datax_small] = re_organize_segMask(seg_mask_small, datax_big)
    label_small = unique(seg_mask_small(:));
    label_small(label_small == 0) = [];
    seg_mask_small_idx = label2idx(seg_mask_small);
    seg_mask_small_idx = seg_mask_small_idx(:);
    keys_label_small = cell(length(label_small),1);
    for i = 1:length(label_small)
        keys_label_small{i} = datax_big(label_small(i));
    end
    [unique_keys, ia, ic] = unique(keys_label_small);
    unique_values = cell(length(unique_keys),1);
    seg_mask_small_idx2 = cell(length(unique_keys),1);
    seg_mask_small = seg_mask_small.*0;
    for i = 1:length(unique_keys)
        unique_values{i}=  i;
        seg_mask_small_idx2{i} = cell2mat(seg_mask_small_idx(label_small(ic == i)));
        seg_mask_small(seg_mask_small_idx2{i}) = i;
    end
    datax_small = containers.Map(unique_keys,unique_values);




    





















end