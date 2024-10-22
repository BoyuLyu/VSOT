function [F1score, precision, recall, adjjaccardx] = level_1_performance_testing_helper(gt, target, surfaceTri, triArea)
%% F1 score measure in terms of each triangle
    TP_id = find((gt == 2) & (target == 2));
    FP_id = find((gt == 1) & (target == 2));
    FN_id = find((gt == 2) & (target == 1));
    F1score = (2*sum(triArea(TP_id)))/(sum(2*triArea(TP_id)) + sum(triArea(FP_id)) + sum(triArea(FN_id)));
%%
    PP_id = find(target == 2);
    P_id = find(gt == 2);
    precision = sum(triArea(TP_id))/ sum(triArea(PP_id));
    recall = sum(triArea(TP_id))/ sum(triArea(P_id));
%% mean IOU
    % basic idea is 2*|A & B|/(|A| + |B|)
    C = 0;
    U = 0;
    % measure the error with respect to each dendrite spine
    % first group gt into clusters
    % then group target surface into clusters
    big_graph_list = [surfaceTri(:,1:2), [1:size(surfaceTri)]'; surfaceTri(:,2:3), [1:size(surfaceTri)]'; ...
        surfaceTri(:,[3,1]), [1:size(surfaceTri)]'];
    big_graph_list(:,1:2) = sort(big_graph_list(:,1:2),2);
    big_graph_list = sortrows(big_graph_list,[1,2]);
    face_graph_list = [big_graph_list(1:2:end,3), big_graph_list(2:2:end,3)];
    
    gt_spines_fid = find(gt == 2);
    gt_spines_fid_binary = gt == 2;
    gt_spines_graph_list = face_graph_list(ismember(face_graph_list(:,1), gt_spines_fid) & ismember(face_graph_list(:,2), gt_spines_fid),:);
    gt_spines_graph_g = graph(gt_spines_graph_list(:,1), gt_spines_graph_list(:,2));
    gt_spines_bin = conncomp(gt_spines_graph_g);
    if(length(gt_spines_bin) < length(gt_spines_fid_binary))
        gt_spines_bin = [gt_spines_bin, zeros(1, length(gt_spines_fid_binary) - length(gt_spines_bin))];
        gt_spines_bin = gt_spines_bin(:).*double(gt_spines_fid_binary);
    else
        gt_spines_bin = gt_spines_bin(:).*double(gt_spines_fid_binary);
    end
    gt_spines_bin_idx = label2idx(gt_spines_bin);

    gt_spines_bin_idx_len = cellfun(@length, gt_spines_bin_idx);
    gt_spines_bin_idx(gt_spines_bin_idx_len == 0) = [];


    target_spines_fid = find(target == 2);
    target_spines_fid_binary = target == 2;
    target_spines_graph_list = face_graph_list(ismember(face_graph_list(:,1), target_spines_fid) & ismember(face_graph_list(:,2), target_spines_fid), :);
    target_spines_graph_g = graph(target_spines_graph_list(:,1), target_spines_graph_list(:,2));
    target_spines_bin = conncomp(target_spines_graph_g);
    if(length(target_spines_bin) < length(target_spines_fid_binary))
        target_spines_bin = [target_spines_bin, zeros(1, length(target_spines_fid_binary) - length(target_spines_bin))];
        target_spines_bin = target_spines_bin(:).*double(target_spines_fid_binary);
    else
        target_spines_bin = target_spines_bin(:).*double(target_spines_fid_binary);
    end
    target_spines_bin_idx = label2idx(target_spines_bin);

    target_searched_list = [];
    for i = 1:length(gt_spines_bin_idx)
        tmpid = gt_spines_bin_idx{i};
        label_target = mode(target_spines_bin(tmpid));
        if(label_target ~= 0)
            Cid = intersect(tmpid, target_spines_bin_idx{label_target});
            Uid = union(tmpid, target_spines_bin_idx{label_target});
            C = C + sum(triArea(Cid));
            U = U + sum(triArea(Uid));
            target_searched_list = [target_searched_list;label_target];
        else
            U = U + sum(triArea(tmpid));
        end
    end
    target_spines_bin_idx(target_searched_list) = [];
    target_spines_bin_idx_len = cellfun(@length, target_spines_bin_idx);
    target_spines_bin_idx(target_spines_bin_idx_len == 0) = [];
    for  i = 1:length(target_spines_bin_idx)
        tmpid = target_spines_bin_idx{i};
        U = U + sum(triArea(tmpid));

    end
    adjjaccardx = C/U;







end