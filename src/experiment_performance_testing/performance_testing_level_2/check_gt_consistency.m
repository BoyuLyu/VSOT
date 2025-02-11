function out_record = check_gt_consistency(off_folder, ann1, ann2, ann3)
    
    list = dir(fullfile(ann1, '*.cut.txt'));
    namex = {list.name};
    for i = 1:length(namex)
        tmp = namex{i};
        tmp = strsplit(tmp, '.');
        tmp = tmp{1};
        namex{i} = tmp;
    end

    IOU_measure = zeros(length(namex),3);

    for i = 1:length(namex)
% read in each .off file and check the cut index of all the three annotators
        name_file = namex{i};
        [Pts,Tri] = read_off(fullfile(off_folder, [name_file , '.off']));
        Tri = Tri';
        Pts = Pts';
        ann1_cut = readtable(fullfile(ann1, [name_file, '.cut.txt']));
        ann1_cut = table2array(ann1_cut);
        ann2_cut = readtable(fullfile(ann2, [name_file, '.cut.txt']));
        ann2_cut = table2array(ann2_cut);
        ann3_cut = readtable(fullfile(ann3, [name_file, '.cut.txt']));
        ann3_cut = table2array(ann3_cut);

        ann_all = {ann1_cut, ann2_cut, ann3_cut};
        column_id = 0;
        vectorMF3_1 = Pts(Tri(:,3),:) - Pts(Tri(:,1),:);
        vectorMF3_2 = Pts(Tri(:,3),:) - Pts(Tri(:,2),:);
        ss_tmp = cross(vectorMF3_1, vectorMF3_2,2);
        area_each_face = 1/2*(sqrt(ss_tmp(:,1).^2 + ss_tmp(:,2).^2 + ss_tmp(:,3).^2));
        load(fullfile(ann1, [name_file, '.face_label.mat'])); % head_neck_label1
        load(fullfile(ann2, [name_file, '.face_label.mat'])); % head_neck_label2
        load(fullfile(ann3, [name_file, '.face_label.mat'])); % head_neck_label3
        face_cut = {head_neck_label1, head_neck_label2, head_neck_label3};
        % % obtain the cut faces for each annotator
        % face_cut = cell(3,1);
        % edge_graph = [Tri(:,1), Tri(:,2), [1:size(Tri,1)]';Tri(:,1), Tri(:,3),[1:size(Tri,1)]';Tri(:,2), Tri(:,3),[1:size(Tri,1)]'];
        % edge_graph(:,1:2) = sort(edge_graph(:,1:2), 2);
        % edge_graph = sortrows(edge_graph, [1,2]);
        % for m = 1:3
        %     gt_cut = ann_all{m};
        %     face_label = zeros(size(Tri,1),1);
        %     cut_curve_edges = [gt_cut(1:end),[gt_cut(2:end);gt_cut(1)]];
        %     cut_curve_edges = sort(cut_curve_edges, 2);
        %     edge_graph_2 = edge_graph;
        %     edge_graph_2(ismember(edge_graph_2(:,1:2), cut_curve_edges, 'rows'),:) = [];
        %     face_G = graph(edge_graph_2(1:2:end,3), edge_graph_2(2:2:end,3));
        %     [bins, binsizes] = conncomp(face_G);
        %     [~, sortedID] = sort(binsizes, 'descend');
        %     if(length(binsizes) == 1)
        %         warning('incomplete separation %s',name_file);
        %     else
        %         part1 = find(bins == sortedID(1));
        %         part2 = find(bins == sortedID(2));
        %         face_label(part1) = 1;
        %         face_label(part2) = 2;
        %         face_cut{m} = face_label;
        %     end

        % end
        for m = 1:2
            for n = (m+1):3
                column_id = column_id + 1;
                % cut_1 = ann_all{m};
                % cut_2 = ann_all{n};
                % dist_min_1 = zeros(length(cut_1),1);
                % for k = 1:length(dist_min_1)
                %     dist_min_1(k) = min(vecnorm(Pts(cut_1(k),:) - Pts(cut_2,:),2,2));
                % end
                % dist_min_2 = zeros(length(cut_2),1);
                % for k = 1:length(dist_min_2)
                %     dist_min_2(k) = min(vecnorm(Pts(cut_2(k),:) - Pts(cut_1,:),2,2));
                % end
                % distance_2_center(i,column_id) = mean([dist_min_1;dist_min_2]);
                % check IOU of the two cuts

                % not sure whether the labels are correlated use the larger IOU as the output
                % IOU1 = sum(face_1 == face_2 & face_1 > 0 & face_2 > 0)/sum(face_1 > 0 | face_2 > 0);
                % face_2_2 = face_2;
                % face_2_2(face_2 == 1) = 2;
                % face_2_2(face_2 == 2) = 1;
                % IOU2 = sum(face_1 == face_2_2 & face_1 > 0 & face_2_2 > 0)/sum(face_1 > 0 | face_2_2 > 0);
                % IOU_measure(i,column_id) = max([IOU1, IOU2]);
                IOU_measure(i,column_id) = sum(area_each_face(face_cut{m} == face_cut{n} & face_cut{m} >0 & face_cut{n} > 0))/sum(area_each_face(face_cut{m} >0 | face_cut{n} > 0));
            end
        end
    end

    % Create a table with the names and the distance matrix
    out_record = array2table(IOU_measure, 'VariableNames', ...
        {'Ann1_Ann2', 'Ann1_Ann3', 'Ann2_Ann3'});
    out_record.Name = namex';
    
    % Reorder the columns to have 'Name' as the first column
    out_record = [out_record(:, end), out_record(:, 1:end-1)];


end
