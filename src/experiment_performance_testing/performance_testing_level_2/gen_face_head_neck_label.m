function gen_face_head_neck_label(offFolder, an1, an2, an3, tifFolder)
    % for the two faces generated after splitting along the cut, 
    % in any part, if the shortest distance to the dendritic shaft is smaller, the part will be denoted as 1 (neck part), otherwise 2 (head part)
    % the final output will be a list of labels for each face

    list1 = dir(fullfile(an1, '*.cut.txt'));
    name1 = {list1.name};
    for i = 1:length(name1)
        tmp = name1{i};
        tmp = strsplit(tmp, '.');
        tmp = tmp{1};
        name1{i} = tmp;
    end

    for i = 1:length(name1)
        name_file = name1{i};
        labelx = tiffreadVolume(fullfile(tifFolder, [name_file, '.tif']));
        labelxidx = label2idx(labelx);
        [lenx, leny, lenz] = size(labelx);
        % [SpineCoorx,SpineCoory, SpineCoorz]  = ind2sub([lenx, leny, lenz], labelxidx{2});
        [ShaftCoorx,ShaftCoory, ShaftCoorz]  = ind2sub([lenx, leny, lenz], labelxidx{1});
        ShaftCoorx = ShaftCoorx * 2;
        ShaftCoory = ShaftCoory * 2;
        ShaftCoorz = ShaftCoorz * 5;
        [Pts,Tri] = read_off(fullfile(offFolder, [name_file , '.off']));
        Tri = Tri';
        Pts = Pts';
        Tri_center = (Pts(Tri(:,1),:) + Pts(Tri(:,2),:) + Pts(Tri(:,3),:))/3;
        edge_graph = [Tri(:,1), Tri(:,2), [1:size(Tri,1)]';Tri(:,1), Tri(:,3),[1:size(Tri,1)]';Tri(:,2), Tri(:,3),[1:size(Tri,1)]'];
        edge_graph(:,1:2) = sort(edge_graph(:,1:2), 2);
        edge_graph = sortrows(edge_graph, [1,2]);
        % split the surface along the cut
        ann1_cut = readtable(fullfile(an1, [name_file, '.cut.txt']));
        ann1_cut = table2array(ann1_cut);
        ann2_cut = readtable(fullfile(an2, [name_file, '.cut.txt']));
        ann2_cut = table2array(ann2_cut);
        ann3_cut = readtable(fullfile(an3, [name_file, '.cut.txt']));
        ann3_cut = table2array(ann3_cut);
        ann_all = {ann1_cut, ann2_cut, ann3_cut};
        head_neck_label = zeros(length(Tri_center),3);
        for m = 1:3
            gt_cut = ann_all{m};
            face_label = zeros(size(Tri,1),1);
            cut_curve_edges = [gt_cut(1:end),[gt_cut(2:end);gt_cut(1)]];
            cut_curve_edges = sort(cut_curve_edges, 2);
            edge_graph_2 = edge_graph;
            edge_graph_2(ismember(edge_graph_2(:,1:2), cut_curve_edges, 'rows'),:) = [];
            face_G = graph(edge_graph_2(1:2:end,3), edge_graph_2(2:2:end,3));
            [bins, binsizes] = conncomp(face_G);
            [~, sortedID] = sort(binsizes, 'descend');
            if(length(binsizes) == 1)
                warning('incomplete separation %s',name_file);
            else
                part1 = find(bins == sortedID(1));
                part2 = find(bins == sortedID(2));
                % check the distance to the shaft
                dist2part1 = sqrt((Tri_center(part1,1)' - ShaftCoorx).^2 + (Tri_center(part1,2)' - ShaftCoory).^2 + (Tri_center(part1,3)' - ShaftCoorz).^2);
                dist2part2 = sqrt((Tri_center(part2,1)' - ShaftCoorx).^2 + (Tri_center(part2,2)' - ShaftCoory).^2 + (Tri_center(part2,3)' - ShaftCoorz).^2);
                if(min(dist2part1(:)) > min(dist2part2(:)))
                    % 1 is neck, 2 is head
                    head_neck_label(part1, m) = 2;
                    head_neck_label(part2, m) = 1;
                else
                    head_neck_label(part1, m) = 1;
                    head_neck_label(part2, m) = 2;
                end
            end
        end
        head_neck_label1 = head_neck_label(:,1);
        save(fullfile(an1, [name_file,'.face_label.mat']), 'head_neck_label1');
        head_neck_label2 = head_neck_label(:,2);
        save(fullfile(an2, [name_file,'.face_label.mat']), 'head_neck_label2');
        head_neck_label3 = head_neck_label(:,3);
        save(fullfile(an3, [name_file,'.face_label.mat']), 'head_neck_label3');
    %     figure; trisurf(Tri(head_neck_label(:,1) == 1,:), Pts(:,1), Pts(:,2), Pts(:,3), 'Facecolor', 'red'); ...
    % hold on; trisurf(Tri(head_neck_label(:,1) == 2,:), Pts(:,1), Pts(:,2), Pts(:,3), 'Facecolor', 'blue');
    %     hold on; plot3(Pts(gt_cut,1), Pts(gt_cut,2), Pts(gt_cut,3), 'LineWidth',10)
    end









end