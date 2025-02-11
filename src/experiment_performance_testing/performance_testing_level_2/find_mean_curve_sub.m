function [curve_vertex, face_fin_label] = find_mean_curve_sub(Tri, Pts, name_file, face_fin_label, part1, part2)



    Tri_center = (Pts(Tri(:,1),:) + Pts(Tri(:,2),:) + Pts(Tri(:,3),:))/3;

    % edge_graph = [Tri(:,1), Tri(:,2), [1:size(Tri,1)]';Tri(:,1), Tri(:,3),[1:size(Tri,1)]';Tri(:,2), Tri(:,3),[1:size(Tri,1)]'];
    % edge_graph(:,1:2) = sort(edge_graph(:,1:2), 2);
    % edge_graph = sortrows(edge_graph, [1,2]);

    % remained faces
    part_center = setdiff(1:size(Tri,1), [part1;part2]);

    % classify the remained faces to two parts
    for m = 1:length(part_center)
        dist_part1 = vecnorm(Tri_center(part_center(m),:) - Tri_center(part1,:),2,2);
        dist_part2 = vecnorm(Tri_center(part_center(m),:) - Tri_center(part2,:),2,2);
        if(min(dist_part1) < min(dist_part2))
            face_fin_label(part_center(m)) = 1;
        else
            face_fin_label(part_center(m)) = 2;
        end
    end

    % obtain the cut between the two parts
    edge_all = [[Tri(:,1), Tri(:,2), [1:size(Tri,1)]']; [Tri(:,1), Tri(:,3), [1:size(Tri,1)]'];[Tri(:,3), Tri(:,2), [1:size(Tri,1)]']];
    edge_all(:,1:2) = sort(edge_all(:,1:2), 2);
    edge_all = sortrows(edge_all,[1,2]);
    edge_all(:,3) = face_fin_label(edge_all(:,3));
    a1 = [1:2:size(edge_all,1)];
    a2 = [2:2:size(edge_all,1)];
    id0 = find(edge_all(a1,3) ~= edge_all(a2,3));
    edge_all = edge_all(a1(id0),:);
    edge_all = edge_all(:,1:2);
    edge_all = sort(edge_all, 2);
    edge_all = unique(edge_all, 'rows');
    % re-order the vertices to make it a loop
    graph_loop = graph(edge_all(:,1), edge_all(:,2));
    allCycles = allcycles(graph_loop);
    cycle_length = cellfun(@length, allCycles);
    [~, idx] = max(cycle_length);
    curve_vertex = allCycles{idx};
    % sanity check if the curve is closed and can separate the two parts
    edge_graph = [Tri(:,1), Tri(:,2), [1:size(Tri,1)]';Tri(:,1), Tri(:,3),[1:size(Tri,1)]';Tri(:,2), Tri(:,3),[1:size(Tri,1)]'];
    edge_graph(:,1:2) = sort(edge_graph(:,1:2), 2);
    edge_graph = sortrows(edge_graph, [1,2]);
    edge_graph_2 = edge_graph;
    curve_vertex = curve_vertex(:);
    cut_curve_edges = [curve_vertex(1:end),[curve_vertex(2:end);curve_vertex(1)]];
    cut_curve_edges = sort(cut_curve_edges, 2);
    edge_graph_2(ismember(edge_graph_2(:,1:2), cut_curve_edges, 'rows'),:) = [];
    face_G = graph(edge_graph_2(1:2:end,3), edge_graph_2(2:2:end,3));
    [bins, binsizes] = conncomp(face_G);
    if(length(binsizes) == 1)
        warning('incomplete separation %s',name_file);
        cut_vertex = [];
        face_fin_label = [];
    end





















end