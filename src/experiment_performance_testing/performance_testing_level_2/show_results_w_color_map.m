function show_results_w_color_map( offFolder, result_path,annotationFolder_gt_curves)
    listx = dir([offFolder, '/*.off']);
    tamada_result_folder = fullfile(result_path, 'tamada_result_total/');
    our_method_cut_result_folder = fullfile(result_path, 'vsot_result_total');
    ofer_output_folder = fullfile(result_path, 'ofer_total/');
    dorkenwalk_result_folder = fullfile(result_path,'dorkenwald_total/' );
    color_result_folder = fullfile(result_path, 'color_map_result');
    if(~exist(color_result_folder, 'dir'))
        mkdir(color_result_folder);
    end

    colormapx = [55,126,184, 255; 228,26,28,255; 77,175,74, 255];
    for j = 1:length(listx)
        % j = 8 for testing
        namex = listx(j).name;
        namex = namex(1:end-4);
        disp(namex)


        [Pts,Tri] = read_off(fullfile(offFolder, [namex , '.off']));
        Tri = Tri';
        Pts = Pts';
        node_colorMap = zeros(size(Pts,1),4);
        if(exist(fullfile(annotationFolder_gt_curves, [namex,'.cut.txt']), 'file'))
            gt_cut = readtable(fullfile(annotationFolder_gt_curves, [namex,'.cut.txt']));
            gt_cut = table2array(gt_cut);
            face_gt = load(fullfile(annotationFolder_gt_curves, [namex,'.face_label.mat']));
            face_gt = face_gt.head_neck_label;
            face_gt = face_gt(:);
            face_part1 = find(face_gt == 1);
            face_part2 = find(face_gt == 2);
            vert_part1 = Tri(face_part1(:),:);
            vert_part1 = unique(vert_part1(:));
            vert_part2 = Tri(face_part2(:),:);
            vert_part2 = unique(vert_part2(:));
            cutCycle = intersect(vert_part1, vert_part2);
            node_colorMap(vert_part1(:),:) = repmat(colormapx(2,:), [length(unique(vert_part1(:))),1]);
            node_colorMap(vert_part2(:),:) = repmat(colormapx(3,:), [length(unique(vert_part2(:))), 1]);

            fid = fopen([fullfile(color_result_folder, [namex, 'gt','_color_ds.off'])],'wt');
            if( fid==-1 )
                error('Can''t open the file.');
            end
                        % header
            fprintf(fid, 'COFF\n');
            fprintf(fid, '%d %d 0\n', size(Pts,1), size(Tri,1));
            % write the points & faces
            fprintf(fid, '%f %f %f %d %d %d %d\n', [Pts, node_colorMap]');
            fprintf(fid, '3 %d %d %d\n', Tri'-1);
            fclose(fid);

            try 
            if(exist(fullfile(tamada_result_folder, [namex,'.cut.txt']), 'file'))

                tamada_cut = readtable(fullfile(tamada_result_folder, [namex,'.cut.txt']));
                tamada_cut = table2array(tamada_cut);
                face_label = load(fullfile(tamada_result_folder, [namex,'.face_label.mat']));
                face_label = face_label.head_neck_label;
                face_label = face_label(:);
                face_part1 = find(face_label == 1);
                face_part2 = find(face_label == 2);
                vert_part1 = Tri(face_part1(:),:);
                vert_part1 = unique(vert_part1(:));
                vert_part2 = Tri(face_part2(:),:);
                vert_part2 = unique(vert_part2(:));
                cutCycle = intersect(vert_part1, vert_part2);
                node_colorMap(vert_part1(:),:) = repmat(colormapx(2,:), [length(unique(vert_part1(:))),1]);
                node_colorMap(vert_part2(:),:) = repmat(colormapx(3,:), [length(unique(vert_part2(:))), 1]);
                fid = fopen([fullfile(color_result_folder, [namex,'tamada','_color_ds.off'])],'wt');
                if( fid==-1 )
                    error('Can''t open the file.');
                end
                % header
                fprintf(fid, 'COFF\n');
                fprintf(fid, '%d %d 0\n', size(Pts,1), size(Tri,1));
                % write the points & faces
                fprintf(fid, '%f %f %f %d %d %d %d\n', [Pts, node_colorMap]');
                fprintf(fid, '3 %d %d %d\n', Tri'-1);
                fclose(fid);

            end
            
            catch ME
                continue;
            end
            try 
    
            if(exist(fullfile(our_method_cut_result_folder,[namex,'.cut.txt']), 'file'))
    
                our_method_cut = readtable(fullfile(our_method_cut_result_folder,[namex,'.cut.txt']));
                our_method_cut = table2array(our_method_cut);
                face_label = load(fullfile(our_method_cut_result_folder, [namex,'.face_label.mat']));
                face_label = face_label.head_neck_label;
                face_label = face_label(:);
                face_part1 = find(face_label == 1);
                face_part2 = find(face_label == 2);
                vert_part1 = Tri(face_part1(:),:);
                vert_part1 = unique(vert_part1(:));
                vert_part2 = Tri(face_part2(:),:);
                vert_part2 = unique(vert_part2(:));
                cutCycle = intersect(vert_part1, vert_part2);
                node_colorMap(vert_part1(:),:) = repmat(colormapx(2,:), [length(unique(vert_part1(:))),1]);
                node_colorMap(vert_part2(:),:) = repmat(colormapx(3,:), [length(unique(vert_part2(:))), 1]);
                fid = fopen([fullfile(color_result_folder, [namex,'vsot','_color_ds.off'])],'wt');
                if( fid==-1 )
                    error('Can''t open the file.');
                end
                % header
                fprintf(fid, 'COFF\n');
                fprintf(fid, '%d %d 0\n', size(Pts,1), size(Tri,1));
                % write the points & faces
                fprintf(fid, '%f %f %f %d %d %d %d\n', [Pts, node_colorMap]');
                fprintf(fid, '3 %d %d %d\n', Tri'-1);
                fclose(fid);
            end
            catch ME
                continue;
            end
            try 
                if(exist(fullfile(ofer_output_folder,[namex,'.cut.txt']), 'file'))
        
                    ofer_method_cut = readtable(fullfile(ofer_output_folder,[namex,'.cut.txt']));
                    ofer_method_cut = table2array(ofer_method_cut);
                    face_label = load(fullfile(ofer_output_folder, [namex,'.face_label.mat']));
                    face_label = face_label.head_neck_label;
                    face_label = face_label(:);
                    face_part1 = find(face_label == 1);
                    face_part2 = find(face_label == 2);
                    vert_part1 = Tri(face_part1(:),:);
                    vert_part1 = unique(vert_part1(:));
                    vert_part2 = Tri(face_part2(:),:);
                    vert_part2 = unique(vert_part2(:));
                    cutCycle = intersect(vert_part1, vert_part2);
                    node_colorMap(vert_part1(:),:) = repmat(colormapx(2,:), [length(unique(vert_part1(:))),1]);
                    node_colorMap(vert_part2(:),:) = repmat(colormapx(3,:), [length(unique(vert_part2(:))), 1]);
                    fid = fopen([fullfile(color_result_folder, [namex,'ofer','_color_ds.off'])],'wt');
                    if( fid==-1 )
                        error('Can''t open the file.');
                    end
                    % header
                    fprintf(fid, 'COFF\n');
                    fprintf(fid, '%d %d 0\n', size(Pts,1), size(Tri,1));
                    % write the points & faces
                    fprintf(fid, '%f %f %f %d %d %d %d\n', [Pts, node_colorMap]');
                    fprintf(fid, '3 %d %d %d\n', Tri'-1);
                    fclose(fid);
                end
            catch ME
                continue;
            end
            try 
                if(exist(fullfile(dorkenwalk_result_folder,[namex,'.cut.txt']), 'file'))
        
                    dorkenwald_method_cut = readtable(fullfile(dorkenwalk_result_folder,[namex,'.cut.txt']));
                    dorkenwald_method_cut = table2array(dorkenwald_method_cut);
                    face_label = load(fullfile(dorkenwalk_result_folder, [namex,'.face_label.mat']));
                    face_label = face_label.head_neck_label;
                    face_label = face_label(:);
                    face_part1 = find(face_label == 1);
                    face_part2 = find(face_label == 2);
                    vert_part1 = Tri(face_part1(:),:);
                    vert_part1 = unique(vert_part1(:));
                    vert_part2 = Tri(face_part2(:),:);
                    vert_part2 = unique(vert_part2(:));
                    cutCycle = intersect(vert_part1, vert_part2);
                    node_colorMap(vert_part1(:),:) = repmat(colormapx(2,:), [length(unique(vert_part1(:))),1]);
                    node_colorMap(vert_part2(:),:) = repmat(colormapx(3,:), [length(unique(vert_part2(:))), 1]);
                    fid = fopen([fullfile(color_result_folder, [namex,'dorkenwald','_color_ds.off'])],'wt');
                    if( fid==-1 )
                        error('Can''t open the file.');
                    end
                    % header
                    fprintf(fid, 'COFF\n');
                    fprintf(fid, '%d %d 0\n', size(Pts,1), size(Tri,1));
                    % write the points & faces
                    fprintf(fid, '%f %f %f %d %d %d %d\n', [Pts, node_colorMap]');
                    fprintf(fid, '3 %d %d %d\n', Tri'-1);
                    fclose(fid);
                end
            catch ME
                continue;
            end

        end
    end
