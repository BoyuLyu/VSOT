%% test the methods for level-1 segmentation 
branchlist = ["D5_Apical_Spines","D5_Branch_1","D5_Branch_2","D5_Branch_3","D5_Branch_4","D5_Branch_5"];
for j = 1:length(branchlist)
    % disp(j)
    current_branch = branchlist{j};
    summary_result_folder = [data_path,'/performance_testing_level_1_segmentation/results/summary_of_results/', current_branch];
    rootFolder_0 = [data_path, '/performance_testing_level_1_segmentation/dataset/', current_branch];
    % current_branch = 'D5_Apical_Spines';
    % rootFolder = ['/work/boyu/EM_astrocyte/test_segmentation_samples/dendrite_spine_segmentation/spinetools_and_our_and_method_3_biotech_result/', current_branch];
    color_result_folder = [data_path, '/performance_testing_level_1_segmentation/results/color_results/', current_branch];
    if ~exist(color_result_folder, 'dir')
        mkdir(color_result_folder)
    end
    [Pts,Tri] = read_off(fullfile(rootFolder_0, ['dendrite.off']));
    % bin_img = tiffreadVolume(fullfile(rootFolder, [current_branch, '_dendrite_volume.tif.tif'])) > 0;

    Tri = Tri';
    Pts = Pts';
    colormapx = [55,126,184, 255; 228,26,28,255; 77,175,74, 255];
    node_colorMap = zeros(size(Pts,1),4);
    % face_centers = (Pts(Tri(:,1),:) + Pts(Tri(:,2),:)  + Pts(Tri(:,3),:) )/3;
    % dict_points = find(bin_img(:) == 1);
    % [lenx, leny, lenz] = size(bin_img);
    % [dict_pointsx, dict_pointsy, dict_pointsz] = ind2sub([lenx, leny, lenz], dict_points);

    % dict_pointsxyz = [dict_pointsx(:), dict_pointsy(:), dict_pointsz(:)];
    % idx = knnsearch(dict_pointsxyz, face_centers, 'K', 1);


    label_face_gt = table2array(readtable(fullfile(summary_result_folder, 'gt_face_classification_result.txt')));
    label_face_us = table2array(readtable(fullfile(summary_result_folder, 'our_face_classification_result.txt')));
    label_face_neurd = table2array(readtable(fullfile(summary_result_folder, 'neurd_face_classification_result.txt')));
    label_face_spinetool = table2array(readtable(fullfile(summary_result_folder, 'spinetool_face_classification_result.txt')));
    label_face_morph = table2array(readtable(fullfile(summary_result_folder, 'morph_face_classification_result.txt')));

    face_part1 = find(label_face_gt == 1);
    face_part2 = find(label_face_gt == 2);
    vert_part1 = Tri(face_part1(:),:);
    vert_part1 = unique(vert_part1(:));
    vert_part2 = Tri(face_part2(:),:);
    vert_part2 = unique(vert_part2(:));
    cutCycle = intersect(vert_part1, vert_part2);
    node_colorMap(vert_part1(:),:) = repmat(colormapx(1,:), [length(unique(vert_part1(:))),1]);
    node_colorMap(vert_part2(:),:) = repmat(colormapx(2,:), [length(unique(vert_part2(:))), 1]);
    % node_colorMap(vert_part2(:),:) = repmat(colormapx(1,:), [length(unique(vert_part2(:))),1]);
    fid = fopen([fullfile(color_result_folder, ['gt','_color_ds.off'])],'wt');
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

    %write gt_surface color


    % generate our own result of surface 
    % us_shaft = tiffreadVolume(fullfile(rootFolder,"shaft_segmentation.tif")) > 0;
    % se1 = strel('sphere',1);
    % us_shaft = imdilate(us_shaft, se1);
    % us_spine = bin_img - us_shaft > 0;
    % us_combined = us_shaft*1 + us_spine*2;
    % label_face_us = us_combined(dict_points(idx));
    face_part1 = find(label_face_us == 1);
    face_part2 = find(label_face_us == 2);
    vert_part1 = Tri(face_part1(:),:);
    vert_part1 = unique(vert_part1(:));
    vert_part2 = Tri(face_part2(:),:);
    vert_part2 = unique(vert_part2(:));
    cutCycle = intersect(vert_part1, vert_part2);
    node_colorMap(vert_part1(:),:) = repmat(colormapx(1,:), [length(unique(vert_part1(:))),1]);
    node_colorMap(vert_part2(:),:) = repmat(colormapx(2,:), [length(unique(vert_part2(:))), 1]);
    % node_colorMap(vert_part2(:),:) = repmat(colormapx(1,:), [length(unique(vert_part2(:))),1]);
    fid = fopen([fullfile(color_result_folder, ['our_method','_color_ds.off'])],'wt');
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


    % organize the result from spinetools
    % read in all the points that are detected as spines
    % sptools_result_path = fullfile(rootFolder, 'spinetools_result/');
    % all_file_name = dir(fullfile(sptools_result_path, '*.off'));
    % points_spine = [];
    % for i = 1:length(all_file_name)
    %     [Pts2,Tri2] = read_off(fullfile(sptools_result_path, all_file_name(i).name));
    %     Tri2 = Tri2';
    %     Pts2 = Pts2';
    %     points_spine = [points_spine; Pts2];
    % end
    % vert_idx = knnsearch(Pts, points_spine, 'K', 1);
    % vert_part1 = unique(vert_idx);
    % vert_part2 = setdiff((1:size(Pts)), vert_part1);




    face_part1 = find(label_face_spinetool == 1);
    face_part2 = find(label_face_spinetool == 2);
    vert_part1 = Tri(face_part1(:),:);
    vert_part1 = unique(vert_part1(:));
    vert_part2 = Tri(face_part2(:),:);
    vert_part2 = unique(vert_part2(:));
    node_colorMap(vert_part1(:),:) = repmat(colormapx(1,:), [length(unique(vert_part1(:))),1]);
    node_colorMap(vert_part2(:),:) = repmat(colormapx(2,:), [length(unique(vert_part2(:))), 1]);
    fid = fopen([fullfile(color_result_folder, ['spinetools_method','_color_ds.off'])],'wt');
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






    face_part1 = find(label_face_morph == 1);
    face_part2 = find(label_face_morph == 2);
    vert_part1 = Tri(face_part1(:),:);
    vert_part1 = unique(vert_part1(:));
    vert_part2 = Tri(face_part2(:),:);
    vert_part2 = unique(vert_part2(:));
    node_colorMap(vert_part1(:),:) = repmat(colormapx(1,:), [length(unique(vert_part1(:))),1]);
    node_colorMap(vert_part2(:),:) = repmat(colormapx(2,:), [length(unique(vert_part2(:))), 1]);
    fid = fopen([fullfile(color_result_folder, ['morph_method','_color_ds.off'])],'wt');
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



    face_part1 = find(label_face_neurd == 1);
    face_part2 = find(label_face_neurd == 2);
    vert_part1 = Tri(face_part1(:),:);
    vert_part1 = unique(vert_part1(:));
    vert_part2 = Tri(face_part2(:),:);
    vert_part2 = unique(vert_part2(:));
    node_colorMap(vert_part1(:),:) = repmat(colormapx(1,:), [length(unique(vert_part1(:))),1]);
    node_colorMap(vert_part2(:),:) = repmat(colormapx(2,:), [length(unique(vert_part2(:))), 1]);
    fid = fopen([fullfile(color_result_folder, ['neurd_method','_color_ds.off'])],'wt');
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