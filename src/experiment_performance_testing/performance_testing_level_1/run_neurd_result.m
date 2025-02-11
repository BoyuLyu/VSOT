%% test the methods for level-1 segmentation 
clear
current_branch = 'D5_Branch_1';
rootFolder = ['/work/boyu/EM_astrocyte/test_segmentation_samples/dendrite_spine_segmentation/spinetools_and_our_and_method_3_biotech_result/',current_branch];
[Pts,Tri] = read_off(fullfile(rootFolder, ['dendrite.off']));
bin_img = tiffreadVolume(fullfile(rootFolder, [current_branch, '_dendrite_volume.tif.tif'])) > 0;

Tri = Tri';
Pts = Pts';
colormapx = [255, 64, 64, 255;64, 255, 255,255;64, 255, 64, 255];
node_colorMap = zeros(size(Pts,1),4);
face_centers = (Pts(Tri(:,1),:) + Pts(Tri(:,2),:)  + Pts(Tri(:,3),:) )/3;
dict_points = find(bin_img(:) == 1);
[lenx, leny, lenz] = size(bin_img);
[dict_pointsx, dict_pointsy, dict_pointsz] = ind2sub([lenx, leny, lenz], dict_points);

dict_pointsxyz = [dict_pointsx(:), dict_pointsy(:), dict_pointsz(:)];
idx = knnsearch(dict_pointsxyz, face_centers, 'K', 1);
gt_shaft = tiffreadVolume(fullfile(rootFolder,[current_branch, '_shaft_volume.tif.tif'])) > 0;
gt_spine = bin_img - gt_shaft > 0;
gt_spine_roi = bwlabeln(gt_spine);
disp(max(gt_spine_roi(:)))
neurd_root_folder = '/work/boyu/EM_astrocyte/test_segmentation_samples/dendrite_spine_segmentation/neurd_results';
neurd_result_path = fullfile(neurd_root_folder, current_branch);
all_file_name = dir(fullfile(neurd_result_path, '*.off'));
points_spine = [];
for i = 1:length(all_file_name)
    [Pts2,Tri2] = read_off(fullfile(neurd_result_path, all_file_name(i).name));
    Tri2 = Tri2';
    Pts2 = Pts2';
    points_spine = [points_spine; Pts2];
end
vert_idx = knnsearch(Pts, points_spine, 'K', 1);
vert_part1 = unique(vert_idx);
vert_part2 = setdiff((1:size(Pts)), vert_part1);

node_colorMap(vert_part1(:),:) = repmat(colormapx(1,:), [length(unique(vert_part1(:))),1]);
% node_colorMap(cutCycle,:) = repmat(colormapx(2,:), [length(unique(cutCycle)), 1]);
node_colorMap(vert_part2(:),:) = repmat(colormapx(3,:), [length(unique(vert_part2(:))),1]);
fid = fopen([fullfile(rootFolder, ['neurd_method','_color_ds.off'])],'wt');
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