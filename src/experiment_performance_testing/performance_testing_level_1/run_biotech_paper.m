%% biotechnology paper that uses erosion and then dilation to obtain the dendrite spines
% implementation in MATLAB
clear
current_branch = 'D5_Apical_Spines';
rootFolder = ['../data/performance_testing_level_1_segmentation/results/bio_tech_paper_result/',current_branch];
bin_img = tiffreadVolume(fullfile(rootFolder,[current_branch, '_dendrite_volume.tif.tif'])) > 0;
se0 = strel('sphere', 3);
bin_img2 = imopen(bin_img, se0);
bin_img2_roi = bwlabeln(bin_img2);
bin_img2_roi_idx = label2idx(bin_img2_roi);
bin_img2 = bin_img2.*0;
len_cell = cellfun(@length, bin_img2_roi_idx);
bin_img2(bin_img2_roi_idx{(len_cell == max(len_cell))}) = 1;
se1 = strel('sphere', 1);
bin_img3 = imdilate(bin_img2, se1);
% figure; volshow(bin_img3)
spine_region = (bin_img - bin_img3) > 0;
spine_region_roi = bwlabeln(spine_region);
spine_region_roi_idx = label2idx(spine_region_roi);
spine_region_roi_idx(cellfun(@length, spine_region_roi_idx) < 3) = [];
spine_region_roi_idx = spine_region_roi_idx(:);
new_combined = zeros(size(spine_region));
new_combined = double(bin_img);
new_combined(cell2mat(spine_region_roi_idx)) = 2;
% figure; volshow(new_combined)
[Pts,Tri] = read_off(fullfile(rootFolder, ['dendrite.off']));
Tri = Tri';
Pts = Pts';
colormapx = [255, 64, 64, 255;64, 255, 255,255;64, 255, 64, 255];
node_colorMap = zeros(size(Pts,1),4);
face_centers = (Pts(Tri(:,1),:) + Pts(Tri(:,2),:)  + Pts(Tri(:,3),:) )/3;
%find the nearest faces to t
dict_points = find(bin_img(:) == 1);
[lenx, leny, lenz] = size(bin_img);
[dict_pointsx, dict_pointsy, dict_pointsz] = ind2sub([lenx, leny, lenz], dict_points);

dict_pointsxyz = [dict_pointsx(:), dict_pointsy(:), dict_pointsz(:)];
idx = knnsearch(dict_pointsxyz, face_centers, 'K', 1);
label_face = new_combined(dict_points(idx));
face_part1 = find(label_face == 1);
face_part2 = find(label_face == 2);
vert_part1 = Tri(face_part1(:),:);
vert_part1 = unique(vert_part1(:));
vert_part2 = Tri(face_part2(:),:);
vert_part2 = unique(vert_part2(:));
cutCycle = intersect(vert_part1, vert_part2);
node_colorMap(vert_part1(:),:) = repmat(colormapx(3,:), [length(unique(vert_part1(:))),1]);
node_colorMap(cutCycle,:) = repmat(colormapx(2,:), [length(unique(cutCycle)), 1]);
node_colorMap(vert_part2(:),:) = repmat(colormapx(1,:), [length(unique(vert_part2(:))),1]);
fid = fopen([fullfile(rootFolder, ['method_3','_color_ds.off'])],'wt');
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

