function find_mean_curve(selectedID, out_record, offFolder, an1, an2, an3, gt_folder)
% for the cutIDs they are consistent, we can use the mean curve as the final output
refined_folder = fullfile(gt_folder, 'refined_cutID');
if(~exist(refined_folder, 'dir'))
    mkdir(refined_folder)
end
names = out_record.Name;
for i = 1:size(selectedID, 1)
    %It goes back to the surface (select the faces )

    name_file = names{i};
    load(fullfile(an1, [name_file, '.face_label.mat'])); % head_neck_label1
    load(fullfile(an2, [name_file, '.face_label.mat'])); % head_neck_label2
    load(fullfile(an3, [name_file, '.face_label.mat'])); % head_neck_label3
    face_label = {head_neck_label1, head_neck_label2, head_neck_label3};
    
    [Pts,Tri] = read_off(fullfile(offFolder, [name_file , '.off']));
    Tri = Tri';
    Pts = Pts';

    face_fin_label = zeros(size(Tri,1),1);
    if(sum(selectedID(i,2:end) == 3))
        part1 = find(face_label{1} == 1 & face_label{2} == 1 & face_label{3} == 1);
        face_fin_label(part1) = 1;
        part2 = find(face_label{1} == 2 & face_label{2} == 2 & face_label{3} == 2);
        face_fin_label(part2) = 2;
        [curve_vertex, face_fin_label] = find_mean_curve_sub(Tri, Pts, name_file, face_fin_label, part1, part2);

    else
        two_annotator = find(selectedID(i,2:end) == 1);
        part1 = find(face_label{two_annotator(1)} == 1 & face_label{two_annotator(2)} == 1);
        face_fin_label(part1) = 1;
        part2 = find(face_label{two_annotator(1)} == 2 & face_label{two_annotator(2)} == 2);
        face_fin_label(part2) = 2;
        [curve_vertex, face_fin_label] = find_mean_curve_sub(Tri, Pts,name_file, face_fin_label, part1, part2);

    end
    % save the refined cutID
    writetable(table(curve_vertex), fullfile(refined_folder, [name_file, '.cut.txt']));
    % save the face_label
    head_neck_label = face_fin_label;
    save(fullfile(refined_folder, [name_file, '.face_label.mat']), 'head_neck_label'); % neck is 1 head is 2


end