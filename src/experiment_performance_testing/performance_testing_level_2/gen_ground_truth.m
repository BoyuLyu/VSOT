function gen_ground_truth(gt_folder, an1, an2,an3,volume_w_shaft_folder)
% remove any samples that are not in all the annotation folders
% gt_folder = ''
an1 = fullfile(gt_folder, an1);
an2 = fullfile(gt_folder, an2);
an3 = fullfile(gt_folder, an3);
volume_w_shaft_folder = fullfile(gt_folder, volume_w_shaft_folder);
list1 = dir(fullfile(an1, '*.cut.txt'));
name1 = {list1.name};
for i = 1:length(name1)
    tmp = name1{i};
    tmp = strsplit(tmp, '.');
    tmp = tmp{1};
    name1{i} = tmp;
end
list2 = dir(fullfile(an2, '*.cut.txt'));
name2 = {list2.name};
for i = 1:length(name2)
    tmp = name2{i};
    tmp = strsplit(tmp, '.');
    tmp = tmp{1};
    name2{i} = tmp;
end
list3 = dir(fullfile(an3, '*.cut.txt'));
name3 = {list3.name};
for i = 1:length(name3)
    tmp = name3{i};
    tmp = strsplit(tmp, '.');
    tmp = tmp{1};
    name3{i} = tmp;
end
list4 = dir(fullfile(volume_w_shaft_folder, '*.tif'));
name4 = {list4.name};
for i = 1:length(name4)
    tmp = name4{i};
    tmp = strsplit(tmp, '.');
    tmp = tmp{1};
    name4{i} = tmp;
end
% find the intersection of the names
name_shared = name1(ismember(name1, name2));
name_shared = name_shared(ismember(name_shared, name3));
name_shared = name_shared(ismember(name_shared, name4));
% remove the samples that are not in all the folders
name_difference_1 = setdiff(name1, name_shared);
name_difference_2 = setdiff(name2, name_shared);
name_difference_3 = setdiff(name3, name_shared);
name_difference_4 = setdiff(name4, name_shared);
for i = 1:length(name_difference_1)
    delete(fullfile(an1, [name_difference_1{i}, '.cut.txt']));
    delete(fullfile(an1, [name_difference_1{i}, '.color.off']));
end
for i = 1:length(name_difference_2)
    delete(fullfile(an2, [name_difference_2{i}, '.cut.txt']));
    delete(fullfile(an2, [name_difference_2{i}, '.color.off']));
end
for i = 1:length(name_difference_3)
    delete(fullfile(an3, [name_difference_3{i}, '.cut.txt']));
    delete(fullfile(an3, [name_difference_3{i}, '.color.off']));
end
for i = 1:length(name_difference_4)
    delete(fullfile(volume_w_shaft_folder, [name_difference_4{i}, '.tif']));
end

