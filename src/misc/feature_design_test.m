spine_save_folder = '/work/boyu/EM_astrocyte/materials_for_paper_VSOT/test_samples_for_feature_design';
coordinate_output_folder = fullfile(spine_save_folder,'coor_result');
spine_head_neck_save_folder = fullfile(spine_save_folder,'head_neck_segmentation');
if(~exist(spine_head_neck_save_folder, 'dir'))
    mkdir(spine_head_neck_save_folder);
end

if(~exist(coordinate_output_folder, 'dir'))
    mkdir(coordinate_output_folder);
end

% test features 

curvature_scale = 30;
resx_s = 2;
resy_s = 2;
resz_s = 5;
% 11_232_7_3 has neck!!
comSeg.level_2_segmentation_main_speedup_output_cutVex(spine_save_folder, spine_head_neck_save_folder, coordinate_output_folder, curvature_scale, resx_s, resy_s, resz_s,1200,1);
listx = dir(fullfile(spine_head_neck_save_folder, '*.tif'));
names = {listx.name};
score_fin = [];
for i = 1:length(names)
    tmp = strsplit(names{i}, '.');
    
    if(exist(fullfile(spine_save_folder, [tmp{1}, '.off']), 'file'))
        spineHNROI = tiffreadVolume(fullfile(spine_head_neck_save_folder, names{i}));
        all_roi = bwlabeln(spineHNROI > 0);
        all_roi_idx = label2idx(all_roi);
        all_roi_idx = all_roi_idx(:);
        target_indicators = cellfun(@(x) sum(spineHNROI(x) == 3), all_roi_idx);
        target_mask = zeros(size(spineHNROI), 'uint8');
        target_mask(cell2mat(all_roi_idx(target_indicators > 0))) = 1;
        spineHNROI = spineHNROI.*target_mask;
        [lenx, leny, lenz] = size(spineHNROI);
        [~,headMeanRadiusx,neckLengthx, ~, neckMeanRadiusx, neckRadiusSTD, neckRadiusMedian, spheremap] = structQuant.genDendriteSpineScore_nocleft([], lenx, leny, lenz,spineHNROI, 16,16, 40);
        % headVolumex = headVolumex*16*16*40/1000/1000/1000;
        score_fin = [score_fin; [i, headMeanRadiusx,neckLengthx,neckMeanRadiusx, neckRadiusSTD, neckRadiusMedian, spheremap]];
    end
end
% nan_label = true(size(score_fin,1),1);
% for i = 2:7
% 
%     nan_label = nan_label&~isnan(score_fin(:,i))&~isinf(score_fin(:,i));
% 
% 
% end
score1 = zeros(size(score_fin, 1),1);
score2 = zeros(size(score_fin,1),1);
score3 = zeros(size(score_fin,1),1);
score4 = zeros(size(score_fin,1),1);
score5 = zeros(size(score_fin,1),1);
for i = 1:size(score_fin,1)
    score1(i) = 1/score_fin(i,2);
    if(score_fin(i,3) == 0)
        score2(i) = 0;
        score3(i) = 0;
        score4(i) = 1;
        score5(i) = 0;
    else
        score2(i) = score_fin(i,3)/(score_fin(i,4)^2);
        score3(i) = score_fin(i,5);
        score4(i) = exp(-score_fin(i,2)/score_fin(i,4));
        score5(i) = score_fin(i,3)/score_fin(i,2);
    end
end

k = 4; % number of clusters other than the stubby ones
% id_stubby = find(score2 == 0);
% id_non_stubby = find(score2 ~= 0);
% score1_remained = score1(score2 ~= 0);
% score1_remained = score1_remained(:);
% score2_remained = score2(score2 ~= 0);
% score2_remained = score2_remained(:);
score1 = normalize(score1, 'zscore');
score2 = normalize(score2, 'zscore');
score3 = normalize(score3, 'zscore');
score4 = normalize(score4, 'zscore');
score5 = normalize(score5, 'zscore');
% score5 = normalize(score_fin(:,2));
% score6 = normalize(score_fin(:,3));
scorematrix = [score1(:), score2(:), score3(:), score4(:),score5(:)];
[idx, C] = kmeans(scorematrix,k);
[coeff, score, latent, tsquared, explained] = pca(scorematrix);
figure;
gscatter(score(:, 1), score(:, 2), idx);
xlabel('PC1');
ylabel('PC2');
title('K-means Clustering of Spine Features');

% plot the box plot of head radius for each group
% boxplot(score_fin(:,2), idx);
Y = tsne(scorematrix);
figure;
gscatter(Y(:,1), Y(:,2), idx);
title('t-SNE of geometry features'); xlabel('Dimension 1'); ylabel('Dimension 2')

figure; boxplot(score_fin(:,2), idx);title('Head radius of each group');xlabel('Group ID'); ylabel('Radius(nm)')
figure; boxplot(score_fin(:,5), idx);title('Neck length std of each group')
figure; boxplot(score_fin(:,3)/1000, idx);title('Neck length of each group');xlabel('Group ID'); ylabel('Length(\mum)')
figure; boxplot(score_fin(:,2)./score_fin(:,4), idx);title('Ratio between head-neck radius');xlabel('Group ID'); ylabel('Ratio')

% add the stubby ones back
idx_group = label2idx(idx);
% idx_all = zeros(size(score1));
% idx_all(id_stubby) = 5;
% for i = 1:k
%     idx_all(id_non_stubby(idx_group{i})) = i;
% end
% idx_group = label2idx(idx_all);

% figure;gscatter(score1(:), score2(:), idx_all)
% xlabel('head resistance'); ylabel('neck resistance')
for i = 1:k
    name_group1 = names(score_fin(idx_group{i},1));
    folderx = fullfile(spine_save_folder, ['group_', num2str(i)]);
    if ~exist(folderx, 'dir')
        mkdir(folderx);
    else
        rmdir(folderx, 's');
        mkdir(folderx);
    end
    % % Randomly select 36 samples
    num_samples = min(36, length(name_group1));
    rand_indices = randperm(length(name_group1), num_samples);
    selected_samples = name_group1(rand_indices);
    
    % % Plot the selected samples in a 6x6 grid
    fig = figure;
    for j = 1:num_samples
        tmp = strsplit(selected_samples{j}, '.');
        [Pts, Tri] = read_off(fullfile(spine_save_folder, [tmp{1}, '.off']));
        Tri = Tri';
        Pts = Pts';
        subplot(6, 6, j);
        trisurf(Tri, Pts(:, 1), Pts(:, 2), Pts(:, 3), 'FaceColor', 'interp', 'EdgeColor', 'none');

    end
    title(['group_',num2str(i)])
    % saveas(fig, fullfile(folderx, 'random_samples_grid.png'));
    parfor j = 1:length(name_group1)
        tmp = strsplit(name_group1{j}, '.');
        [Pts, Tri] = read_off(fullfile(spine_save_folder, [tmp{1}, '.off']));
        Tri = Tri';
        Pts = Pts';
        fig = figure('visible', 'off');
        trisurf(Tri, Pts(:, 1), Pts(:, 2), Pts(:, 3), 'FaceColor', 'interp', 'EdgeColor', 'none');
        saveas(fig, fullfile(folderx, [name_group1{j}, '.png']));
    end

end
