rootFolder = 'D:\Guoqiang_lab\VSOT_samples\VSOT_data_paired_with_codes\data\performance_testing_level_2_segmentation\large_dataset_3_annotator_L2345_w_stubby\gt_final_800\subtypes';
cut_index_foldr = fullfile(rootFolder, 'cut_index');
offFolder = fullfile(rootFolder, 'mushroom');
save_plot_folder = fullfile(rootFolder, 'cut_cycle_plot/mushroom');
if(~exist(save_plot_folder, 'dir'))
    mkdir(save_plot_folder);
end
files = dir(fullfile(cut_index_foldr, '*.cut.txt'));
filesname = {files.name};
for i = 1:length(filesname)
    namex = filesname{i};
    namex = strsplit(namex, '.');
    namex = namex{1};
    if(exist(fullfile(offFolder, [namex , '.off']), 'file'))

        [Pts,Tri] = read_off(fullfile(offFolder, [namex , '.off']));
        Tri = Tri';
        Pts = Pts';
        colormapx = [255, 0, 0, 255;255, 215, 0,255;0, 255, 0, 255];
        node_colorMap = repmat(colormapx(2,:),length(Pts),1);
        curve_vertex = readtable(fullfile(cut_index_foldr, [namex, '.cut.txt']));
        curve_vertex = table2array(curve_vertex);
        node_colorMap(curve_vertex,:) = repmat(colormapx(1,:),length(curve_vertex),1) ;
        fid = fopen([save_plot_folder,'/', namex,'.color.off'],'wt');
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

end