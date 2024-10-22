function examine_plot_w_speedup(before_speedup_folder, after_speedup_folder, vis_save_folder)

if(~exist(vis_save_folder, 'dir'))
    mkdir(vis_save_folder);
end
% Load the files in the before_speedup_foldet
before_tif_files = dir(fullfile(before_speedup_folder, '*.tif'));
after_tif_files = dir(fullfile(after_speedup_folder, '*.tif'));

for i = 1:length(before_tif_files)
        tmp_mask1 = tiffreadVolume(fullfile(before_speedup_folder, before_tif_files(i).name));
        fig = figure('Visible', 'off');
        % Create the first subplot
        subplot(2, 3, 1);
        imagesc(max(tmp_mask1, [], 3));
        title('front view (Z)');

        % Create the second subplot
        subplot(2, 3, 2);
        imagesc(squeeze(max(tmp_mask1, [], 2)));
        title('side view (Y)');

        % Create the third subplot
        subplot(2, 3, 3);
        imagesc(squeeze(max(tmp_mask1, [], 1)));
        title('top view (X)');

        tmp_mask2 = tiffreadVolume(fullfile(after_speedup_folder, before_tif_files(i).name));

        % Create the fourth subplot
        subplot(2, 3, 4);
        imagesc(max(tmp_mask2, [], 3));
        title('front view (Z)');

        % Create the fifth subplot
        subplot(2, 3, 5);
        imagesc(squeeze(max(tmp_mask2, [], 2)));
        title('side view (Y)');

        % Create the sixth subplot
        subplot(2, 3, 6);
        imagesc(squeeze(max(tmp_mask2, [], 1)));
        title('top view (X)');
        % Save the figure to a file
        saveas(fig, fullfile(vis_save_folder, [before_tif_files(i).name, '.jpg']));
        
        % Close the figure
        close(fig);
end


end