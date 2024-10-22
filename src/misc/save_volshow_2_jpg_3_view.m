function save_volshow_2_jpg(folder, filename)
    tmp_mask = tiffreadVolume(fullfile(folder, filename));
    fig = figure('Visible', 'off');
    % Create the first subplot
    subplot(1, 3, 1);
    imagesc(max(tmp_mask, [],3));
    title('front view  (Z)');
    subplot(1, 3, 2);
    % Create the second subplot
    imagesc(squeeze(max(tmp_mask, [],2)));
    title('side view (Y)');
     subplot(1, 3, 3);           
    % Create the third subplot
    imagesc(squeeze(max(tmp_mask, [],1)));
    title('top view (X)');
    
    % Save the figure to a file
    saveas(fig, fullfile(folder, [filename, '.jpg']));
    
    % Close the figure
    close(fig);
end
