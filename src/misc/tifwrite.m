function tifwrite(outLabel, ImName)
ImName = char(ImName);
if size(outLabel,4)==1
    imwrite(outLabel(:,:,1),[ImName,'.tif']);
    for i = 2:size(outLabel,3)
        imwrite(outLabel(:,:,i),[ImName,'.tif'],'WriteMode','append');
    end
else
    imwrite(squeeze(outLabel(:,:,1,:)),[ImName,'.tif']);
    for i = 2:size(outLabel,4)
        %disp(i);
        imwrite(squeeze(outLabel(:,:,i,:)),[ImName,'.tif'],'WriteMode','append');
    end
end
end