function map_volume_label_2_surface(label, Pts, Tri, output, output_mat_num)
[Pts,Tri] = read_off('/work/boyu/EM_astrocyte/test_segmentation_samples/annotation/11_324_12_4.off');
Tri  =Tri';
Pts = Pts';

labelx = tiffreadVolume('/work/boyu/EM_astrocyte/test_segmentation_samples/dendrite_spine_ours_result/11_324_12_4_hn_result.tif');
labelxidx = label2idx(labelx);
[lenx, leny, lenz] = size(labelx);
[headCoorx,headCoory, headCoorz]  = ind2sub([lenx, leny, lenz], labelxidx{2});
[neckCoorx,neckCoory, neckCoorz]  = ind2sub([lenx, leny, lenz], labelxidx{3});
Tri_center = [mean([Pts(Tri(:,1), 1),Pts(Tri(:,2), 1),Pts(Tri(:,3), 1)],2)/2, ...
    mean([Pts(Tri(:,1), 2),Pts(Tri(:,2), 2),Pts(Tri(:,3), 2)],2)/2,...
    mean([Pts(Tri(:,1), 3),Pts(Tri(:,2), 3),Pts(Tri(:,3), 3)],2)/5];

% labelID = output(2,:);
% labelID = labelID{1};
% C = strsplit(labelID,' ');
% labelIDx = str2double(C);
% figure;trisurf(Tri,Pts(:,1),Pts(:,2),Pts(:,3),'Facecolor','red','FaceAlpha',0.1);


node_colorMap = zeros(length(Pts),4);
node_colorMap2 = zeros(length(Pts),1);
colormapx = [255, 0, 0, 255;0, 0, 255, 255;0, 255, 0, 255];
for i = 1:size(Pts,1)
    distHead = sqrt((Pts(i,1)/2 - headCoorx).^2 + (Pts(i,2)/2 - headCoory).^2 + (Pts(i,3)/5 - headCoorz).^2);
    distNeck = sqrt((Pts(i,1)/2 - neckCoorx).^2 + (Pts(i,2)/2 - neckCoory).^2 + (Pts(i,3)/5 - neckCoorz).^2);
    if(min(distHead) > min(distNeck))
        node_colorMap(i,:) = colormapx(1,:);
    else
        node_colorMap(i,:) = colormapx(2,:);
    end

end

% for i = 1:length(output_mat_num)
%     for j = 1:3
%         if(node_colorMap(Tri(i,j)) >0 && node_colorMap(Tri(i,j)) ~= output_mat_num(i))
%             node_colorMap2(Tri(i,j)) = 1;
%         elseif( node_colorMap(Tri(i,j)) == 0)
%             node_colorMap(Tri(i,j),:) = colormapx(output_mat_num(i) + 1,:);
%         end
%     end
% end
figure; trisurf(Tri, Pts(:,1), Pts(:,2), Pts(:,3))
figure; trisurf(Tri, Pts(:,1), Pts(:,2), Pts(:,3), 'FaceVertexCData',node_colorMap2,'Facecolor','interp') 
% 
% outputColoder = repmat([192,192,192,255], size(Pts,1), 1);
% outputColoder(node_colorMap2 == 1,:) = repmat([255, 0, 0, 255], sum(node_colorMap2 == 1),1);

fid = fopen('/work/boyu/EM_astrocyte/test_segmentation_samples/dendrite_spine_ours_result/13_212_2_5_hn_result_with2Label.off','wt');
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