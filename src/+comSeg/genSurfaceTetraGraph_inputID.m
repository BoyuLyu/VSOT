function [face_graph0, uniqueXYZ,facesx] = genSurfaceTetraGraph_inputID(test_region2ID,all_cc_points,lenx, leny, lenz)


xxshift = [];
yyshift = [];
zzshift = [];
for i = -1:1
    for j = -1:1
        for k = -1:1
            if((i==0 && j==0 && k~=0))
               for m = [-1,1]
                   for n = [-1,1]
                      xxshift = [xxshift, [0, 0, m*0.5, 0,     0, 0, m*0.5, m*0.5]]; 
                      yyshift = [yyshift, [0, 0, n*0.5, n*0.5, 0, 0, n*0.5, 0]];
                      zzshift = [zzshift, [0, k, 0.5*k, 0.5*k ,0, k, 0.5*k, 0.5*k]];
                   end
               end
            elseif(i==0 && j~=0 && k == 0)
               for m = [-1,1]
                   for n = [-1,1]
                      xxshift = [xxshift, [0, 0, m*0.5, 0,     0, 0, m*0.5, m*0.5]]; 
                      yyshift = [yyshift, [0, j, 0.5*j, 0.5*j, 0, j, 0.5*j, 0.5*j]];
                      zzshift = [zzshift, [0, 0, n*0.5, n*0.5 ,0, 0, n*0.5, 0]];
                   end
               end                
            elseif(i~=0 && j == 0 && k == 0)
               for m = [-1,1]
                   for n = [-1,1]
                      xxshift = [xxshift, [0, i, 0.5*i, 0.5*i, 0, i, 0.5*i, 0.5*i]]; 
                      yyshift = [yyshift, [0, 0, m*0.5, 0,     0, 0, m*0.5, m*0.5]];
                      zzshift = [zzshift, [0, 0, n*0.5, n*0.5 ,0, 0, n*0.5, 0]];
                   end
               end            
            end
        end
    end
end


[allPointx,allPointy,allPointz] = ind2sub([lenx, leny, lenz],test_region2ID);

% curSurface = [];
% vertices = [];
% allPointx = id0XYZ(:,1);
% allPointy = id0XYZ(:,2);
% allPointz = id0XYZ(:,3);
id0XYZ = [allPointx(:),allPointy(:),allPointz(:)];
allPointxSS = allPointx(:)' + xxshift(:);
allPointySS = allPointy(:)' + yyshift(:);
allPointzSS = allPointz(:)' + zzshift(:);
allPointXYZ = [allPointxSS(:), allPointySS(:), allPointzSS(:)];
idRM = find(allPointXYZ(:,1) < 1 | allPointXYZ(:,1) > lenx | allPointXYZ(:,2) < 1| allPointXYZ(:,2) > leny | allPointXYZ(:,3) < 1| allPointXYZ(:,3) > lenz);
groupx = unique(floor((idRM - 1)/4) + 1);
idRM_groupID = zeros(4*length(groupx), 1);
for i = 1:length(groupx)
    idRM_groupID((4*(i-1)+1) : 4*i) = [((groupx(i) - 1)*4 + 1):(4*groupx(i))];
end
allPointXYZ(idRM_groupID,:) = [];
[uniqueXYZ, ia, ic] = unique(allPointXYZ, 'rows');
newlyadded_interger = uniqueXYZ(uniqueXYZ(:,1) == floor(uniqueXYZ(:,1)) & uniqueXYZ(:,2) == floor(uniqueXYZ(:,2)) & uniqueXYZ(:,3) == floor(uniqueXYZ(:,3)), :);
newlyadded_interger = newlyadded_interger(~ismember(newlyadded_interger, id0XYZ, 'rows'),:);
newlyadded_interger_neix = newlyadded_interger(:,1)' + xxshift(:);
newlyadded_interger_neiy = newlyadded_interger(:,2)' + yyshift(:);    
newlyadded_interger_neiz = newlyadded_interger(:,3)' + zzshift(:);
newlyadded_interger_neiXYZ = [newlyadded_interger_neix(:), newlyadded_interger_neiy(:), newlyadded_interger_neiz(:)];
newlyadded_interger_neiXYZ = unique(newlyadded_interger_neiXYZ, 'rows');

newlyadded_interger_neiXYZ = newlyadded_interger_neiXYZ(~ismember(newlyadded_interger_neiXYZ, [all_cc_points(:,1), all_cc_points(:,2), all_cc_points(:,3)], 'rows'), :);

int_id_rm = find(ismember(uniqueXYZ,newlyadded_interger_neiXYZ, 'rows' ));
mappingx = zeros(size(uniqueXYZ,1), 1);
index0 = [1:size(uniqueXYZ,1)];
index0 = setdiff(index0, int_id_rm);
mappingx(index0) = 1:length(index0);
uniqueXYZ(int_id_rm,:) = [];
tera_x = reshape(ic, 4,[]);
tera_x = tera_x';
tera_x = sort(tera_x, 2);
tera_x = unique(tera_x, 'rows');
tera_x = mappingx(tera_x);
tera_x(tera_x(:,1) == 0| tera_x(:,2) == 0| tera_x(:,3) == 0| tera_x(:,4) ==0,:) = [];
facesx = [tera_x(:, [1,2,3]), [1:size(tera_x,1)]'; tera_x(:,[1,3,4]), [1:size(tera_x,1)]';tera_x(:, [1,2,4]), [1:size(tera_x,1)]';tera_x(:, [2,3,4]), [1:size(tera_x,1)]'];
facesx(:,1:3) = sort(facesx(:,1:3), 2);
facesx = sortrows(facesx, [1:3]);
[~, ia_face,ic_face] = unique(facesx(:, 1:3), 'rows');
countsx = accumarray(ic_face, 1);
id_single_tri = ia_face(countsx == 1);
face_graph0 = facesx(id_single_tri, 1:3);% the outer surface of the tetrahedron graph
% [face_graph0, uniqueXYZ] = refineEdge(face_graph0, uniqueXYZ, xxshift, yyshift, zzshift ,lenx, leny, lenz);

end