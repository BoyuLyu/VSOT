function [face_graph0, uniqueXYZ] = genSurfaceTetraGraph_new(test_region2)

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


allPointxyz = [0,0,0;
    1,0,0;
    0,1,0;
    1,1,0;
    0,0,1;
    1,0,1;
    0,1,1;
    1,1,1;];
allPointx = allPointxyz(:,1);
allPointy = allPointxyz(:,2);
allPointz = allPointxyz(:,3);

allPointxSS = allPointx(:)' + xxshift(:);
allPointySS = allPointy(:)' + yyshift(:);
allPointzSS = allPointz(:)' + zzshift(:);
allPointXYZ = [allPointxSS(:), allPointySS(:), allPointzSS(:)];

idRM = find(allPointXYZ(:,1) < 0 | allPointXYZ(:,1) > 1 | allPointXYZ(:,2) < 0| allPointXYZ(:,2) > 1 | allPointXYZ(:,3) < 0| allPointXYZ(:,3) > 1);
groupx = unique(floor((idRM - 1)/4) + 1);
idRM_groupID = zeros(4*length(groupx), 1);
for i = 1:length(groupx)
    idRM_groupID((4*(i-1)+1) : 4*i) = [((groupx(i) - 1)*4 + 1):(4*groupx(i))];
end
allPointXYZ(idRM_groupID,:) = [];

% [allPointx,allPointy,allPointz] = ind2sub([lenx, leny, lenz],test_region2ID);


% in each tetrahedron further split into five small tetrhedrons
newNodesXYZ = zeros(size(allPointXYZ,1)*5, 3);
for i = 1:(size(allPointXYZ,1)/4)
    tmpxyz = allPointXYZ((i-1)*4+1:i*4,:);
    pointAB = tmpxyz(floor(tmpxyz(:,1)) == tmpxyz(:,1) & floor(tmpxyz(:,2)) == tmpxyz(:,2) & floor(tmpxyz(:,3)) == tmpxyz(:,3),:);
    pointA = pointAB(1,:);
    pointB = pointAB(2,:);
    pointC = tmpxyz(floor(tmpxyz(:,1)) ~= tmpxyz(:,1) & floor(tmpxyz(:,2)) ~= tmpxyz(:,2) & floor(tmpxyz(:,3)) ~= tmpxyz(:,3),:);
    pointD = setdiff(tmpxyz, [pointAB; pointC], 'rows');
    pointCenter = (pointC + (pointA + pointB)/2)/2; % the center point of ABC
    pointLeft = (pointC - pointA)*2/3 + pointA;
    pointRight = (pointC - pointB)*2/3 + pointB;
    t1 = [pointA; pointB; pointD; pointCenter];
    t2 = [pointA; pointCenter; pointLeft; pointD];
    t3 = [pointB; pointCenter; pointRight; pointD];
    t4 = [pointD; pointC; pointCenter; pointLeft];
    t5 = [pointD; pointC; pointCenter; pointRight];
    newNodesXYZ((i-1)*20+1 : i*20,:) = [t1;t2;t3;t4;t5];
end

[uniqueXYZ, ia, ic] = unique(newNodesXYZ, 'rows');
tera_x = reshape(ic, 4,[]);
tera_x = tera_x';
tera_x = sort(tera_x, 2);
tera_x = unique(tera_x, 'rows');
% tera_x = mappingx(tera_x);
facesx = [tera_x(:, [1,2,3]), [1:size(tera_x,1)]'; tera_x(:,[1,3,4]), [1:size(tera_x,1)]';tera_x(:, [1,2,4]), [1:size(tera_x,1)]';tera_x(:, [2,3,4]), [1:size(tera_x,1)]'];
facesx(:,1:3) = sort(facesx(:,1:3), 2);
facesx = sortrows(facesx, [1:3]);
nodeGraph = [];
for i = 1:3
    for j = (i+1):4
        nodeGraph = [nodeGraph; [tera_x(:,i), tera_x(:,j)]];
        
    end
end
nodeGraph = sort(nodeGraph, 2);
nodeGraph = unique(nodeGraph, 'rows');
G = graph(nodeGraph(:,1), nodeGraph(:,2));
figure;
p = plot(G, 'LineWidth', 2);
p.XData = uniqueXYZ(:,1);
p.YData = uniqueXYZ(:,2);
p.ZData = uniqueXYZ(:,3);

[~, ia_face,ic_face] = unique(facesx(:, 1:3), 'rows');
countsx = accumarray(ic_face, 1);
id_face_single = ia_face(countsx == 1);
facesx2 = facesx;
facesx2(id_face_single,:) = [];
facesx2(:,1:3) = sort(facesx2(:,1:3), 2);
facesx2 = sortrows(facesx2, [1:3]);
tera_graph = [facesx2(1:2:end-1,4), facesx2(2:2:end,4)];

G2 = graph(tera_graph(:,1), tera_graph(:,2));
tetra_center = zeros(size(tera_x,1),3);
for i = 1:size(tera_x,1)
    tetra_center(i,:) = mean( uniqueXYZ(tera_x(i,:),:),1);
    
    
    
end
figure;
q = plot(G2, 'LineWidth',2);
q.XData = tetra_center(:,1);
q.YData = tetra_center(:,2);
q.ZData = tetra_center(:,3);



hold on;











end