function [surfaceCellHead, nodeCellHead,surfaceCellNeck,nodeCellNeck] = genCutSurface_speedup_head_neck(cutVexcell, mF2, headp, neckp, Vnieuwx,Vnieuw, test_region2, allTetrax, allVertex)
    %%
    
    xxshift0 = zeros(3,3,3);
    yyshift0 = zeros(3,3,3);
    zzshift0 = zeros(3,3,3);
    for i = -1:1
        for j = -1:1
            for k = -1:1
                xxshift0(i+2, j+2, k+2) = i;
                yyshift0(i+2, j+2, k+2) = j;
                zzshift0(i+2, j+2, k+2) = k;
            end
        end
    end
    
    
    [lenx, leny, lenz] =size(test_region2);
    face_p = find(mF2(:,1) == headp | mF2(:,2) == headp|mF2(:,3) == headp);
    face_n = find(mF2(:,1) == neckp | mF2(:,2) == neckp|mF2(:,3) == neckp);
    surfaceCellHead = cell(length(cutVexcell),1);
    nodeCellHead = cell(length(cutVexcell),1); 
    surfaceCellNeck = cell(length(cutVexcell),1);
    nodeCellNeck = cell(length(cutVexcell),1);
    for i = 1:length(cutVexcell)
        try
        % i = 75
    %         disp(i)
            cutVex = cutVexcell{i};
            if(length(cutVex) < 10)
                cutVexcell{i} = nan;
            elseif(any(cutVex == headp) || any(cutVex == neckp))
                cutVexcell{i} = nan;
            else
                bbNode = Vnieuw(cutVex,:);
                xlim_low = max(floor(min(bbNode(:,1))) - 3,1);
                xlim_upp = min(ceil(max(bbNode(:,1))) + 3, lenx);
                ylim_low = max(floor(min(bbNode(:,2))) - 3, 1);
                ylim_upp = min(ceil(max(bbNode(:,2))) + 3, leny);
                zlim_low = max(floor(min(bbNode(:,3))) - 3,1);
                zlim_upp = min(ceil(max(bbNode(:,3))) + 3,lenz);
                bbNode(:,1) = bbNode(:,1) - xlim_low + 1;
                bbNode(:,2) = bbNode(:,2) - ylim_low + 1;
                bbNode(:,3) = bbNode(:,3) - zlim_low + 1;
                id_allVertex = [1:size(allVertex,1)];
                id_allVertex = id_allVertex(allVertex(:,1) >= xlim_low & allVertex(:,1) <= xlim_upp & ...
                    allVertex(:,2) >= ylim_low & allVertex(:,2) <= ylim_upp & allVertex(:,3) >= zlim_low & allVertex(:,3) <= zlim_upp);
                id0XYZ = allVertex(allVertex(:,1) >= xlim_low & allVertex(:,1) <= xlim_upp & ...
                    allVertex(:,2) >= ylim_low & allVertex(:,2) <= ylim_upp & allVertex(:,3) >= zlim_low & allVertex(:,3) <= zlim_upp, :);
                id0XYZ(:,1) = id0XYZ(:,1) - xlim_low + 1;
                id0XYZ(:,2) = id0XYZ(:,2) - ylim_low + 1;
                id0XYZ(:,3) = id0XYZ(:,3) - zlim_low + 1;            
                %% find a smaller bounding box based on pca and then map back to the original scale
                [coeff,score_w,latent] = pca(bbNode);            
    %             bbNode2 =  bbNode - mean(bbNode,1);
                bbx_rotated = (id0XYZ - mean(bbNode,1))*coeff;
                id_allVertex(bbx_rotated(:,1)> max(ceil(score_w(:,1))) +1 | bbx_rotated(:,2)> max(ceil(score_w(:,2)))+1 | bbx_rotated(:,3)> max(ceil(score_w(:,3)))+1|...
                bbx_rotated(:,1)< min(floor(score_w(:,1)))-1 | bbx_rotated(:,2)< min(floor(score_w(:,2)))-1 | bbx_rotated(:,3)< min(floor(score_w(:,3)))-1) = [];
                id0XYZ(bbx_rotated(:,1)> max(ceil(score_w(:,1))) +1 | bbx_rotated(:,2)> max(ceil(score_w(:,2)))+1 | bbx_rotated(:,3)> max(ceil(score_w(:,3)))+1|...
                bbx_rotated(:,1)< min(floor(score_w(:,1)))-1 | bbx_rotated(:,2)< min(floor(score_w(:,2)))-1 | bbx_rotated(:,3)< min(floor(score_w(:,3)))-1,:) = [];
                smallTetrax = allTetrax(ismember(allTetrax(:,1), id_allVertex)|ismember(allTetrax(:,2), id_allVertex)|ismember(allTetrax(:,3), id_allVertex)|ismember(allTetrax(:,4), id_allVertex),:);
                id_allVertex = unique(smallTetrax(:));
                id0XYZ = allVertex(id_allVertex,:);
                id0XYZ(:,1) = id0XYZ(:,1) - xlim_low + 1;
                id0XYZ(:,2) = id0XYZ(:,2) - ylim_low + 1;
                id0XYZ(:,3) = id0XYZ(:,3) - zlim_low + 1;  
                smallFacesX = [smallTetrax(:, [1,2,3]), [1:size(smallTetrax,1)]'; smallTetrax(:,[1,3,4]), [1:size(smallTetrax,1)]';smallTetrax(:, [1,2,4]), [1:size(smallTetrax,1)]';smallTetrax(:, [2,3,4]), [1:size(smallTetrax,1)]'];
                id_allVertex_mapping = zeros(max(id_allVertex),1);
                id_allVertex_mapping(id_allVertex) = 1:length(id_allVertex);
                smallFacesX(:,1:3) = id_allVertex_mapping(smallFacesX(:,1:3));
                smallFacesX(:,1:3) = sort(smallFacesX(:,1:3), 2);     
                smallFacesX = sortrows(smallFacesX, [1:3]);
                [cutMF, cutvertices] = comSeg.findMinCut_v5_speedup(id0XYZ, smallFacesX,bbNode,lenx, leny, lenz); % min-cut
                cutvertices(:,1) = (cutvertices(:,1) + xlim_low - 1);
                cutvertices(:,2) = (cutvertices(:,2) + ylim_low - 1);
                cutvertices(:,3) = (cutvertices(:,3) + zlim_low - 1);
                [headHalf_surface] = comSeg.splitAlongCycle_speedup(mF2, cutVex,face_p); % graph connected components
                [neckHalf_surface] = comSeg.splitAlongCycle_speedup(mF2, cutVex,face_n); % graph connected components
                % combine the vertices from the half of the surfce and the cut
                % surface
                [combinedSurfaceHead, verticesFinHead] = comSeg.combineSurface(cutMF, cutvertices, headHalf_surface, Vnieuw, Vnieuw(cutVex,:));
                [combinedSurfaceNeck, verticesFinNeck] = comSeg.combineSurface(cutMF, cutvertices, neckHalf_surface, Vnieuw, Vnieuw(cutVex,:));
                if(size(combinedSurfaceHead, 1) > 2000)
                    [combinedSurfaceHead,verticesFinHead] = reducepatch(combinedSurfaceHead,verticesFinHead,2000);
    
                end
                surfaceCellHead{i} = combinedSurfaceHead;
                nodeCellHead{i} = verticesFinHead;
                if(size(combinedSurfaceNeck, 1) > 2000)
                    [combinedSurfaceNeck,verticesFinNeck] = reducepatch(combinedSurfaceNeck,verticesFinNeck,2000);
    
                end
                surfaceCellNeck{i} = combinedSurfaceNeck;
                nodeCellNeck{i} = verticesFinNeck;
    
            end
        catch ME
            surfaceCellHead{i} = [];
            nodeCellHead{i} = [];
            surfaceCellNeck{i} = [];
            nodeCellNeck{i} = [];
    %              warning('problem happens at %s', num2str(i))
            continue;
        end
    end
    
    end