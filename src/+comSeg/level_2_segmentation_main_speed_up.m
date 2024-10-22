function level_2_segmentation_main_speed_up(spine_save_folder, spine_head_neck_save_folder, outputFolder, curvature_scale, resx, resy, resz,minpoints,smooth_kernel)
    % curvature_scale: the scale used to calculate the curvature score (default 30)
    % resx, resy, resz can be set as the greatest common divisor of the original resolution (16x16x40 -> 2x2x5)
    % resx = 2, resy = 2, resz = 5;
    % rather than computing the two scores and then combined in the end we select the candidate points by selecting from the surface score
    % then calculate the geometry score for each point and then combine the two scores
    files_x = dir(fullfile(spine_save_folder, '*.tif'));
    if(~exist(spine_head_neck_save_folder, "dir"))
        mkdir(spine_head_neck_save_folder)
    end
    tic;
    coor_head_cell = cell(length(files_x), 1);
    coor_neck_cell = cell(length(files_x), 1);
    parfor k = 1:length(files_x)
        disp(k)
        try
            spineROI = tiffreadVolume(fullfile(spine_save_folder,[num2str(k),'.tif']));
            shaft_mask = spineROI == 1;
            test_region = spineROI == 2;
            test_region_roi = bwlabeln(test_region);
            test_region_roi_idx = label2idx(test_region_roi);
            len_region = cellfun(@length, test_region_roi_idx);
            test_region = test_region_roi == (find(len_region == max(len_region),1));
            test_region0 = test_region;
            % test_region0 = test_region_roi == (find(len_region == max(len_region),1));
            % head_volume = zeros(size(test_region));
            % neckx = zeros(size(test_region));
            [lenxtmp, lenytmp, lenztmp] = size(test_region);
            % curID1 = find(test_region(:) == 1);
            % [curID1x, curID1y, curID1z] = ind2sub([lenxtmp, lenytmp, lenztmp], curID1);
            % test_regionID = find(test_region(:) == 1);  
            %% generate the surface using the original scale
            se = strel('sphere', smooth_kernel);
            test_region = imdilate(test_region, se);
            test_region = imclose(test_region, se);
    %%=========================================================================================================
            % generate the neighbor matrix for the new computational unit
            [xxshift, yyshift, zzshift] = comSeg.genShiftMatrixTetra;
            % generate the surface as well as the tetrahedron based on our designed computational unit
            % the vertices of the surfaces are located at the grid points/ center of grid points in the original volume
            [mF2, allVertex, allTetrax] = comSeg.genSurfaceTetraGraph_v3(test_region, xxshift, yyshift, zzshift);
            Vnieuw = allVertex; % allVertex contains all the vertices in the volume
            %make sure the vertices matrix are all on the surface of the volume
            vertex_unique = unique(mF2(:));
            mappingx = zeros(max(vertex_unique),1);
            mappingx(vertex_unique(:)) = [1:length(vertex_unique(:))];
            mF2 = mappingx(mF2);
            Vnieuw = Vnieuw(vertex_unique(:),:);
            Vnieuwx = [resx*Vnieuw(:,1),resy*Vnieuw(:,2),resz*Vnieuw(:,3)];% convert the vertices to the real coordinate system
            mF2 = comSeg.checkWindingOrder(mF2); % make sure the winding order is correct for faster computation in the smoothing step
    %           figure;trisurf(mF2,Vnieuw(:,1),Vnieuw(:,2),Vnieuw(:,3),'Facecolor','red','FaceAlpha',0.1);
            % #60 & #284
            [Vnieuw2_out,mF2_out] = taubin_smooth_tiny_mesh_mex(Vnieuwx, mF2, 0.8, 0.53, 40); % mex file compiled from the c++ code for speedup
            mappingx = zeros(size(Vnieuwx,1),1);
            mF2_t = mF2';
            [uniquex,ia,ic]= unique(mF2_t(:),'stable');
            mappingx(uniquex) = [1:length(mappingx)];
            Vnieuw2 = Vnieuw2_out(mappingx,:);
            % idSelected = 1:size(Vnieuwx,1);
            % figure;trisurf(mF2,Vnieuw2(:,1),Vnieuw2(:,2),Vnieuw2(:,3),'Facecolor','red','FaceAlpha',0.1);
            [mF3,Vnieuw3] = reducepatch(mF2,Vnieuwx,min(size(mF2,1),4000));% reduce the number of the vertices for calculation for speedup
            [~, idSelected] = ismember(Vnieuw3, Vnieuwx,'rows'); % start using the smaller set of the vertex to generate the two scores
            %generate the curvature score
%                     
    %%=========================================================================================================
            % Search for the tip of the head and tip of the neck
            % the tip of the neck is obtained by searching for the contact site with the dendrite shaft
            % the tip of the head is obtained by searching for the furthest point from the tip of the neck
    
            se = strel('sphere', 1);
            contactSite = test_region.* imdilate(shaft_mask, se);
            contactID = find(contactSite(:) == 1);    
            [contactIDx, contactIDy, contactIDz] = ind2sub([lenxtmp, lenytmp, lenztmp], contactID);
            nearest_vertex = round(Vnieuw);
            nearest_vertexID = sub2ind([lenxtmp, lenytmp, lenztmp], nearest_vertex(:,1), nearest_vertex(:,2), nearest_vertex(:,3));
            % idrm1 = find(test_region(nearest_vertexID) == 0);
            geoDist2Starting = comSeg.getGeoDist_server(test_region, contactID,resx, resy, resz);
            % distCircleCenter = zeros(size(Vnieuw,1),1);
            distVexStartingPoint = zeros(size(Vnieuw,1),1);
            for j = 1:size(nearest_vertex,1)
                distVexStartingPoint(j) = geoDist2Starting(nearest_vertex(j,1), nearest_vertex(j,2), nearest_vertex(j,3));
            end
            distVexStartingPoint(distVexStartingPoint == 1.0000e+10)= nan;
            neckpoint = find(distVexStartingPoint == min(distVexStartingPoint),1); % the vertex id that close to the starting point
            headpoint = find(distVexStartingPoint == max(distVexStartingPoint), 1); % the vertex id of the furthest point 
            % assign certain points as invalid
            distVexStartingPoint(distVexStartingPoint == 0) = 1;
            ccHead = find(distVexStartingPoint > quantile(distVexStartingPoint,0.98));
            ccHead_all_nei = mF2(ismember(mF2(:,1), ccHead)| ismember(mF2(:,2), ccHead)|ismember(mF2(:,3), ccHead),:);
            ccHead = unique(ccHead_all_nei(:));
            ccNeck = find(distVexStartingPoint <= quantile(distVexStartingPoint,0.05)); % all the faces related to these two will be categorized as real head and neck    
            ccNeck_all_nei = mF2(ismember(mF2(:,1), ccNeck)| ismember(mF2(:,2), ccNeck)|ismember(mF2(:,3), ccNeck),:);
            ccNeck = unique(ccNeck_all_nei(:));

    %%=========================================================================================================
    % Calculate the curvature score using the adaptive scale
            idSelected = setdiff(idSelected, [ccHead(:);ccNeck(:)]);
            [newCurvature] = comSeg.genCurvatureFace4Scale_speedup(mF2, Vnieuw2, curvature_scale, idSelected); % 48 cores 2.16 s
            % figure;trisurf(mF2,Vnieuw2(:,1),Vnieuw2(:,2),Vnieuw2(:,3),'FaceVertexCData',newCurvature,'Facecolor','interp');   

    %%=========================================================================================================
            % select the candidate points based on the surface score and then calculate the geometry score
            newCurvature(ccHead) = nan;
            newCurvature(ccNeck) = nan;
            [~,p] = sort(newCurvature);
            % newCurvature_order = newCurvature;
            % newCurvature_order(p) = 1:length(p);
            
            idSelected = p(1:min(minpoints,sum(~isnan(newCurvature))));
            % newCurvature_order = newCurvature_order.*0;
            % newCurvature_order(idSelected) = 1:length(idSelected);
            
    %% start the process calculating the geometry score by obtain the optimal segmentation of the head and neck for each point
            % find the shortest separating cycle for each point
            cutVex2_speedup = comSeg.searchSSpathCutGraph_v2_speedup(Vnieuwx, idSelected, mF2,headpoint,neckpoint, ccHead, ccNeck);
            % for the cut paths that have repetitive points, remove
            % the points between the repetitive ones
            cutVexcellUnique2 = comSeg.cleanReptitve_v2(cutVex2_speedup);
    %%=========================================================================================================
            % find the optimal cut for each point based on the cycle using the computational unit we designed
            [surfaceCellHead, nodeCellHead] = comSeg.genCutSurface_speedup(cutVexcellUnique2(idSelected), mF2, headpoint,neckpoint,Vnieuwx,Vnieuw, test_region, allTetrax, allVertex);
            sphereMap = zeros(length(idSelected),1) + nan;
            for j = 1:length(sphereMap)
                if(~isempty(surfaceCellHead{j}))
                    verticesFin = nodeCellHead{j};
                    combinedSurface = surfaceCellHead{j}; % this contains the surface of the head part after segmentation
                    vectorMF3_1 = verticesFin(combinedSurface(:,3),:) - verticesFin(combinedSurface(:,1),:);
                    vectorMF3_2 = verticesFin(combinedSurface(:,3),:) - verticesFin(combinedSurface(:,2),:);
                    ss_tmp = cross(vectorMF3_1, vectorMF3_2,2);
                    ss = 1/2*sum(sqrt(ss_tmp(:,1).^2 + ss_tmp(:,2).^2 + ss_tmp(:,3).^2));
                    volume = abs(comSeg.calVolfromMesh(combinedSurface, verticesFin));                                               
                    r1 = (volume/(4*pi/3))^(1/3);
                    r2 = sqrt(ss/(4*pi));
                    sphereMap(j) = abs(r1 - r2)/r2;
                end
            end
            sphereMap(sphereMap(:) == 0) = nan;
    %%=========================================================================================================
    % Calculate the centerline width score for the cylinder map in the 
            [pathidxyz, width_eachp] = comSeg.calWidthfromVol(test_region, Vnieuw(neckpoint,:), Vnieuw(headpoint,:), resx, resy, resz); % width is from the tip of the shaft to the tip of the head
            % 
            % centerline_dist_2_start = cumsum(sqrt(sum((pathidxyz(2:end,:) - pathidxyz(1:end-1,:)).^2, 2)));
            % centerline_dist_2_start = [0;centerline_dist_2_start(:)];
            % ratio_neck_head = zeros(length(width),1) + nan;
            % for i = 2:length(ratio_neck_head)-1
            %     ratio_neck_head(i) = width(i+1) - width(i);
            % end
            width_max_point = find(width_eachp == max(width_eachp), 1,'last');

            cylinderMap = zeros(length(idSelected),1) + nan;%indicates whether it is at the head part that 
            for j = 1:length(cylinderMap)
                cur_idxyz = Vnieuwx(idSelected(j),:);
                nearest_centerlinep = find(sum((pathidxyz - cur_idxyz).^2,2) == min(sum((pathidxyz - cur_idxyz).^2,2)));
                if(nearest_centerlinep <= width_max_point)
                    cylinderMap(j) = 0;

                end
            end
            % geo_score = normalize(sphereMap, 'range') + normalize(cylinderMap, 'range');
            % idSelected = idSelected(~isnan(cylinderMap));
            curvaturemap = newCurvature(idSelected);
            curvaturemap(isnan(cylinderMap)) = nan;
            sphereMap(isnan(cylinderMap)) = nan;
    %%=========================================================================================================
            % combine the curvature score and the geometry score
            [~,p] = sort(curvaturemap);
            curvatureOrder = curvaturemap;
            curvatureOrder(p) = 1:length(p);
            [~,p] = sort(sphereMap); %nans are at the last several positions
            geometryOrder = sphereMap;
            geometryOrder(p) = 1:length(p); 
            geometryOrder(isnan(sphereMap)) = nan;
            ccScore = curvatureOrder + geometryOrder + abs(curvatureOrder - geometryOrder);
            % ccScore = curvatureOrder + 2*geometryOrder;
            ccScore(ismember(idSelected, ccHead)) = nan;
            ccScore(ismember(idSelected, ccNeck)) = nan;
            id_point = idSelected(find(ccScore == min(ccScore),1));
    % %%=========================================================================================================
    %         % based on the combined score search for the best cycle on the surface for the sepraration
    %         % find all the cycles, which one has the minimum mean
    %         % score
    %         cutVex3 = comSeg.searchSSpathCutGraph_v2_speedup(Vnieuw3, [1:size(Vnieuw3,1)], mF3,[],[], find(ismember(idSelected, ccHead)), find(ismember(idSelected, ccNeck)));
    
    % %                                         cutVex3 = searchSSpathCutGraph_v2(Vnieuw3, mF3,[],[], find(ismember(idSelected, ccHead)), find(ismember(idSelected, ccNeck)));
    %         mean_scoreAll = cellfun(@(c) mean(ccScore(c)), cutVex3);
    %         id_point = find(mean_scoreAll == min(mean_scoreAll));
    %         coarse_cycle_path = cutVex3{id_point};
    %         % for the same point, the cycle found in the coarse resolution
    %         % graph can be different from the one found in the
    %         % fine graph so a match is needed
    % %%=========================================================================================================
    %         % Run the final segmentation
    %         real_id_point = comSeg.matchSScycle(Vnieuwx,cutVexcellUnique2(idSelected(coarse_cycle_path)), Vnieuw3, coarse_cycle_path);
            [head_volume, cutvex_fin] = comSeg.genCutVolume_speedup(cutVexcellUnique2, id_point,Vnieuw, test_region, allTetrax, allVertex,resx, resy, resz);
            neckx = double(test_region) - head_volume;
            neckx(neckx(:) < 0) = 0;
            neckx_idx = bwconncomp(neckx);
            neckx_idx_cell = neckx_idx.PixelIdxList;
            if(~isempty(neckx_idx_cell))
                size_cc = cellfun(@length, neckx_idx_cell);
                neckID = neckx_idx_cell{size_cc == max(size_cc)};
                neckx = zeros(size(neckx));
                neckx(neckID) = 1;
                if(isempty(intersect(neckID, contactID)))
                    tmp = neckx;
                    neckx = head_volume;
                    head_volume = tmp;
                end
                spineHead=double(head_volume)*3 + double(neckx)*2 + double(shaft_mask);
            else
                spineHead = double(head_volume)*3 + double(shaft_mask);
            end
            tifwrite(uint8(spineHead), fullfile(spine_head_neck_save_folder, num2str(k)));
            % whole_volume = imerode((head_volume*2 + neckx), se);
            % 
            % 
            % [mF2_0, Vnieuw_0, allTetrax_0] = comSeg.genSurfaceTetraGraph_v3(test_region0, xxshift, yyshift, zzshift);
            % %         allVertex = Vnieuw; % contains all the vertex in the volume
            %         %make sure the vertices matrix are all on the surface 
            % 
            % vertex_unique_0 = unique(mF2_0(:));
            % mappingx = zeros(max(vertex_unique_0),1);
            % mappingx(vertex_unique_0(:)) = [1:length(vertex_unique_0(:))];
            % mF2_0 = mappingx(mF2_0);
            % Vnieuw_0 = Vnieuw_0(vertex_unique_0(:),:);
            % % Not for computation, more to ensure the latter comparison to be reliable so that no extra points are included
            % 
            % surface_0_id = sub2ind([lenxtmp,lenytmp, lenztmp], Vnieuw_0(:,1), Vnieuw_0(:,2), Vnieuw_0(:,3));
            % id_head = find(whole_volume(:) == 2);
            % id_neck = find(whole_volume(:) == 1);
            % bd_head_neck = utils.regionGrowxyz3D(id_head, 1, lenxtmp, lenytmp, lenztmp);
            % bd_head_neck = intersect(bd_head_neck, id_neck);
            % bd_head_neck = intersect(bd_head_neck, surface_0_id);
            % [bd_head_neckx, bd_head_necky,bd_head_neckz] = ind2sub([lenxtmp , lenytmp, lenztmp], bd_head_neck);
            % coor_cut_tmp = zeros(length(cutvex_fin),3);
            % for j = 1:length(cutvex_fin)
            %     dist_gd = sqrt((bd_head_neckx - Vnieuw(cutvex_fin(j),1)).^2 + (bd_head_necky - Vnieuw(cutvex_fin(j),2)).^2  + (bd_head_neckz - Vnieuw(cutvex_fin(j),3)).^2 );
            %     id_small_grid_point = find(dist_gd == min(dist_gd),1);
            %     coor_cut_tmp(j,:)  = [2*bd_head_neckx(id_small_grid_point), 2*bd_head_necky(id_small_grid_point), 5*bd_head_neckz(id_small_grid_point)];
            % end
            % [coor_headx,coor_heady,coor_headz]  = ind2sub([lenxtmp , lenytmp, lenztmp], id_head);
            % coor_head = [coor_headx(:)*2, coor_heady(:)*2, coor_headz(:)*5];
            % [coor_neckx, coor_necky, coor_neckz] = ind2sub([lenxtmp , lenytmp, lenztmp], id_neck);   
            % coor_neck = [coor_neckx(:)*2, coor_necky(:)*2, coor_neckz(:)*5];
            % coor_head_cell{k} = coor_head;
            % coor_neck_cell{k} = coor_neck;
            % save(fullfile(outputFolder, [namex,'_head.mat']), 'coor_head');
            % save(fullfile(outputFolder, [namex,'_neck.mat']), 'coor_neck');
        catch ME
            % warning('Problem occured at %s',[num2str(k)]);
            continue;
        end
    end
    % if(~isempty(outputFolder))
    %     for k = 1:length(coor_neck_cell)
    %         if(~isempty(coor_neck_cell{k}))
    %             coor_head = coor_head_cell{k};
    %             coor_neck = coor_neck_cell{k};
    %             file_name = files_x(k).name;
    %             namex = file_name(1:end-4);
    %             save(fullfile(outputFolder, [namex,'_head.mat']), 'coor_head');
    %             save(fullfile(outputFolder, [namex,'_neck.mat']), 'coor_neck');
    %         else
    %             disp(files_x(k).name)
    %         end
    %     end
    % end
    toc;
    end