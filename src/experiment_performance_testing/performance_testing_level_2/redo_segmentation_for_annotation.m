
clear
rootFolder = '/work/boyu/EM_astrocyte/astro_11_28_16_16_40/';
neuronListFolder = '/work/boyu/EM_astrocyte/astro_11_28_64_64_80/';
outputFolder= '/work/boyu/EM_astrocyte/test_segmentation_samples/our_segmentation_result_cut_coordinates';
addpath('../resources/curvatures/')
addpath('../resources/ImprovedSurfaceSmooth/')
addpath('../resources/CurvatureEstimation')
addpath('/home/boyu/Documents/edt_mex/edt_mex/edt_mex')
addpath('../resources/TAUBIN/TAUBIN/')
addpath('/home/boyu/Documents/tinymesh_mex/tiny_mesh_mex')
addpath('/home/boyu/Documents/iso2mesh/')
addpath('../resources/data_to_off/')
addpath('../resources/data_to_off/')
addpath('/home/boyu/Documents/graph_related/graph_mex/')
xxshift = zeros(3,3);
yyshift = zeros(3,3);
for i = -1:1
    for j = -1:1
        xxshift((i+2), (j+2)) = i;
        yyshift((i+2), (j+2)) = j;
    end
end
offFolder = '/work/boyu/EM_astrocyte/test_segmentation_samples/gt_300/surface_off_300';
listx = dir([offFolder, '/*.off']);


id_all = zeros(length(listx), 6); %astrocyte ID, ix, iy, iz, neuronID, dendrite spine ID

for i = 1:length(listx)
    namex = listx(i).name;
    namex_split = strsplit(namex, '.');
    id_str = namex_split{1};
    id_four = strsplit(id_str, '_');
    id_all(i,1) = str2double(id_four{1});
    folder_id = id_four{2};
    id_all(i,2) = str2double(folder_id(1));    
    id_all(i,3) = str2double(folder_id(2));   
    id_all(i,4) = str2double(folder_id(3));
    id_all(i,5) = str2double(id_four{3});
    id_all(i,6) = str2double(id_four{4});
end
coor_cut = cell(length(listx),1);
% for mm = 1:length(listx)
for mm = 1:length(listx)
    namex = listx(mm).name;
    try
        nn = id_all(mm,1);
        ix = id_all(mm,2);
        iy = id_all(mm,3);
        iz = id_all(mm,4);
        m = id_all(mm, 5);
        i = id_all(mm, 6);
    
        disp(nn)
        fullSegfolder_root = [rootFolder,'astro_', num2str(nn), '_minnie65'];
        neuronList = [neuronListFolder, 'astro_', num2str(nn), '_minnie65/', 'top20_neuron_id_no_soma.txt'];
        opts = delimitedTextImportOptions("NumVariables", 1);
    
        % Specify range and delimiter
        opts.DataLines = [1, Inf];
        opts.Delimiter = ",";
        
        % Specify column names and types
        opts.VariableNames = "VarName1";
        opts.VariableTypes = "uint64";
        
        % Specify file level properties
        opts.ExtraColumnsRule = "ignore";
        opts.EmptyLineRule = "read";    
        neuronList_str = table2array(readtable(neuronList, opts));
        [xxshift, yyshift, zzshift] = genShiftMatrixTetra;
        fullSegfolder = fullfile(fullSegfolder_root, num2str(ix), num2str(iy),num2str(iz));
        curpsID = neuronList_str(m);
        disp([ix, iy, iz])
        targetFolder = fullSegfolder;
        spineROI = tiffreadVolume(fullfile(targetFolder,['spine_branched_ROI_',num2str(curpsID),'.tif']));
        spineMask0 = tiffreadVolume(fullfile(targetFolder,['spine_new_no_soma_',num2str(curpsID),'.tif']));
        
        shaft_mask = spineMask0 == 1;
        spineMask0 = [];
        spineROIidx = label2idx(spineROI);
        [lenx, leny, lenz] = size(spineROI);
        spine_summary = table2array(readtable(fullfile(targetFolder,['spine_branched_',num2str(curpsID),'_branch_summary.csv'])));
    
        [lenx, leny, lenz] = size(spineROI);
        test_region = spineROI == i;
    %                     disp(i)
        curID0 = spineROIidx{i};
        [idtmpx, idtmpy, idtmpz] = ind2sub([lenx, leny, lenz], curID0);
        bbx = [max(min(idtmpx) -5, 1), min(max(idtmpx) + 5, lenx);max(min(idtmpy) - 5, 1), min(max(idtmpy)+5, leny);max(min(idtmpz) - 5, 1), min(max(idtmpz)+5, lenz)];
        test_region = test_region(bbx(1,1): bbx(1,2),bbx(2,1):bbx(2,2),bbx(3,1): bbx(3,2));
        head_volume = zeros(size(test_region));
        neckx = zeros(size(test_region));
        [lenxtmp, lenytmp, lenztmp] = size(test_region);
        curID1 = find(test_region(:) == 1);
        [curID1x, curID1y, curID1z] = ind2sub([lenxtmp, lenytmp, lenztmp], curID1);
        test_regionID = find(test_region(:) == 1);


        [mF2_0, Vnieuw_0, allTetrax_0] = genSurfaceTetraGraph_v3(test_region, xxshift, yyshift, zzshift);
%         allVertex = Vnieuw; % contains all the vertex in the volume
        %make sure the vertices matrix are all on the surface 
    
        vertex_unique_0 = unique(mF2_0(:));
        mappingx = zeros(max(vertex_unique_0),1);
        mappingx(vertex_unique_0(:)) = [1:length(vertex_unique_0(:))];
        mF2_0 = mappingx(mF2_0);
        Vnieuw_0 = Vnieuw_0(vertex_unique_0(:),:);
%         Vnieuwx = [2*Vnieuw(:,1),2*Vnieuw(:,2),5*Vnieuw(:,3)];% this matrix only contains the vertices on the surface
        se = strel('sphere', 1);
        test_region = imdilate(test_region, se);
        test_region = imclose(test_region, se);
        [mF2, Vnieuw, allTetrax] = genSurfaceTetraGraph_v3(test_region, xxshift, yyshift, zzshift);
        allVertex = Vnieuw; % contains all the vertex in the volume
        %make sure the vertices matrix are all on the surface 
    
        vertex_unique = unique(mF2(:));
        mappingx = zeros(max(vertex_unique),1);
        mappingx(vertex_unique(:)) = [1:length(vertex_unique(:))];
        mF2 = mappingx(mF2);
        Vnieuw = Vnieuw(vertex_unique(:),:);
        Vnieuwx = [2*Vnieuw(:,1),2*Vnieuw(:,2),5*Vnieuw(:,3)];% this matrix only contains the vertices on the surface
        mF2 = checkWindingOrder(mF2);
        [ Vnieuw2_out,mF2_out] = taubin_smooth_tiny_mesh_mex( Vnieuwx, mF2, 0.8, 0.53, 40); 
        mF2_out2 = reshape(mF2_out, 3, [])' + 1;
        mappingx = zeros(size(Vnieuwx,1),1);
        mF2_t = mF2';
        [uniquex,ia,ic]= unique(mF2_t(:),'stable');
        mappingx(uniquex) = [1:length(mappingx)];
        Vnieuw2 = Vnieuw2_out(mappingx,:);
    
        [mF3,Vnieuw3] = reducepatch(mF2,Vnieuwx,min(size(mF2,1), 8000));%????
        [~, idSelected] = ismember(Vnieuw3, Vnieuwx,'rows'); % start using the smaller set of the vertex to generate the two scores
        [newCurvature] = genCurvatureFace4Scale_speedup(mF2, Vnieuw2, 30, idSelected); % 48 cores 2.16 s
        shaftMaskx = shaft_mask(bbx(1,1): bbx(1,2),bbx(2,1):bbx(2,2),bbx(3,1):bbx(3,2));
        [lenxtmp, lenytmp, lenztmp] = size(test_region);
        se = strel('sphere', 1);
        contactSite = test_region.* imdilate(shaftMaskx, se);
        contactID = find(contactSite(:) == 1);    
        [contactIDx, contactIDy, contactIDz] = ind2sub([lenxtmp, lenytmp, lenztmp], contactID);
        nearest_vertex = round(Vnieuw);
        nearest_vertexID = sub2ind([lenxtmp, lenytmp, lenztmp], nearest_vertex(:,1), nearest_vertex(:,2), nearest_vertex(:,3));
        idrm1 = find(test_region(nearest_vertexID) == 0);
        geoDist2Starting = getGeoDist_server(test_region, contactID,2,2,5);
        distCircleCenter = zeros(size(Vnieuw,1),1);
        distVexStartingPoint = zeros(size(Vnieuw,1),1);
        for j = 1:size(nearest_vertex,1)
            
            distVexStartingPoint(j) = geoDist2Starting(nearest_vertex(j,1), nearest_vertex(j,2), nearest_vertex(j,3)); %geodesic distance to the head or neck
        end
        distVexStartingPoint(distVexStartingPoint == 1.0000e+10)= nan;
        neckpoint = find(distVexStartingPoint == min(distVexStartingPoint),1); % the vertex id that close to the starting point
        headpoint = find(distVexStartingPoint == max(distVexStartingPoint), 1); % the vertex id of the furthest point 
        % assign certain points as invalid
        distVexStartingPoint(distVexStartingPoint == 0) = 1;
        ccHead = find(distVexStartingPoint > quantile(distVexStartingPoint,0.98));
        ccNeck = find(distVexStartingPoint <= quantile(distVexStartingPoint,0.1)); % all the faces related to these two will be categorized as real head and neck               
        %% find the shortest cycle
    
        cutVex2_speedup = searchSSpathCutGraph_v2_speedup(Vnieuwx, idSelected, mF2,headpoint,neckpoint, ccHead, ccNeck);
        cutVexcellUnique2 = cleanReptitve_v2(cutVex2_speedup);
        % for the cut paths that have repetitive points, remove
        % the points between the repetitive ones
        %% obtain the sphericity score 
        [surfaceCell, nodeCell] = genCutSurface_speedup(cutVexcellUnique2(idSelected), mF2, headpoint,neckpoint,Vnieuwx,Vnieuw, test_region, allTetrax, allVertex);
        sphereMap = zeros(length(idSelected),1) + nan;
    
        parfor j = 1:length(sphereMap)
    %                                              disp(j)
            if(~isempty(surfaceCell{j}))
                verticesFin = nodeCell{j};
                combinedSurface = surfaceCell{j};
                vectorMF3_1 = verticesFin(combinedSurface(:,3),:) - verticesFin(combinedSurface(:,1),:);
                vectorMF3_2 = verticesFin(combinedSurface(:,3),:) - verticesFin(combinedSurface(:,2),:);
                ss_tmp = cross(vectorMF3_1, vectorMF3_2,2);
                ss = 1/2*sum(sqrt(ss_tmp(:,1).^2 + ss_tmp(:,2).^2 + ss_tmp(:,3).^2));
                volume = abs(calVolfromMesh(combinedSurface, verticesFin));
                volumeCompare2(j) = volume;                                                
                r1 = (volume/(4*pi/3))^(1/3);
                r2 = sqrt(ss/(4*pi));
                sphereMap(j) = abs(r1 - r2)/r2;
            end
        end
    
        sphereMap(sphereMap == 0) = nan;
        %%
        % sphereScoreMap should be acted as a weight applied on
        % the curvature score
        [~,p] = sort(newCurvature(idSelected));
        curvatureOrder = newCurvature(idSelected);
        curvatureOrder(p) = 1:length(p);
        [~,p] = sort(sphereMap);
        geometryOrder = sphereMap;
        geometryOrder(p) = 1:length(p); 
        ccScore = curvatureOrder + geometryOrder + abs(curvatureOrder - geometryOrder);
        ccScore(ismember(idSelected, ccHead)) = nan;
        ccScore(ismember(idSelected, ccNeck)) = nan;
        % find all the cycles, which one has the minimum mean
        % score
        cutVex3 = searchSSpathCutGraph_v2_speedup(Vnieuw3, [1:size(Vnieuw3,1)], mF3,[],[], find(ismember(idSelected, ccHead)), find(ismember(idSelected, ccNeck)));
        mean_scoreAll = cellfun(@(c) mean(ccScore(c)), cutVex3);
        id_point = find(mean_scoreAll == min(mean_scoreAll));
        coarse_cycle_path = cutVex3{id_point};
        % for the same point, the cycle found in the coarse
        % graph can be very different from the one found in the
        % fine graph so a match is needed
        
        real_id_point = matchSScycle(Vnieuwx,cutVexcellUnique2(idSelected(coarse_cycle_path)), Vnieuw3, coarse_cycle_path);
        [head_volume, cutvex_fin] = genCutVolume_speedup(cutVexcellUnique2, idSelected(real_id_point),Vnieuw, test_region, allTetrax, allVertex);
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
        end
        se = strel('sphere', 1);

        whole_volume = imerode((head_volume*2 + neckx), se);
        surface_0_id = sub2ind([lenxtmp,lenytmp, lenztmp], Vnieuw_0(:,1), Vnieuw_0(:,2), Vnieuw_0(:,3));
        id_head = find(whole_volume(:) == 2);
        id_neck = find(whole_volume(:) == 1);
        bd_head_neck = regionGrowxyz3D(id_head, 1, lenxtmp, lenytmp, lenztmp);
        bd_head_neck = intersect(bd_head_neck, id_neck);
        bd_head_neck = intersect(bd_head_neck, surface_0_id);
        [bd_head_neckx, bd_head_necky,bd_head_neckz] = ind2sub([lenxtmp , lenytmp, lenztmp], bd_head_neck);
        coor_cut_tmp = zeros(length(cutvex_fin),3);
        for j = 1:length(cutvex_fin)
            dist_gd = sqrt((bd_head_neckx - Vnieuw(cutvex_fin(j),1)).^2 + (bd_head_necky - Vnieuw(cutvex_fin(j),2)).^2  + (bd_head_neckz - Vnieuw(cutvex_fin(j),3)).^2 );
            id_small_grid_point = find(dist_gd == min(dist_gd),1);
            coor_cut_tmp(j,:)  = [2*bd_head_neckx(id_small_grid_point), 2*bd_head_necky(id_small_grid_point), 5*bd_head_neckz(id_small_grid_point)];
        end
        [coor_headx,coor_heady,coor_headz]  = ind2sub([lenxtmp , lenytmp, lenztmp], id_head);
        coor_head = [coor_headx(:)*2, coor_heady(:)*2, coor_headz(:)*5];
        [coor_neckx, coor_necky, coor_neckz] = ind2sub([lenxtmp , lenytmp, lenztmp], id_neck);   
        coor_neck = [coor_neckx(:)*2, coor_necky(:)*2, coor_neckz(:)*5];

%         bb_candidates = cutvexCell(id_smaller_surface_bbID);
%         len_bb_candidates = cellfun(@length, bb_candidates);
%         bb_candidates(len_bb_candidates == 0) = [];
%         len_bb_candidates(len_bb_candidates == 0) = [];
%         bd_head_neck_cycle = bb_candidates{find(len_bb_candidates == min(len_bb_candidates),1)};
% %         figure;trisurf(Tri,Pts(:,1),Pts(:,2),Pts(:,3),'Facecolor','red','FaceAlpha',0.1);
% %         hold on;plot3(bd_head_neck_x*2, bd_head_neck_y*2, bd_head_neck_z*5, '.'); 
%         bd_head_neck_coor = [Vnieuw_0(bd_head_neck_cycle,1)*2, Vnieuw_0(bd_head_neck_cycle,2)*2, Vnieuw_0(bd_head_neck_cycle,3)*5];
%         coor_cut_tmp = bd_head_neck_coor;
        save(fullfile(outputFolder, [namex,'_head.mat']), 'coor_head');
        save(fullfile(outputFolder, [namex,'_neck.mat']), 'coor_neck');
        
    catch ME
            warning('Problem occured at %s',[num2str(mm)]);
            continue;
    end
end

