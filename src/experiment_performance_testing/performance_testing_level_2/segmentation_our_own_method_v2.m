clear


addpath('../resources/curvatures/')
addpath('../resources/ImprovedSurfaceSmooth/')
addpath('../resources/CurvatureEstimation')
addpath('/home/boyu/Documents/edt_mex/edt_mex/edt_mex')
addpath('../resources/TAUBIN/TAUBIN/')
addpath('/home/boyu/Documents/iso2mesh/')
addpath('../resources/data_to_off/')
addpath('/home/boyu/Documents/src_mex/mex_EM_analysis/mex_EM_analysis')

offFolder = '/work/boyu/EM_astrocyte/test_segmentation_samples/gt_300/surface_off_300';
annotationFolder = '/work/boyu/EM_astrocyte/test_segmentation_samples/gt_300/annotation_json_300';
annotationFolder_gt_curves = '/work/boyu/EM_astrocyte/test_segmentation_samples/gt_300/cut_cycle_ID_300';
listx = dir([offFolder, '/*.off']);
our_method_coordinate_folder = '/work/boyu/EM_astrocyte/test_segmentation_samples/our_segmentation_result_cut_coordinates';
our_method_cut_result_folder = '/work/boyu/EM_astrocyte/dendrite_segmentation_peer_methods/our_method/result';
parfor j = 1:length(listx)
    
    namex = listx(j).name;
    namex = namex(1:end-4);
    disp(namex)


    [Pts,Tri] = read_off(fullfile(offFolder, [namex , '.off']));
    Tri = Tri';
    Pts = Pts';

% 
% figure;trisurf(Tri,Pts(:,1),Pts(:,2),Pts(:,3),'Facecolor','red','FaceAlpha',0.1);
% 
% figure;trisurf(mF2,Vnieuwx(:,1),Vnieuwx(:,2),Vnieuwx(:,3),'FaceVertexCData',ccScore_fin,'Facecolor','interp');
% hold on;trisurf(mF2,Vnieuwx(:,1),Vnieuwx(:,2),Vnieuwx(:,3), 'Facecolor','blue','FaceAlpha',0.1)

%% find a cycle in the mesh for the annotation that is close to the cycle in the tetrahedra mesh
% for each vertex in the tetra-mesh ,find its closest vertex on the smooth
% mesh
    if(exist(fullfile(our_method_coordinate_folder,[namex,'.off_head.mat']), 'file'))
        aa = load(fullfile(our_method_coordinate_folder,[namex,'.off_head.mat']));
        coor_head_tmp = aa.coor_head;
        bb = load(fullfile(our_method_coordinate_folder,[namex,'.off_neck.mat']));
        coor_neck_tmp = bb.coor_neck;

        face_center_idxyz = (Pts(Tri(:,1),:)+ Pts(Tri(:,2),:)+ Pts(Tri(:,3),:))/3;
        % check which face center should belong to the head part which should
        % belong to the neck part
        %surface coordinates are multiplied by [2,2,5]
        face_label = zeros(size(face_center_idxyz,1),1);
%         face_center_idxyz = face_center_idxyz./([2,2,5]);
        for m = 1:size(face_center_idxyz,1)
            dist2head = min(vecnorm(face_center_idxyz(m,:) - coor_head_tmp, 2,2));
            dist2neck = min(vecnorm(face_center_idxyz(m,:) - coor_neck_tmp, 2,2));
            if(dist2head < dist2neck)
                face_label(m) = 1;
            else
                face_label(m) = 2;
            end
        end
        node_colorMap = zeros(length(Pts),4);
        node_colorMap2 = zeros(length(Pts),1);
        colormapx = [255, 0, 0, 255;0, 255, 0, 255;0, 0, 255, 255];
        edge_all = [[Tri(:,1), Tri(:,2), [1:size(Tri,1)]']; [Tri(:,1), Tri(:,3), [1:size(Tri,1)]'];[Tri(:,3), Tri(:,2), [1:size(Tri,1)]']];
        edge_all(:,1:2) = sort(edge_all(:,1:2), 2);
        edge_all = sortrows(edge_all,[1,2]);
        edge_all(:,3) = face_label(edge_all(:,3));
        a1 = [1:2:size(edge_all,1)];
        a2 = [2:2:size(edge_all,1)];
        id0 = find(edge_all(a1,3) ~= edge_all(a2,3));
        node_list = edge_all(a1(id0),1:2);
        node_list = unique(node_list(:));
        writematrix(node_list(:), fullfile(our_method_cut_result_folder,[namex,'_cut_v2.txt']));
    end
    
    % 
%     node_colorMap2 = zeros(size(Pts,1),1);
%     node_colorMap2(node_list(:)) = 1;
% figure; trisurf(Tri, Pts(:,1), Pts(:,2), Pts(:,3), 'FaceVertexCData',node_colorMap2,'Facecolor','interp') 

end
