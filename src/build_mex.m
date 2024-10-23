function build_mex(resource_folder)
% read in one example file to ensure the mex functions are compiled correctly

%%=========================================================================
spine_path = fullfile(resource_folder, 'mex_build_test_file/test.tif');
spine = tiffreadVolume(spine_path);

[lenx, leny, lenz] = size(spine);

cd(fullfile(resource_folder,'/edt_mex/edt_mex/edt_mex'))
mex edt_mex.cpp
resx = 2;
resy = 2;
resz = 5;
tmp_den_1d = logical(spine(:) >0);

try
    tmp_dist_1d = edt_mex(tmp_den_1d, lenx, leny, lenz, resx,resy,resz);
    tmp_dist = reshape(tmp_dist_1d, size(spine));
    disp('success')
catch ME
    warning('edt mex function failed to run');
end
%%=========================================================================

cd(fullfile(resource_folder, '/graph_related/graph_mex'))

mex solve_min_cut.cpp

terminalWeights=[
16,0;
13,0;
0,20;
0,4;
];
pInPtr=[
1,2,10.9,4.5;
1,3,12.1,-1.5;
2,3,-1.5,9.3;
2,3,14.2,0;
1,4,0,7.5];
try
    outx = solve_min_cut(terminalWeights, pInPtr);
    disp('success')
catch ME
    warning('min-cut mex function failed to run');
end

%%=========================================================================
cd(fullfile(resource_folder, '/src_mex/mex_EM_analysis/mex_EM_analysis'))

mex imchaferDist3D.cpp

maskstart = zeros(lenx, leny, lenz);
spine_id = find(spine(:) > 0);
maskstart(spine_id(1)) = 1;
maskstart_1d = logical(maskstart(:));
try
    scoreIJ = imchaferDist3D(tmp_den_1d,maskstart_1d, lenx, leny, lenz,resx,resy,resz);
    disp('success')
catch ME
    warning('imchaferDist3D mex function failed to run');
end

%%=========================================================================
cd(fullfile(resource_folder,'/tinymesh_mex/tiny_mesh_mex'))
path1 = fullfile(resource_folder,'/tinymesh/src/tinymesh');
path2 = fullfile(resource_folder,'/tinymesh/src/tinymesh/ext/eigen');
path3 = fullfile(resource_folder,'/tinymesh/build/lib');
eval(['mex CXXFLAGS=''$CXXFLAGS -std=c++14'' taubin_smooth_tiny_mesh_mex.cpp -I',path1,' -I',path2,' -L',path3,' -ltinymesh']);

[spine_idx, spine_idy,spine_idz] = ind2sub(size(spine), find(spine(:) > 0));
spine_idx = 2.*spine_idx(:);
spine_idy = 2.*spine_idy(:);
spine_idz = 5.*spine_idz(:);
newidTMP = [spine_idx(:), spine_idy(:), spine_idz(:)];
bbk = boundary(spine_idx,spine_idy,spine_idz,1);
try
    [Vnieuw2_out,mF2_out] = taubin_smooth_tiny_mesh_mex(newidTMP, bbk, 0.8, 0.53, 40); % mex file compiled from the c++ code for speedup
    disp('success')
catch ME
    warning('taubin smooth mex function failed to run');
end



% mex CXXFLAGS='$CXXFLAGS -std=c++14' taubin_smooth_tiny_mesh_mex.cpp -I/work/boyu/EM_astrocyte/VSOT_code/VSOT/resources/mex_build_test/tinymesh/src/tinymesh ...
%     -I/work/boyu/EM_astrocyte/VSOT_code/VSOT/resources/mex_build_test/tinymesh/src/tinymesh/ext/eigen ...
%     -L/work/boyu/EM_astrocyte/VSOT_code/VSOT/resources/mex_build_test/tinymesh/build/lib ...
%     -ltinymesh

end