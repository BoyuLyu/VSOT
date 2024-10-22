function [geoDist2Starting] = getGeoDist_server(test_region2,startingid,dx, dy, dz)
% addpath('/home/boyu/Documents/src_mex/mex_EM_analysis/mex_EM_analysis')
[lenx, leny, lenz] = size(test_region2);
test_region1D = logical(test_region2(:));
startx = zeros(lenx, leny, lenz,'logical');
startx(startingid) = 1;
startx1D = startx(:);
dist1 = imchaferDist3D(test_region1D,startx1D,lenx, leny, lenz, dx, dy, dz);
geoDist2Starting = reshape(dist1, lenx, leny, lenz);
end