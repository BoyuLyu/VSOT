
tifResultPath = 'D:\Guoqiang_lab\VSOT_samples\VSOT_data_paired_with_codes\data\performance_testing_level_2_segmentation\large_dataset_3_annotator_L2345_w_stubby\segmentation_results\vsot_result_tif_volume';


files = dir(fullfile(tifResultPath, "*.tif"));
headVolumex_all = zeros(length(files), 1);
headMeanRadiusx_all = zeros(length(files), 1);
neckLengthx_all = zeros(length(files), 1);
neckSectionx_all = zeros(length(files), 1);
neckMeanRadiusx_all = zeros(length(files), 1);
neckRadiusSTD_all = zeros(length(files), 1);
for i = 1:length(files)
    try
        filename = files(i).name;
        filename = strsplit(filename, '.');
        filename = filename{1};
        disp(filename)
        spineHNROI = tiffreadVolume(fullfile(tifResultPath, files(i).name));
        [lenx, leny, lenz] = size(spineHNROI);
        resx = 16;
        resy = 16;
        resz = 40;
        [headVolumex,headMeanRadiusx,neckLengthx, neckSectionx, neckMeanRadiusx, neckRadiusSTD] = structQuant.genDendriteSpineScore_nocleft([], lenx, leny, lenz,spineHNROI, resx,resy, resz);
        headVolumex_all(i) = headVolumex;
        headMeanRadiusx_all(i) = headMeanRadiusx;
        neckLengthx_all(i) = neckLengthx;
        neckSectionx_all(i) = neckSectionx;
        neckMeanRadiusx_all(i) = neckMeanRadiusx;
        neckRadiusSTD_all(i) = neckRadiusSTD;
    catch ME
        continue;
    end
end


