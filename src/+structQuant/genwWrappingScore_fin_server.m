function [perimeterRatio,perimeterWeightedWrappingArea,postSynapseTouchingArea,postSynapseTouchingRatio,preSynapseTouchingArea,preSynapseTouchingRatio,...
    headNeckTouchingArea, headNeckTouchingRatio] = genwWrappingScore_fin_server(cleftID, lenx, leny, lenz, maskAstro, segMaskFull,postSynNewID, pre_label, spineID, spineHNROI, resx, resy, resz)
[cleftIDx,cleftIDy, cleftIDz] = ind2sub([lenx, leny, lenz], cleftID);
cleftROI = zeros(lenx, leny, lenz);
cleftROI(cleftID) = 1;
centerx = round(mean(cleftIDx));
centery = round(mean(cleftIDy));
centerz = round(mean(cleftIDz));
perimeterRatio = 0;
perimeterWeightedWrappingArea = 0;
postSynapseTouchingArea = 0;
postSynapseTouchingRatio = 0;
preSynapseTouchingArea = 0;
preSynapseTouchingRatio = 0;
headNeckTouchingArea = [0, 0];
headNeckTouchingRatio = [0, 0];

if(sum(maskAstro(:)) ==0)
    perimeterRatio = 0;
    perimeterWeightedWrappingArea = 0;
    postSynapseTouchingArea = 0;
    postSynapseTouchingRatio = 0;
    preSynapseTouchingArea = 0;
    preSynapseTouchingRatio = 0;
    headNeckTouchingArea = [0, 0];
    headNeckTouchingRatio = [0, 0];
else
%     segMasktmp = segMaskFull(max(centerx - 50,1):min(centerx+50, lenx), max(centery-50,1):min(centery+50, leny), max(centerz-20,1): min(centerz+20, lenz));
    se2 = strel('sphere', 1);
    maskPre = double(segMaskFull == pre_label);
    maskPost = zeros(lenx, leny, lenz);
    maskPost(spineID) = 1;
    maskPre_big = imdilate(maskPre, se2);
    maskPost_big = imdilate(maskPost, se2);
    maskCenter = maskPre_big.*maskPost_big;

    mask_spine_big = maskPost_big.*spineHNROI;
    mask_spine_big = imdilate(mask_spine_big, se2).*maskPost_big - spineHNROI;
    mask_spine_big = mask_spine_big.*(1 - maskCenter); % the boundary of the 
    mask_spine_big(mask_spine_big<0) = 0;
    if(sum(mask_spine_big(:) == 3) > 0 && sum(mask_spine_big(:) == 2) > 0)
        mask_spine_head_id = find(mask_spine_big(:) == 3);
        mask_spine_neck_id = find(mask_spine_big(:) == 2);
        headNeckTouchingArea = [sum(maskAstro(mask_spine_head_id)), sum(maskAstro(mask_spine_neck_id))];
        headNeckTouchingRatio = [sum(maskAstro(mask_spine_head_id))/ length(mask_spine_head_id), sum(maskAstro(mask_spine_neck_id))/ length(mask_spine_neck_id)];
    end
    maskPre = maskPre(max(centerx - 50,1):min(centerx+50, lenx), max(centery-50,1):min(centery+50, leny), max(centerz-20,1): min(centerz+20, lenz));
    maskAstro = maskAstro(max(centerx - 50,1):min(centerx+50, lenx), max(centery-50,1):min(centery+50, leny), max(centerz-20,1): min(centerz+20, lenz));
    maskPre_big = maskPre_big(max(centerx - 50,1):min(centerx+50, lenx), max(centery-50,1):min(centery+50, leny), max(centerz-20,1): min(centerz+20, lenz));
    maskPost_big = maskPost_big(max(centerx - 50,1):min(centerx+50, lenx), max(centery-50,1):min(centery+50, leny), max(centerz-20,1): min(centerz+20, lenz));
    maskPost = maskPost(max(centerx - 50,1):min(centerx+50, lenx), max(centery-50,1):min(centery+50, leny), max(centerz-20,1): min(centerz+20, lenz));
    maskCenter = maskPre_big.*maskPost_big;
    maskCenter_id = find(maskCenter(:) == 1);
    maskPre_dist = structQuant.getGeoDist_server(maskPre_big - maskPre,maskCenter_id, resx, resy, resz);
    pre_surface = double(maskPre_dist > 0 & maskPre_dist < 500);
    pre_surface_id = find(pre_surface(:) == 1);
    maskPost_dist = structQuant.getGeoDist_server(maskPost_big - maskPost,maskCenter_id, resx, resy, resz);
    post_surface = double(maskPost_dist > 0 & maskPost_dist < 500);
    post_surface_id = find(post_surface(:) == 1);
    

    
    maskCenter = maskCenter.*(1 - maskPre);
    maskCenter = maskCenter.*(1 - maskPost); % the perimeter of the cleft
    tmpid = find(maskCenter(:) == 1);


    if(isempty(tmpid))
        perimeterRatio = nan;
        perimeterWeightedWrappingArea = nan;
        postSynapseTouchingArea = nan;
        postSynapseTouchingRatio = nan;
        preSynapseTouchingArea = nan;
        preSynapseTouchingRatio = nan;
        headNeckTouchingArea = nan;
        headNeckTouchingRatio = nan;
    else
        perimeterRatio = sum(maskAstro(tmpid))/length(tmpid);
        postSynapseTouchingArea = sum(maskAstro(pre_surface_id));
        postSynapseTouchingRatio = postSynapseTouchingArea/length(pre_surface_id);
        preSynapseTouchingArea = sum(maskAstro(post_surface_id));
        preSynapseTouchingRatio = preSynapseTouchingArea/ length(post_surface_id);
        pre_surface_idx = pre_surface_id(maskAstro(pre_surface_id) == 1);
        post_surface_idx = post_surface_id(maskAstro(post_surface_id) == 1);
        allDistAstro = [maskPre_dist(pre_surface_idx); maskPost_dist(post_surface_idx)];
        perimeterWeightedWrappingArea = sum(exp(-18.*(allDistAstro/1000).^2));
    end
end
    










end