function [single_synapse_struct, double_synapse_struct] = ...
    gen_quantification_spine_sub_function(spine_mask,cleft_mask_small,astro_mask_small,seg_mask_small,postSynNewID, resx,resy, resz)
    spine_mask_id = find(spine_mask(:) > 1);
    spineHNROI = spine_mask;
    cleftROI = bwlabeln(cleft_mask_small);
    cleftROIidx = label2idx(cleftROI);
    labelx = cleftROI(spine_mask_id);
    labelx(labelx==0)=[];
    cleftSize = cellfun(@length, cleftROIidx);
    cleftSize = cleftSize(:);
    [C,ia,ic] = unique(labelx);
    a_counts = accumarray(ic,1);
    value_counts = [C, a_counts];
    value_counts((a_counts./cleftSize(C) < 0.2) | cleftSize(C) < 10,:) = []; % avoid including the clefts that are only partially touching the spine
    single_synapse_struct = struct();
    double_synapse_struct = struct();
    [lenx, leny, lenz] = size(spine_mask);
    % certain spine mask contains an extra part of the dendrite, remove it
    all_roi = bwlabeln(spineHNROI > 0);
    all_roi_idx = label2idx(all_roi);
    all_roi_idx = all_roi_idx(:);
    target_indicators = cellfun(@(x) sum(spineHNROI(x) == 3), all_roi_idx);
    target_mask = zeros(size(spineHNROI), 'uint8');
    target_mask(spine_mask_id(target_indicators > 0)) = 1;
    spineHNROI = spineHNROI.*target_mask;

    if(size(value_counts,1) == 0||size(value_counts,1) >=2)
        % for certain spines without any cleft touching save only the
        % quantification related to the dendrite spine
       synpatic_double_single_indicator = size(value_counts,1);
       cleftID = [];
       single_synapse_struct.synpatic_double_single_indicator = synpatic_double_single_indicator;
        [headVolumex,headMeanRadiusx,neckLengthx, neckSectionx, neckMeanRadiusx] = structQuant.genDendriteSpineScore_nocleft(spine_mask_id, lenx, leny, lenz,spineHNROI, resx,resy, resz);
        single_synapse_struct.headVolumex = headVolumex;
        single_synapse_struct.headMeanRadiusx = headMeanRadiusx;
        single_synapse_struct.neckLengthx = neckLengthx;
        single_synapse_struct.neckSectionx = neckSectionx;
        single_synapse_struct.neckMeanRadiusx = neckMeanRadiusx;
        spineID = spine_mask_id;
        [spineIDcoorx, spineIDcoory,spineIDcoorz] = ind2sub([lenx, leny, lenz], spineID(1));
        single_synapse_struct.spineCoordinate = [spineIDcoorx,spineIDcoory,spineIDcoorz];
        single_synapse_struct.singleSynapticCleftSize = length(cleftID);
        single_synapse_struct.sinsperimeterRatio = nan;
        single_synapse_struct.sinsperimeterWeightedWrappingArea = nan;
        single_synapse_struct.sinspostSynapseTouchingArea = nan;
        single_synapse_struct.sinspostSynapseTouchingRatio = nan;
        single_synapse_struct.sinspreSynapseTouchingArea = nan;
        single_synapse_struct.sinspreSynapseTouchingRatio = nan;
        single_synapse_struct.singleSynHeadNeckTouchingArea = nan;
        single_synapse_struct.singleSynHeadNeckTouchingRatio = nan;    
    
    elseif(size(value_counts,1) == 1)
    %                         disp(i)
        cleftID = cleftROIidx{value_counts(1,1)};
        %                         outputMask(cleftID) = 7 + value_counts(1,1); % 7 is to compensate for the difference between the file and the neuroglancer
        pre_post_labelAll = seg_mask_small(cleftID);
        pre_post_labelAll(pre_post_labelAll == postSynNewID) = [];
        pre_post_labelAll(pre_post_labelAll == 0) = [];
        pre_post_label = mode(pre_post_labelAll);
        pre_post_label(pre_post_label == 0) = [];
        if(length(pre_post_label) == 1)
            pre_label = pre_post_label;
        %                             outputMask(seg_mask_smallidx{pre_label}) = pre_label;
            synpatic_double_single_indicator = 1;
            single_synapse_struct.synpatic_double_single_indicator = synpatic_double_single_indicator;
            [headVolumex,headMeanRadiusx,neckLengthx, neckSectionx, neckMeanRadiusx] = structQuant.genDendriteSpineScore_server(cleftID, spine_mask_id, lenx, leny, lenz,spineHNROI, resx, resy, resz);
            single_synapse_struct.headVolumex = headVolumex;
            single_synapse_struct.headMeanRadiusx = headMeanRadiusx;
            single_synapse_struct.neckLengthx = neckLengthx;
            single_synapse_struct.neckSectionx = neckSectionx;
            single_synapse_struct.neckMeanRadiusx = neckMeanRadiusx;
            [perimeterRatio,perimeterWeightedWrappingArea,postSynapseTouchingArea,postSynapseTouchingRatio,preSynapseTouchingArea,preSynapseTouchingRatio,...
                headNeckTouchingArea, headNeckTouchingRatio] = structQuant.genwWrappingScore_fin_server(cleftID, lenx, leny, lenz, astro_mask_small, seg_mask_small,postSynNewID, pre_label, spine_mask_id, spineHNROI, resx, resy, resz);
            single_synapse_struct.singleSynapticCleftSize = length(cleftID);
            single_synapse_struct.sinsperimeterRatio = perimeterRatio;
            single_synapse_struct.sinsperimeterWeightedWrappingArea = perimeterWeightedWrappingArea;
            single_synapse_struct.sinspostSynapseTouchingArea = postSynapseTouchingArea;
            single_synapse_struct.sinspostSynapseTouchingRatio = postSynapseTouchingRatio;
            single_synapse_struct.sinspreSynapseTouchingArea = preSynapseTouchingArea;
            single_synapse_struct.sinspreSynapseTouchingRatio = preSynapseTouchingRatio;
            single_synapse_struct.singleSynHeadNeckTouchingArea = headNeckTouchingArea;
            single_synapse_struct.singleSynHeadNeckTouchingRatio = headNeckTouchingRatio;
        %                             singleSynpostSynapseTouchingArea = [singleSynpostSynapseTouchingArea;postContactSize];
        %                             singleSynpreSynapseTouchingArea = [singleSynpreSynapseTouchingArea;preContactSize];
        %                             singleSynID = [singleSynID;7 + value_counts(1,1)];
        %                             singleSynHeadVolume = [singleSynHeadVolume; sum(spineHNROI(spineROIidx{i}) == 3)];
            spineID = spine_mask_id;
            [spineIDcoorx, spineIDcoory,spineIDcoorz] = ind2sub([lenx, leny, lenz], spineID(1));
            single_synapse_struct.spineCoordinate = [spineIDcoorx,spineIDcoory,spineIDcoorz];
        else
            warning('single synpase pre root ID not right')
        end
    elseif(size(value_counts,1) == 2)
        %                         disp('found one!')
        cleftID1 = cleftROIidx{value_counts(1,1)};
        %                         outputMask(cleftID1) = 7 + value_counts(1,1);
        pre_post_label1All = seg_mask_small(cleftID1);
        pre_post_label1All(pre_post_label1All == postSynNewID) = [];
        pre_post_label1All(pre_post_label1All == 0) = [];
        pre_post_label1 = mode(pre_post_label1All);
        pre_post_label1(pre_post_label1 == 0) = [];
        cleftID2 = cleftROIidx{value_counts(2,1)};
        %                         outputMask(cleftID2) = 7 + value_counts(2,1);
        pre_post_label2All = seg_mask_small(cleftID2);
        pre_post_label2All(pre_post_label2All == postSynNewID) = [];
        pre_post_label2All(pre_post_label2All == 0) = [];
        pre_post_label2 = mode(pre_post_label2All);
        pre_post_label2(pre_post_label2 == 0) = [];
            if(length(pre_post_label1) == 1 && length(pre_post_label2) == 1)
                synpatic_double_single_indicator = 2;
                double_synapse_struct.synpatic_double_single_indicator = synpatic_double_single_indicator;
                pre_label1 = pre_post_label1;
                pre_label2 = pre_post_label2;
                if(length(cleftID1) > length(cleftID2))
                    [headVolumex,headMeanRadiusx,neckLengthx, neckSectionx, neckMeanRadiusx] = structQuant.genDendriteSpineScore_server(cleftID1, spine_mask_id, lenx, leny, lenz,spineHNROI, resx, resy, resz);
                else
                    [headVolumex,headMeanRadiusx,neckLengthx, neckSectionx, neckMeanRadiusx] = structQuant.genDendriteSpineScore_server(cleftID2, spine_mask_id, lenx, leny, lenz,spineHNROI, resx, resy, resz);
                end

                double_synapse_struct.headVolumex = headVolumex;
                double_synapse_struct.headMeanRadiusx = headMeanRadiusx;
                double_synapse_struct.neckLengthx = neckLengthx;
                double_synapse_struct.neckSectionx = neckSectionx;
                double_synapse_struct.neckMeanRadiusx = neckMeanRadiusx;
                [perimeterRatio1,perimeterWeightedWrappingArea1,postSynapseTouchingArea1,postSynapseTouchingRatio1,preSynapseTouchingArea1,preSynapseTouchingRatio1,...
                    headNeckTouchingArea, headNeckTouchingRatio] = structQuant.genwWrappingScore_v3_server(cleftID1, lenx, leny, lenz, astro_mask_small, seg_mask_small,postSynNewID, pre_label1,spine_mask_id, spineHNROI, resx, resy, resz);
                [perimeterRatio2,perimeterWeightedWrappingArea2,postSynapseTouchingArea2,postSynapseTouchingRatio2,preSynapseTouchingArea2,preSynapseTouchingRatio2,...
                    ~, ~] = structQuant.genwWrappingScore_v3_server(cleftID2, lenx, leny, lenz, astro_mask_small, seg_mask_small,postSynNewID, pre_label2, spine_mask_id, spineHNROI, resx, resy, resz);
                double_synapse_struct.doubleSynapticCleftSize = [length(cleftID1), length(cleftID2)];
                double_synapse_struct.dousperimeterRatio = [perimeterRatio1, perimeterRatio2];
                double_synapse_struct.dousperimeterWeightedWrappingArea = [perimeterWeightedWrappingArea1, perimeterWeightedWrappingArea2];
                double_synapse_struct.douspostSynapseTouchingArea = [postSynapseTouchingArea1, postSynapseTouchingArea2];
                double_synapse_struct.douspostSynapseTouchingRatio = [postSynapseTouchingRatio1, postSynapseTouchingRatio2];
                double_synapse_struct.douspreSynapseTouchingArea = [preSynapseTouchingArea1, preSynapseTouchingArea2];
                double_synapse_struct.douspreSynapseTouchingRatio = [preSynapseTouchingRatio1, preSynapseTouchingRatio2];
                double_synapse_struct.doubleSynHeadNeckTouchingArea = headNeckTouchingArea;
                double_synapse_struct.doubleSynHeadNeckTouchingRatio = headNeckTouchingRatio;
                tmp_position = zeros(1,2);
                if(mode(spineHNROI(cleftID1)) == 2)
                    tmp_position(1,1) = tmp_position(1,1) + 1;
                elseif(mode(spineHNROI(cleftID1)) == 3)
                    tmp_position(1,2) = tmp_position(1,2) + 1;
                end
                if(mode(spineHNROI(cleftID2)) == 2)
                    tmp_position(1,1) = tmp_position(1,1) + 1;
                elseif(mode(spineHNROI(cleftID2)) == 3)
                    tmp_position(1,2) = tmp_position(1,2) + 1;
                end                            
                double_synapse_struct.doubleSynPosition = tmp_position;
                spineID = spine_mask_id;
                [spineIDcoorx, spineIDcoory,spineIDcoorz] = ind2sub([lenx, leny, lenz], spineID(1));
                double_synapse_struct.spineCoordinate = [spineIDcoorx,spineIDcoory,spineIDcoorz];
            else
                warning('double synpase pre root ID not right')                      
            end


    end
end