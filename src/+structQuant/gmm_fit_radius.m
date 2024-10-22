function [small_group_id, large_group_id] = gmm_fit_radius(dendrite_radius2)
    id2_non_zero_radius = find(dendrite_radius2 ~= 0);
    dendrite_radius2 = dendrite_radius2(id2_non_zero_radius);
    GMModel = fitgmdist(dendrite_radius2,2);
    ratios = GMModel.ComponentProportion;
    mus = GMModel.mu;
    svars = GMModel.Sigma;
    figure; histogram(dendrite_radius2, 50, 'Normalization','pdf');
    distribution1 = ratios(1)*normpdf(sort(dendrite_radius2(:)),mus(1), sqrt(svars(1)));
    distribution2 = ratios(2)*normpdf(sort(dendrite_radius2(:)),mus(2), sqrt(svars(2)));
    hold on;plot(sort(dendrite_radius2(:)), distribution1, 'LineWidth',3);
    hold on; plot(sort(dendrite_radius2(:)), distribution2,'LineWidth',3);
    title('GMM fitting of the dendrite radius')
    if(mus(1) < mus(2))
        distribution1 = ratios(1)*normpdf(sort(dendrite_radius2(:)),mus(1), sqrt(svars(1)));
        distribution2 = ratios(2)*normpdf(sort(dendrite_radius2(:)),mus(2), sqrt(svars(2)));
    else
        distribution2 = ratios(1)*normpdf(sort(dendrite_radius2(:)),mus(1), sqrt(svars(1)));
        distribution1 = ratios(2)*normpdf(sort(dendrite_radius2(:)),mus(2), sqrt(svars(2)));  
    end
    
    [sorted_dendrite_radius2, id_sorted_2] = sort(dendrite_radius2(:));
    split_point = find(distribution1 <= distribution2,1);

    large_group_id = id_sorted_2(split_point:end);
    small_group_id = setdiff(id_sorted_2, large_group_id);
    large_group_id = id2_non_zero_radius(large_group_id);
    small_group_id = id2_non_zero_radius(small_group_id);
end