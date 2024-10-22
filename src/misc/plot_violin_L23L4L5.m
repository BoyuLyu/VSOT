function plot_violin_L23L4L5(singleSynHeadVolume2, singleSynHeadVolume4, singleSynHeadVolume5, type)
    if(type == 'head_volume')
        singleSynHeadVolume2(singleSynHeadVolume2 == 0) = [];
        singleSynHeadVolume2(singleSynHeadVolume2>quantile(singleSynHeadVolume2, 0.99)) = [];
        singleSynHeadVolume4(singleSynHeadVolume4 == 0) = [];
        singleSynHeadVolume4(singleSynHeadVolume4>quantile(singleSynHeadVolume4, 0.99)) = [];
        singleSynHeadVolume5(singleSynHeadVolume5 == 0) = [];
        singleSynHeadVolume5(singleSynHeadVolume5>quantile(singleSynHeadVolume5, 0.99)) = [];
        x = [singleSynHeadVolume2(:);singleSynHeadVolume4(:);singleSynHeadVolume5(:)];
        x = x.*16.*16.*40/10^9;
        y = [zeros(length(singleSynHeadVolume2(:)),1);ones(length(singleSynHeadVolume4(:)),1);ones(length(singleSynHeadVolume5(:)),1).*2];
        
        figure; violinplot(x,y,'ShowMean',true,'ShowMedian',true); title(['volume of dendrite spine head'])
        set(gca,'xtick',[1,2,3],'xticklabel',{'L2', 'L4', 'L5'})
        ylabel('volume(\mu m^3)')
        [h,p1,ci,stats] = ttest2(singleSynHeadVolume2(:),singleSynHeadVolume4(:));
        [h,p2,ci,stats] = ttest2(singleSynHeadVolume4(:),singleSynHeadVolume5(:));
        [h,p3,ci,stats] = ttest2(singleSynHeadVolume2(:),singleSynHeadVolume5(:));
        sigstar({[1,2],[2,3], [1,3]},[p1,p2,p3])
    end


end