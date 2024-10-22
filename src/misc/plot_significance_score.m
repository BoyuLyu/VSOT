function [position_cell, pvalueArray] = plot_significance_score(bar_plot_mean, bar_plot_test, barhandle)
    group_num = length(bar_plot_mean);
    
    variable_num = size(bar_plot_mean{1}, 2);
    position_cell = [];
    pvalueArray = [];
    for i = 1:(group_num - 1)
        for j = (i+1):group_num
            for k = 1:variable_num
                tmp_testi = bar_plot_test{i};
                tmp_testj = bar_plot_test{j};
                x1 = barhandle(i).XEndPoints(k); % X position of bar 1 in group 1
                x2 = barhandle(j).XEndPoints(k); % X position of bar 2 in group 1
                position_cell = [position_cell, {[x1, x2]}];
                [~,p1] = ttest2(tmp_testi(:,k), tmp_testj(:,k));
                if(p1 > 0.05)
                    p1 = nan;
                end
                pvalueArray = [pvalueArray, p1];
            end

        end
    end





end