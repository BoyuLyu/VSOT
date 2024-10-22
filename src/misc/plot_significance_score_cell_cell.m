function [position_cell, pvalueArray] = plot_significance_score_cell_cell(bar_plot_mean, bar_plot_test, barhandle)
    group_num = length(bar_plot_mean);
    
    variable_num = length(bar_plot_mean{1});
    position_cell = [];
    pvalueArray = [];
    for i = 1:(group_num - 1)
        for j = (i+1):group_num
            for k = 1:variable_num
                tmp_testi = bar_plot_test{i, k};
                tmp_testj = bar_plot_test{j, k};
                if(variable_num == 1)
                    x1 = barhandle.XEndPoints(i); % X position of bar 1 in group 1
                    x2 = barhandle.XEndPoints(j); % X position of bar 2 in group 1
                    position_cell = [position_cell, {[x1, x2]}];
                    [~,p1] = ttest2(tmp_testi, tmp_testj);
                    if(p1>=0.05)
                        p1 = nan;
                    end
                    pvalueArray = [pvalueArray, p1];


                else

                    x1 = barhandle(i).XEndPoints(k); % X position of bar 1 in group 1
                    x2 = barhandle(j).XEndPoints(k); % X position of bar 2 in group 1
                    position_cell = [position_cell, {[x1, x2]}];
                    [~,p1] = ttest2(tmp_testi, tmp_testj);
                    if(p1>=0.05)
                        p1 = nan;
                    end
                    pvalueArray = [pvalueArray, p1];
                end
            end

        end
    end





end