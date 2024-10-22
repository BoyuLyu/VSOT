function plot_kde_contour_3(x1, y1, x2, y2, x3, y3)
    x_all = [x1(:); x2(:); x3(:)];
    y_all = [y1(:); y2(:); y3(:)];
    x1_range = linspace(min(x_all), max(x_all), 100);
    y1_range = linspace(min(y_all), max(y_all), 100);
    [X1, Y1] = meshgrid(x1_range, y1_range);    
    density1 = ksdensity([x1(:), y1(:)], [X1(:), Y1(:)]);
    density1 = reshape(density1, size(X1)); 

    x2_range = linspace(min(x_all), max(x_all), 100);
    y2_range = linspace(min(y_all), max(y_all), 100);
    [X2, Y2] = meshgrid(x2_range, y2_range);
    density2 = ksdensity([x2(:), y2(:)], [X2(:), Y2(:)]);
    density2 = reshape(density2, size(X2));

    x3_range = linspace(min(x_all), max(x_all), 100);
    y3_range = linspace(min(y_all), max(y_all), 100);
    [X3, Y3] = meshgrid(x3_range, y3_range);
    density3 = ksdensity([x3(:), y3(:)], [X3(:), Y3(:)]);
    density3 = reshape(density3, size(X3));


    figure;
    contour(X1, Y1, density1, 'LineWidth', 1, 'LineColor', 'r');
    hold on;
    contour(X2, Y2, density2, 'LineWidth', 1,'LineColor', 'g');
    hold on;
    contour(X3, Y3, density3, 'LineWidth', 1, 'LineColor', 'b');
    legend('L2/3', 'L4', 'L5');
    hold off;


end
