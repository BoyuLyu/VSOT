function centerPointsID = find_valid_region(pathIDxyz, idSelected, cutCyclecell, surfaceVertex,resx, resy, resz)
    
    centerPointsID = nan(length(idSelected), 1);
    % check the center of all the polygons that go through the selected points
    % the goal is to ensure the center point of the polygon is also within
    % the limit
    for i = 1:length(idSelected)
        id = idSelected(i);
        if(~isempty(cutCyclecell{id}))
            centerxyz = mean(surfaceVertex(cutCyclecell{id},:),1);
            % check the center of the polygon
            distAll = sqrt((pathIDxyz(:,1) - centerxyz(1)).^2 + (pathIDxyz(:,2) - centerxyz(2)).^2 + (pathIDxyz(:,3) - centerxyz(3)).^2);
            id2 = find(distAll == min(distAll),1);
            centerPointsID(i) = id2;
        end

    end







end