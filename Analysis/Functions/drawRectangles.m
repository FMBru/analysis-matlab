function [out] = drawRectangles(centers,size, angle)
%DRAWRECTANGLES 

    if nargin < 3
        angle = 0; % default
    end
    
    for j=1:length(centers(:,1))
        xc = centers(j,1);
        yc = centers(j,2);
        h = size/2; 
    
        % square vertices for non tilted 
        X = [-h  h  h -h];
        Y = [-h -h  h  h];
        
        % to tilt the squares
        theta = deg2rad(angle);
        R = [cos(theta) -sin(theta); 
             sin(theta)  cos(theta)];
        
        % Rotation
        coords = R * [X; Y];
        Xr = coords(1,:) + xc;
        Yr = coords(2,:) + yc;
        
        plot([Xr Xr(1)], [Yr Yr(1)], 'w-', 'LineWidth', 2);
        %axis equal
    end

    scatter(centers(:,1), centers(:,2), 20, 'filled', 'MarkerFaceColor','white');
end

