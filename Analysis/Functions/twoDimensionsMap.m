function [map, stdMap, sumMap] = twoDimensionsMap(geom_frame, values)
%2DMAP: Creates a 2D map using accumarray (mean values, sum and std) 

% Keep only events inside the bin range
valid = geom_frame.xBin > 0 & geom_frame.yBin > 0;

% Compute mean in each bin
map = accumarray([geom_frame.yBin(valid), geom_frame.xBin(valid)], values(valid), [numel(geom_frame.yEdges)-1, numel(geom_frame.xEdges)-1], @(x) mean(x, 'omitnan'), NaN);

% Compute std in each bin
stdMap = accumarray([geom_frame.yBin(valid), geom_frame.xBin(valid)], values(valid), [numel(geom_frame.yEdges)-1, numel(geom_frame.xEdges)-1], @(x) std(x, 'omitnan'), NaN);

% Compute sum in each bin
sumMap = accumarray([geom_frame.yBin(valid), geom_frame.xBin(valid)], values(valid), [numel(geom_frame.yEdges)-1, numel(geom_frame.xEdges)-1], @(x) sum(x, 'omitnan'), NaN);



end

