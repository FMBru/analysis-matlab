function [geom_frame, table_subset] = geometricalFrame(table,geom, step, opt)
%GEOMETRICALFRAME: creates a 2D geometrical framework around the median of
%the X and Y in the TABLE, dividing the space in a grid with STEP pitch
%returns then the geometrical frame with the edges of the bin and the grid
%and the subset of the events of the table inside the frame

x_centroid = median(table.X);
y_centroid = median(table.Y);

if nargin < 4 || isempty(opt)
    opt = 'median';
end
 
if strcmp(opt, 'center') && isfield(geom, 'x_centroid') && isfield(geom, 'y_centroid') 
    x_centroid = geom.x_centroid;
    y_centroid = geom.y_centroid;
end


%create 2D framework (to be improved)
x_first = x_centroid - geom.x_window_length/2;
x_last = x_centroid + geom.x_window_length/2;
x_step = step;

y_first = y_centroid - geom.y_window_length/2;
y_last = y_centroid + geom.y_window_length/2;
y_step = step;

xEdges = x_first:x_step:x_last;
yEdges = y_first:y_step:y_last;

xCenters = xEdges(1:end-1) + diff(xEdges)/2;
yCenters = yEdges(1:end-1) + diff(yEdges)/2;

[Xgrid, Ygrid] = meshgrid(xCenters, yCenters);

geoCut = table.X > x_first & table.X < x_last & table.Y > y_first & table.Y < y_last;
table_subset = table(geoCut, :);

% Assign each event to a bin
[~, ~, xBin] = histcounts(table_subset.X, xEdges);
[~, ~, yBin] = histcounts(table_subset.Y, yEdges);



%packing this geometrical frame into a structure
geom_frame.xEdges = xEdges;
geom_frame.yEdges = yEdges;
geom_frame.xCenters = xCenters;
geom_frame.yCenters = yCenters;
geom_frame.Xgrid = Xgrid;
geom_frame.Ygrid = Ygrid;
geom_frame.xBin = xBin;
geom_frame.yBin = yBin;
geom_frame.geoCut = geoCut;

end

