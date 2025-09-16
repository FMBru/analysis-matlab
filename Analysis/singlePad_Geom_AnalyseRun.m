%clear all
close all

%% configuration
shouldSaveMat = false;
twalk.en = 1;                      % enable timewalk correction
enable_cut_twalk = true;           %plot the timewalk correction in the sampling area

x_window_length = 30;               %x size of the region in which we are focusing geometrically
y_window_length = 30;               %y size of the region in which we are focusing geometrically

pad_window_length = 10;             %pad size (multipad)

%% table approach
all_events_table = table(trackerX(:), trackerY(:), MM_maxy(:), MCP_maxy(:), e_peak_MM(:), time_diff_sigmoid(:), riseTime(:), 'VariableNames', {'X', 'Y', 'MM_amp', 'MCP_amp', 'e_peak_MM', 'SAT', 'riseTime'});

%create 2D framework (to be improved)
x_first = median(trackerX) - x_window_length/2;
x_last = median(trackerX) + x_window_length/2;
x_step = 1;

y_first = median(trackerY) - y_window_length/2;
y_last = median(trackerY) + y_window_length/2;
y_step = 1;

xEdges = x_first:x_step:x_last;
yEdges = y_first:y_step:y_last;

xCenters = xEdges(1:end-1) + diff(xEdges)/2;
yCenters = yEdges(1:end-1) + diff(yEdges)/2;

[Xgrid, Ygrid] = meshgrid(xCenters, yCenters);

%packing this geometrical frame into a structure
geom_frame.xEdges = xEdges;
geom_frame.yEdges = yEdges;
geom_frame.xCenters = xCenters;
geom_frame.yCenters = yCenters;
geom_frame.Xgrid = Xgrid;
geom_frame.Ygrid = Ygrid;

%select the subset of events (cuts)
MM_ampCut = all_events_table.MM_amp > 0.01*max(MM_maxy) & all_events_table.MM_amp < 0.99*max(MM_maxy);
MCP_ampCut = all_events_table.MCP_amp >= 0.15 & all_events_table.MCP_amp > 0.01*max(MCP_maxy) & all_events_table.MCP_amp < 0.99*max(MCP_maxy);
riseTimeCut = abs(all_events_table.riseTime) < 10;
SATCut = abs(all_events_table.SAT - median(time_diff_sigmoid)) < 5;
squareCut = all_events_table.X > x_first & all_events_table.X < x_last & all_events_table.Y > y_first & all_events_table.Y < y_last;

subset = all_events_table(squareCut & MM_ampCut & MCP_ampCut & riseTimeCut & SATCut, :);

% Assign each event to a bin
[~, ~, xBin] = histcounts(subset.X, xEdges);
[~, ~, yBin] = histcounts(subset.Y, yEdges);

% Keep only events inside the bin range
valid = xBin > 0 & yBin > 0;

% Compute counts map
countMap = accumarray([yBin(valid), xBin(valid)], 1, [numel(yEdges)-1, numel(xEdges)-1], @sum, 0);
plotMap(xEdges, yEdges, countMap);
title('Count map after cuts');

% Compute amplitude maps (only mean amplitude map on the plot)
[ampMap, stdampMap, sumampMap] = twoDimensionsMap(xBin, yBin, xEdges, yEdges, subset.MM_amp);
title('Mean Amplitude Map');


%% find the center of the pad

[pad_x_center,pad_y_center] = findPadCenter(geom_frame, sumampMap, pad_window_length);


%% other maps in the big frame

% Compute mean SAT map
[satMap, stdsatMap, sumsatMap] = twoDimensionsMap(xBin, yBin, xEdges, yEdges, subset.SAT);
title('Mean SAT Map');

% Compute mean e-peak in each bin
e_peak_Map = twoDimensionsMap(xBin, yBin, xEdges, yEdges, subset.e_peak_MM);
title('Mean e-peak Map');

% Compute mean rise time in each bin
riseTime_Map = twoDimensionsMap(xBin, yBin, xEdges, yEdges, subset.riseTime);
title('Mean rise time Map');

% mask = Xgrid <= 75 & Xgrid >= 55 & Ygrid <= 68 & Ygrid >= 55;
% ampMapCut = countMap;
% ampMapCut(~mask) = NaN;
% 
% plotMap(xEdges, yEdges, ampMapCut);

%% selecting a smaller region

% override pad center
% pad_x_center = 75.2;
% pad_y_center = 79.3;

%create 2D framework (to be improved)
x_first = pad_x_center - pad_window_length/2;
x_last = pad_x_center + pad_window_length/2;
x_step = 0.25;

y_first = pad_y_center - pad_window_length/2;
y_last = pad_y_center + pad_window_length/2;
y_step = 0.25;

xEdges = x_first:x_step:x_last;
yEdges = y_first:y_step:y_last;

xCenters = xEdges(1:end-1) + diff(xEdges)/2;
yCenters = yEdges(1:end-1) + diff(yEdges)/2;

[Xgrid, Ygrid] = meshgrid(xCenters, yCenters);

%packing this geometrical frame into a structure
geom_frame.xEdges = xEdges;
geom_frame.yEdges = yEdges;
geom_frame.xCenters = xCenters;
geom_frame.yCenters = yCenters;
geom_frame.Xgrid = Xgrid;
geom_frame.Ygrid = Ygrid;

%select the subset of events (geometrical cuts)
squareCut = subset.X > x_first & subset.X < x_last & subset.Y > y_first & subset.Y < y_last;

pad_table = subset(squareCut, :);

% Assign each event to a bin
[~, ~, xBin] = histcounts(pad_table.X, xEdges);
[~, ~, yBin] = histcounts(pad_table.Y, yEdges);

% Keep only events inside the bin range
valid = xBin > 0 & yBin > 0;

% Compute counts map
countMap = accumarray([yBin(valid), xBin(valid)], 1, [numel(yEdges)-1, numel(xEdges)-1], @sum, 0);
plotMap(xEdges, yEdges, countMap);
title('Count map after cuts');

% Compute amplitude maps (only mean amplitude map on the plot)
[ampMap, stdampMap, sumampMap] = twoDimensionsMap(xBin, yBin, xEdges, yEdges, pad_table.MM_amp);
title('Mean Amplitude Map');


%% other maps in the big frame

% Compute mean SAT map
[satMap, stdsatMap, sumsatMap] = twoDimensionsMap(xBin, yBin, xEdges, yEdges, pad_table.SAT);
title('Mean SAT Map');

% Compute mean e-peak in each bin
e_peak_Map = twoDimensionsMap(xBin, yBin, xEdges, yEdges, pad_table.e_peak_MM);
title('Mean e-peak Map');

% Compute mean rise time in each bin
riseTime_Map = twoDimensionsMap(xBin, yBin, xEdges, yEdges, pad_table.riseTime);
title('Mean rise time Map');
