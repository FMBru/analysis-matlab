function [pad_x_center, pad_y_center] = singlePadMaps(padNumber, sub_frame, subset, geom, store_folder, shouldSave, opt)
%SINGLEPADPLOTS make plots for count map, amplitude map, rise time map, SAT
%map, e-peak map and find the center
if nargin < 7 || isempty(opt)
    opt = 'pad';
end

vars = {'MM_amp', 'SAT', 'e_peak_MM', 'riseTime'};

padNumberString = num2str(padNumber);

if strcmp(opt, 'full')
    for i=1:length(vars)
        vars{i} = [vars{i} padNumberString];
    end
end



% Compute counts map
[~,~,countMap] = twoDimensionsMap(sub_frame, ones(height(subset),1));
plotMap(sub_frame, countMap);
title(['Count map pad' padNumberString]);
if shouldSave
    saveas(gcf,[store_folder '\CountMap_pad' padNumberString '.png'])
end
close

% Compute amplitude maps 
[ampMap, stdampMap, sumampMap] = twoDimensionsMap(sub_frame, subset.(vars{1}));
% Plot the mean map
plotMap(sub_frame, ampMap);
title(['Mean amp map pad' padNumberString]);
if strcmp(opt, 'full')
    hold on
    drawRectangles(geom.centers,9.8,0);
    hold off
end
if shouldSave
    saveas(gcf,[store_folder '\meanAmpMap_pad' padNumberString '.png']);
end
close

% Compute amplitude maps 
ampMCPMap = twoDimensionsMap(sub_frame, subset.MCP_amp);
% Plot the mean map
plotMap(sub_frame, ampMCPMap);
title(['Mean MCP amp map pad' padNumberString]);
if shouldSave
    saveas(gcf,[store_folder '\meanMCPAmpMap_pad' padNumberString '.png']);
end
close

% Compute mean SAT map
[satMap, stdsatMap, sumsatMap] = twoDimensionsMap(sub_frame, subset.(vars{2}));
% Plot the mean map
plotMap(sub_frame, satMap);
title(['Mean SAT map pad' padNumberString]);
if strcmp(opt, 'full')
    hold on
    drawRectangles(geom.centers,9.8,0);
    hold off
end
if shouldSave
    saveas(gcf,[store_folder '\meanSATMap_pad' padNumberString '.png']);
end
close

% Plot the time res map
plotMap(sub_frame, stdsatMap);
title(['Mean time res map pad' padNumberString]);
if strcmp(opt, 'full')
    hold on
    drawRectangles(geom.centers,9.8,0);
    hold off
end
if shouldSave
    saveas(gcf,[store_folder '\timeresMap_pad' padNumberString '.png']);
end
close

% Compute mean e-peak in each bin
e_peak_Map = twoDimensionsMap(sub_frame, subset.(vars{3}));
[pad_x_center,pad_y_center] = findPadCenter(sub_frame, e_peak_Map, geom.pad_window_length, 2.5);
close all;
plotMap(sub_frame,e_peak_Map);
hold on
scatter(pad_x_center, pad_y_center, 20, 'filled', 'MarkerFaceColor','white');
hold off
title(['Mean e-peak map pad' padNumberString]);
if shouldSave
    saveas(gcf,[store_folder '\meanEPEAKMap_pad' padNumberString '.png']);
end
close

% Compute mean rise time in each bin
riseTime_Map = twoDimensionsMap(sub_frame, subset.(vars{4}));
plotMap(sub_frame,riseTime_Map);
title(['Rise time map pad' padNumberString]);
if shouldSave
    saveas(gcf,[store_folder '\meanRiseTimeMap_pad' padNumberString '.png']);
end
close

end

