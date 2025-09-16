function [out] = resolutionStuff(finalTable,geom,folder, quantity, name, unit, addInfo)
%RESOLUTIONSTUFF Summary of this function goes here

[frame, subfinalTable] = geometricalFrame(finalTable, geom, 1, 'center');

% Compute counts map
[~,~,countMap] = twoDimensionsMap(frame, ones(height(subfinalTable),1));
plotMap(frame, countMap);
title(['Count map full' addInfo] );
saveas(gcf,[folder '\countMap-' addInfo '.png']);
close;

[MeanMap, RMSMap] = twoDimensionsMap(frame, subfinalTable.(quantity));

%Plot map
plotMap(frame, MeanMap);
title(['Mean ' name ' map']);
hold on
drawRectangles(geom.centers,9.8,0.4);
hold off
saveas(gcf,[folder '\mean-' name '-Map-' addInfo '.png']);
close

% Plot the time res map
plotMap(frame, RMSMap);
title(['RMS ' name ' map']);
hold on
drawRectangles(geom.centers,9.8,0.4);
hold off
saveas(gcf,[folder '\RMS-' name '-Map-' addInfo '.png']);
close

distribution(RMSMap(:),20,['RMS ' name], unit, folder, addInfo, true);
distribution(subfinalTable.(quantity),100,name,unit, folder, addInfo, true);

end

