function VisualiseMultipadEvent (padCenters,channelsEnabled,dataArray,plotTitle,plotPath,plottingStringDesc,plottingString,minC,maxC,eventTrackerX,eventTrackerY,eventCalculatedX,eventCalculatedY,runID)


dataArray;
%close all
minC;
maxC;

padLength = 10; %mm
padHeight = 10; %mm

h1 = figure

centralCircleSize = 1;
coneCircleSize = 6;

% colorMap = [
%     1.0000    1.0000         0;
%     1.0000    0.9688    0.0313;
%     1.0000    0.9375    0.0625;
%     1.0000    0.9063    0.0938;
%     1.0000    0.8750    0.1250;
%     1.0000    0.8438    0.1563;
%     1.0000    0.8125    0.1875;
%     1.0000    0.7813    0.2188;
%     1.0000    0.7500    0.2500;
%     1.0000    0.7188    0.2813;
%     1.0000    0.6875    0.3125;
%     1.0000    0.6563    0.3438;
%     1.0000    0.6250    0.3750;
%     1.0000    0.5938    0.4063;
%     1.0000    0.5625    0.4375;
%     1.0000    0.5313    0.4688;
%     1.0000    0.5000    0.5000;
%     1.0000    0.4688    0.5313;
%     1.0000    0.4375    0.5625;
%     1.0000    0.4063    0.5938;
%     1.0000    0.3750    0.6250;
%     1.0000    0.3438    0.6563;
%     1.0000    0.3125    0.6875;
%     1.0000    0.2813    0.7188;
%     1.0000    0.2500    0.7500;
%     1.0000    0.2188    0.7813;
%     1.0000    0.1875    0.8125;
%     1.0000    0.1563    0.8438;
%     1.0000    0.1250    0.8750;
%     1.0000    0.0938    0.9063;
%     1.0000    0.0625    0.9375;
%     1.0000    0.0313    0.9688;
%     1.0000         0    1.0000;;
%     1.0000    0.0323    0.9677;
%     1.0000    0.0645    0.9355;;
%     1.0000    0.0968    0.9032;
%     1.0000    0.1290    0.8710;
%     1.0000    0.1613    0.8387;
%     1.0000    0.1935    0.8065;
%     1.0000    0.2258    0.7742;
%     1.0000    0.2581    0.7419;
%     1.0000    0.2903    0.7097;
%     1.0000    0.3226    0.6774;
%     1.0000    0.3548    0.6452;
%     1.0000    0.3871    0.6129;
%     1.0000    0.4194    0.5806;
%     1.0000    0.4516    0.5484;
%     1.0000    0.4839    0.5161;
%     1.0000    0.5161    0.4839;
%     1.0000    0.5484    0.4516;
%     1.0000    0.5806    0.4194;
%     1.0000    0.6129    0.3871;
%     1.0000    0.6452    0.3548;
%     1.0000    0.6774    0.3226;
%     1.0000    0.7097    0.2903;
%     1.0000    0.7419    0.2581;
%     1.0000    0.7742    0.2258;
%     1.0000    0.8065    0.1935;
%     1.0000    0.8387    0.1613;
%     1.0000    0.8710    0.1290;
%     1.0000    0.9032    0.0968;
%     1.0000    0.9355    0.0645;
%     1.0000    0.9677    0.0323;
%     1.0000    1.0000         0];
%
%

colorMap = colormap;
maxDataPoint = max(dataArray)

for r=1:10
    for c=1:10
        padX = (r-1)*padLength;
        padY = (c-1)*padHeight;
        
        padNumber = (10-r+1)+(10-c)*10;
        
        %check if value for hitmap exists
        
        channelID = getChannelForPadNumber(padNumber,runID);
        padNumber;
        channelID;
        
        foundChannel = false;
        for pos = 1:length(channelsEnabled)
            if channelsEnabled(pos) == channelID
                countsChannel = dataArray(pos);
                
                if minC==0 && maxC == 0
                    
                    relativePosition = countsChannel/maxDataPoint;
                else
                    relativePosition = (countsChannel-minC)/(maxC-minC);
                end
                colorMapIndex = round(relativePosition*length(colorMap));
                
                if isnan(colorMapIndex) || colorMapIndex==0 || colorMapIndex<0
                    
                else
                    if colorMapIndex>length(colorMap)
                        colorMapIndex = length(colorMap);
                    end
                    
                    colorMapEntry = colorMap(colorMapIndex,:);
                    
                    rectangle('Position',[(padX-1) (padY-2) padLength padHeight],'FaceColor',[colorMapEntry(1) colorMapEntry(2) colorMapEntry(3)]);
                    
                    str=sprintf(plottingString, countsChannel);
                    text(padX,(padY+2),str,'FontSize',8);
                    
                end
            end
        end
        str = ['' num2str(padNumber) ' - ' plottingStringDesc];
        text(padX,padY,str,'FontSize',6);
        
        
        
        
        rectangle('Position',[(padX-1) (padY-2) padLength padHeight]);
    end
end

%rectangle('Position', [eventTrackerX eventTrackerY 2 2], 'Curvature', [1, 1], 'FaceColor', [0, 1, 0]);
%rectangle('Position', [eventCalculatedX eventCalculatedY 2 2], 'Curvature', [1, 1], 'FaceColor', [1, 0, 0]);



rectangle('Position',[eventTrackerX-centralCircleSize/2 eventTrackerY-centralCircleSize/2 centralCircleSize centralCircleSize],'Curvature',[1,1],'FaceColor', [0, 1, 0, 0.8],'EdgeColor', [0, 1, 0, 0.7]);
rectangle('Position',[eventTrackerX-coneCircleSize/2 eventTrackerY-coneCircleSize/2 coneCircleSize coneCircleSize],'Curvature',[1,1],'FaceColor', [0, 0, 1, 0.3],'EdgeColor', [0, 0, 1, 0.7]);

rectangle('Position',[eventCalculatedX-centralCircleSize/2 eventCalculatedY-centralCircleSize/2 centralCircleSize centralCircleSize],'Curvature',[1,1],'FaceColor', [1, 0, 0, 0.8],'EdgeColor', [1, 0, 0, 0.7]);


xlim([-10 110]);
ylim([-10 110]);

%caxis([minC maxC]);
xlabel('X position (mm)');
ylabel('Y position (mm)');
title(plotTitle);

set(h1,'Position',[100 100 800 800])

saveas(gcf,[plotPath])

end
