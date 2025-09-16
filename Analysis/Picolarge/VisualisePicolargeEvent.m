function VisualisePicolargeEvent(padIDsArray,dataArray,plotTitle,plotPath,plottingStringDesc,plottingString,minC,maxC,trackerX,trackerY)
%close all

padSide = 10;

h1 = figure;

colorMap = colormap;
maxDataPoint = max(dataArray);

pad0Pos = [29 27]; %center of pad 0 determiend in other run aligned with pad0

padCenters= 0.5*[
    0   0;
    -8.66 15;
    -17.32 0;
    -8.66 -15;
    8.66 -15;
    17.32 0;
    8.66 15;
];



for r=1:7 %loop through pads
            padNumber = r-1;

        padCenter = padCenters(r,:);
        padX = padCenter(1) + pad0Pos(1);
        padY = padCenter(2) + pad0Pos(2);
        
        %check if value for hitmap exists
        channelID = padNumber;
        padNumber;
        channelID;
        
                    hexagonOutline(5,padX,padY);            
                    hold on
        
        foundChannel = false;
        for pos = 1:length(padIDsArray)
            
            if padIDsArray(pos) == channelID
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
                    
                    %rectangle('Position',[(padX-1) (padY-2) padLength padHeight],'FaceColor',[colorMapEntry(1) colorMapEntry(2) colorMapEntry(3)]);
                    hexagon(5,padX,padY,[colorMapEntry(1) colorMapEntry(2) colorMapEntry(3)])
                    hold on

                    str=sprintf(plottingString, countsChannel);
                    text(padX-2,(padY+2),str,'FontSize',10);
                end
            end
        end
        str = ['Pad' num2str(padNumber) ' - ' plottingStringDesc];
        text(padX-2,padY,str,'FontSize',10);
end

plotRange = 25;

xlim([pad0Pos(1)-plotRange pad0Pos(1)+plotRange]);
ylim([pad0Pos(2)-plotRange pad0Pos(2)+plotRange]);
caxis([minC maxC]);

xlabel('X position (mm)');
ylabel('Y position (mm)');
title(plotTitle);

centralCircleSize = 1;
coneCircleSize = 6;

rectangle('Position',[trackerX-centralCircleSize/2 trackerY-centralCircleSize/2 centralCircleSize centralCircleSize],'Curvature',[1,1],'FaceColor', [1, 0, 0, 0.8],'EdgeColor', [1, 0, 0, 0.7]);
rectangle('Position',[trackerX-coneCircleSize/2 trackerY-coneCircleSize/2 coneCircleSize coneCircleSize],'Curvature',[1,1],'FaceColor', [0, 0, 1, 0.3],'EdgeColor', [1, 0, 0, 0.7]);

set(h1,'Position',[400 400 700 700])

saveas(gcf,[plotPath])

end
