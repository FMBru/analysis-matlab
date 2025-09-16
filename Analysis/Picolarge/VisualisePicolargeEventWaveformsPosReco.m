function VisualisePicolargeEventWaveformsPosReco(padIDsArray,dataArray,waveformsArray,plotTitle,plotPath,plottingStringDesc,plottingString,minC,maxC,trackerX,trackerY,recoX,recoY)
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
                    %subplot(5,5,[1 2 3 4 6 7 8 9 11 12 13 14 16 17 18 19])
                    subplot(5,5,[7 8 9 12 13 14 17 18 19])
                    %rectangle('Position',[(padX-1) (padY-2) padLength padHeight],'FaceColor',[colorMapEntry(1) colorMapEntry(2) colorMapEntry(3)]);
                    hexagon(5,padX,padY,[colorMapEntry(1) colorMapEntry(2) colorMapEntry(3)])
                    hold on

                    str=sprintf(plottingString, countsChannel);
                    text(padX-2,(padY+2),str,'FontSize',10);
                end
            end
        end
        str = ['Pad' num2str(padNumber) '' plottingStringDesc];
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

rectangle('Position',[trackerX-centralCircleSize/2 trackerY-centralCircleSize/2 centralCircleSize centralCircleSize],'Curvature',[1,1],'FaceColor', [0, 1, 0, 0.8],'EdgeColor', [0, 1, 0, 0.7]);
rectangle('Position',[trackerX-coneCircleSize/2 trackerY-coneCircleSize/2 coneCircleSize coneCircleSize],'Curvature',[1,1],'FaceColor', [0, 0, 1, 0.3],'EdgeColor', [0, 0, 1, 0.7]);

rectangle('Position',[recoX-centralCircleSize/2 recoY-centralCircleSize/2 centralCircleSize centralCircleSize],'Curvature',[1,1],'FaceColor', [1, 0, 0, 0.8],'EdgeColor', [1, 0, 0, 0.7]);

text(5,50, '• Hit position tracker', 'Color', 'g','FontSize',12)
text(5, 47, '• Reconstructed - centroid active pads', 'Color', 'r','FontSize',12)
c = colorbar('eastoutside');
%c.Label.String = 'Signal amplitude (mV)';
% axis equal
                    subplotPosFree = [3 1 11 21 25 15 5];

%% plot waveforms
                for pos = 1:length(waveformsArray)
                       waveformTemp = waveformsArray(pos);
                       subplot(5,5,subplotPosFree(1+waveformTemp.padID))
                        plot(waveformTemp.x,1000*waveformTemp.y);
                        ylim([0 200]);
                        xlim([-1 2]);
                        %xlabel('Time (ns)');
                        %ylabel('Amplitude (mV)');

                        grid on
                        title(['Pad ' int2str(waveformTemp.padID)])  
                        risetimeString = sprintf('Risetime %.1fns',waveformTemp.riseTime);
                        %text(-0.8,180,risetimeString);
                end

%                 subplot(5,5,23)
%                 for pos = 1:length(waveformsArray)
%                        waveformTemp = waveformsArray(pos);
%                         plot(waveformTemp.x,waveformTemp.y); hold on
%                 end
%                  xlim([0 0.1]);
%                title(['All pads '])  


grid on
set(h1,'Position',[200 200 900 900])

saveas(gcf,[plotPath])

end
