close all;
shouldDisplayMaps = false
eventsArray = [];
residualsArrayX = [];
residualsArrayY = [];
eventTrackerXArray = [];
eventTrackerYArray = [];
multiplicityArray = [];
sumAmpArray = [];
check = 0;

%loop through event IDs to assemble events in structures
for eventPos = 1:length(uniqueEventIDs)
    eventTemp.eventID = uniqueEventIDs(eventPos);
    
    shouldPlot = false;
    if mod(eventPos,plotInterval)==0
        shouldPlot = true;
    end

    %search for all events with this ID
    eventIDFilter = eventIDArray==eventTemp.eventID & glbl_cut; %& eventIDArray>10000;
        str_disp=sprintf('Plotting event ID %d', eventTemp.eventID);
        disp(str_disp);

    chIDArrayEvent = padIDArray(eventIDFilter);
    
    if length(chIDArrayEvent)>0
        %valid event
        time_diff_sigmoidEvent = time_diff_sigmoid(eventIDFilter);
        MCP_maxyEvent = MCP_maxy(eventIDFilter);
        MM_maxyEvent = MM_maxy(eventIDFilter);

        

        trackerXEvent = trackerX(eventIDFilter);
        trackerYEvent = trackerY(eventIDFilter);
        MM_DataArrayEvent = MM_DataArray(eventIDFilter);

        %invert X relative to center
        trackerXEvent = -(trackerXEvent-pad0Pos(1))+pad0Pos(1);

%        MCP_DataArrayEvent = MCP_DataArray(eventIDFilter);
        %acceptedSignalArrayEvent = acceptedSignalArray(eventIDFilter);
        numberAcceptedSignalsInEvent = length(MCP_maxyEvent)
        MM_maxySum = sum(MM_maxyEvent);

        if max(MM_maxySum)>0.050 && numberAcceptedSignalsInEvent>1 %show only events with strong signal
                size(MM_DataArrayEvent);
                %VisualisePicolargeEvent(padIDArrayEvent,1000*MM_maxyEvent,['Picolarge event ' num2str(eventTemp.eventID)],[store_folder '\Event' int2str(eventTemp.eventID) '.png'],'Amp','%0.1f mV',minC,maxC,trackerXEvent(1),trackerYEvent(1))
                waveformsArray = [];

                positionRecoX = 0;
                positionRecoY = 0;
                chIDArrayEvent = unique(chIDArrayEvent);
                %loop through hit pads assemble waveforms array
                numChPresent = length(chIDArrayEvent);
                %padIDArrayEvent
                
                chIDArrayEvent;
                padIDArrayEvent = chIDArrayEvent;

                for pos = 1:length(chIDArrayEvent)
                    ampWeight = MM_maxyEvent(pos)/MM_maxySum;
                    chNumber = chIDArrayEvent(pos); %+1;
                    
                    padIndex = getPadForChannelNumberMediumGranularity(chNumber);
                    padIDArrayEvent (pos) = padIndex;

                     if padIndex <=0 
                         check = check + 1;
                         break
                     end
                    padIndex
                    posPadX = padCenters(padIndex,1)
                    posPadY = padCenters(padIndex,2)
                    
                    %add up estimate from weights
                    positionRecoX = positionRecoX+posPadX*ampWeight;
                    positionRecoY = positionRecoY+posPadY*ampWeight;
                    
                    waveformTemp.padID = padIDArrayEvent(pos);
                    waveformTemp.x = MM_DataArrayEvent(pos).waveform(:,1) - MM_DataArrayEvent(pos).sigmoid.timepoint;
                    waveformTemp.riseTime=MM_DataArrayEvent(pos).sigmoid.timepoint90-MM_DataArrayEvent(pos).sigmoid.timepoint10;% extract electrom peak charge
                    
                    waveformTemp.y = MM_DataArrayEvent(pos).waveform(:,2);
                    %waveformTemp.xRef = MCP_DataArrayEvent(pos).sig.xfit;
                    %waveformTemp.yRef = MCP_DataArrayEvent(pos).sig.yfit;
                    waveformsArray = [waveformsArray waveformTemp];
                end
                positionRecoX
                positionRecoY
                chIDArrayEvent

                residualX = positionRecoX-trackerXEvent(1);
                residualY = positionRecoY-trackerYEvent(1);
                residualsArrayX = [residualsArrayX;residualX];
                residualsArrayY = [residualsArrayY;residualY];
                eventTrackerXArray = [eventTrackerXArray;trackerXEvent(1)];
                eventTrackerYArray = [eventTrackerYArray;trackerYEvent(1)];
                multiplicityArray = [multiplicityArray;length(chIDArrayEvent)];
                sumAmpArray = [sumAmpArray;MM_maxySum];

                if shouldDisplayMaps %& saveMapCounter<numberMapsToSave
                    %VisualisePicolargeEventWaveforms(padIDArrayEvent,1000*MM_maxyEvent,waveformsArray,['Picolarge event ' num2str(eventTemp.eventID)],[store_folder '\Event' int2str(eventTemp.eventID) '.png'],'Amp','%0.1f mV',minC,maxC,trackerXEvent(1),trackerYEvent(1))
                    VisualiseMediumGranularityEventWaveformsPosReco(padCenters,padIDArrayEvent,1000*MM_maxyEvent,waveformsArray,[runTitleString ' - Event ' num2str(eventTemp.eventID)],[store_folderEvents '\Event' int2str(eventTemp.eventID) '.png'],'','%0.1f mV',minC,maxC,trackerXEvent(1),trackerYEvent(1),positionRecoX,positionRecoY)
                    hold off;
                    pause(1);
                    close all;
    
                    saveMapCounter = saveMapCounter+1;
                end
        end
    end
end