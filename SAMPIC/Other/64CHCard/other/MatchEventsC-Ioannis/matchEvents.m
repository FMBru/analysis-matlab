timeDifferencesDetectorScint = [];
timeDifferencesDetectorMCP = [];

timestamps = [];

hitNumberPadding = 10;

numberHits = length(hits);

matchedEvents = [];

for i = 1:numberHits
        if mod(i,1000)==0
        i
    end

    %write file
    
    hit = hits(i);
    if hit.ch == 1
        %detector signal
        
        detectorHitTime = hit.cell0Time;
        detectorHitTimeUnix = hit.unixTimestamp;
        timestampEntry = [detectorHitTime detectorHitTimeUnix];
        timestamps = [timestamps; timestampEntry];
        
        smallTriggerHitsAllTimeDiff = [];
        smallTriggerHitsAllTimeDiffMCP = [];
        
        %is detector signal, go through others to find match
        limitEventIndexLow = i-hitNumberPadding;
        limitEventIndexHigh = i+hitNumberPadding;
        if limitEventIndexLow<1
            limitEventIndexLow = 1;
        end
        if limitEventIndexHigh>numberHits
            limitEventIndexHigh = numberHits;
        end
        
        
        for j = limitEventIndexLow:limitEventIndexHigh
            testHit = hits(j);
            if testHit.ch==2 && i ~= j
                testHitTimestamp = testHit.cell0Time;
                testHitTimestampUnix = testHit.unixTimestamp;
                timeDiff = (testHitTimestamp-detectorHitTime);
                diffEntry = [j timeDiff];
                smallTriggerHitsAllTimeDiff = [smallTriggerHitsAllTimeDiff;diffEntry];
                
            end
        end
        
        eventIndexMinDiff = 0;
        eventIndexMinDiffValue = -1;
        for j = 1:size(smallTriggerHitsAllTimeDiff,1)
            if eventIndexMinDiffValue<0
                eventIndexMinDiff = smallTriggerHitsAllTimeDiff(j,1);
                eventIndexMinDiffValue = smallTriggerHitsAllTimeDiff(j,2);
            elseif smallTriggerHitsAllTimeDiff(j,2)<eventIndexMinDiffValue
                eventIndexMinDiff = smallTriggerHitsAllTimeDiff(j,1);
                eventIndexMinDiffValue = smallTriggerHitsAllTimeDiff(j,2);
            end
        end
        
        smallestTimeDiff = eventIndexMinDiffValue;
        
        
        %if smallestTimeDiff<1e15
            timeDifferencesDetectorScint = [timeDifferencesDetectorScint;smallestTimeDiff];
       % end
        
        for k = limitEventIndexLow:limitEventIndexHigh
            testHitMCP = hits(k);
            if testHitMCP.ch==0
                testHitTimestampMCP = testHitMCP.cell0Time;
                testHitTimestampUnixMCP = testHitMCP.unixTimestamp;
                timeDiffMCP = (testHitTimestampMCP-detectorHitTime);
                diffEntryMCP = [k timeDiffMCP];

                smallTriggerHitsAllTimeDiffMCP = [smallTriggerHitsAllTimeDiffMCP;diffEntryMCP];
                
            end
        end
        
        eventIndexMinDiffMCP = 0;
        eventIndexMinDiffValueMCP = -1;
        for k = 1:size(smallTriggerHitsAllTimeDiffMCP,1)
            if eventIndexMinDiffValueMCP<0
                eventIndexMinDiffMCP = smallTriggerHitsAllTimeDiffMCP(k,1);
                eventIndexMinDiffValueMCP = smallTriggerHitsAllTimeDiffMCP(k,2);
            elseif smallTriggerHitsAllTimeDiffMCP(k,2)<eventIndexMinDiffValueMCP
                eventIndexMinDiffMCP = smallTriggerHitsAllTimeDiffMCP(k,1);
                eventIndexMinDiffValueMCP = smallTriggerHitsAllTimeDiffMCP(k,2);
            end
        end
        
        smallestTimeDiffMCP = eventIndexMinDiffValueMCP;
        
        
        %if smallestTimeDiffMCP<1e15
            timeDifferencesDetectorMCP = [timeDifferencesDetectorMCP;smallestTimeDiffMCP];
       % end
        
        %build event
        
        if eventIndexMinDiffMCP>0 %&& eventIndexMinDiff>0 %require also small trig event to be present?
        
            signalEvent = hit;
            mcpEvent = hits(eventIndexMinDiffMCP);

            event.signal = signalEvent;
            event.mcp = mcpEvent;
            if eventIndexMinDiff == 0
                 event.scint = mcpEvent; %no scintillator data present
            else
                scintEvent = hits(eventIndexMinDiff);

                event.scint = scintEvent;
            end

            matchedEvents = [matchedEvents;event];
        end
    end
end