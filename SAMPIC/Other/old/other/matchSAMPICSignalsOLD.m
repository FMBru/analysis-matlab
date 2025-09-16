
referenceDetectorChannel = 0; %channel number of detector used as timing reference - for latency relative to this detector

channelsEnabled=[];

matchedEvents = [];

matchingDUTREFWindowWidth = 5; %max offset in ns between signal in DUT and REF for matching
matchingEventCounterWindowWidth = matchingDUTREFWindowWidth; %max offset in ns between signal in trigger and REF for matching

eventCounterLatency = 4; %ns delay between event counter and REF

%determine enabled channels, create arrays for cut
for i = 1:length(signals)
    if length(channelsEnabled)==0 || any(ismember(channelsEnabled,signals(i).ch)) == 0
        channelsEnabled = [channelsEnabled;signals(i).ch];
    end
    signalChannels (i) = signals(i).ch;
end

for i = 1:length(channelsEnabled)
    channelID = channelsEnabled(i);
    signalSelectionCut = signalChannels == channelID;
    eval(['signals_ch' num2str(channelID) ' = signals(signalSelectionCut);']);
    
    if  channelsEnabled(i) == referenceDetectorChannel
       %is reference channel 
       signals_chRef = signals(signalSelectionCut);
    end
%    else
%       %is dut channel
%       eval(['ch_' num2str(channelsEnabled(i)) '=[];']);
%    end
end

for i = 1:length(signals_chRef)
   refTime (i) = signals_chRef(i).cell0Time;
end




for i = 1:length(channelsEnabled)
   if  channelsEnabled(i) == referenceDetectorChannel
      %is reference channel 
   else
      %is dut channel
      eval(['signalsDUTchannel = signals_ch' num2str(channelsEnabled(i)) ';']);
      timestamps = [];
      timeDiffToRefArray = [];
      for k = 100:length(signalsDUTchannel)
        %loop through all events in DUT to find time differences to
        %reference
        dutSignal = signalsDUTchannel(k);
        timeDifferencesDUTRef = [];
        detectorHitTime = dutSignal.cell0Time;
        detectorHitTimeUnix = dutSignal.unixTimestamp;
        timestampEntry = [detectorHitTime detectorHitTimeUnix];
        timestamps = [timestamps; timestampEntry];
        timeDiffToRef = detectorHitTime-refTime;
       % hist(abs(timeDiffToRef));
       %pause(1);
       % close all
        minTimeDiffToRef = min(abs(timeDiffToRef));
        if minTimeDiffToRef<100
        	timeDiffToRefArray = [timeDiffToRefArray;minTimeDiffToRef];
        end
      end
      
      %determine latency from min time diff histogram between ch and ref
      latency = median (timeDiffToRefArray);
      % latency = -3; %override to set latency manually
     
      figure
      hold on
      hist(timeDiffToRefArray,100);
      xline(latency,'color','red','linewidth',2)
      xlabel('Min time between DUT and REF (ns)');
      ylabel('Events');
      title(['SAMPIC - Minimum time between DUT and REF - CH: ' num2str(channelsEnabled(i)) ' - Run ' run.id]);
      grid on

      pause(1);
      
      
      % latency is determined
      % assembleEvents
      
      minTime = latency-matchingDUTREFWindowWidth;
      maxTime = latency+matchingDUTREFWindowWidth;
      
      matchedEventsTemp = [];
      
      %loop through events of this channel to find matches
      for k = 1:length(signalsDUTchannel)
        dutSignal = signalsDUTchannel(k);
        timeDifferencesDUTRef = [];
        detectorHitTime = dutSignal.cell0Time;
        timeDiff = refTime-detectorHitTime;
        
        timeMatchToRef = timeDiff>minTime & timeDiff<maxTime;

        matchIdx = find(timeMatchToRef==1);
        
        %if length(matchIdx==1) %unambiguous match
        if length(matchIdx==1)  %for debugging permit more and take first
           %found a single match in reference for this dut event -> unambiguois match
           event.dut = signalsDUTchannel(k);
           event.ref = signals_chRef(matchIdx(1));
           event.timeDiff = timeDiff;
           event.eventID = 0;
           matchedEventsTemp = [matchedEventsTemp;event];
        end
        
      end

      eval(['matchedEvents_ch' num2str(channelsEnabled(i)) ' = matchedEventsTemp;']);
      
   end
   
end

channelsEnabled










timeDifferencesDetectorScint = [];
timeDifferencesDetectorMCP = [];

timestamps = [];

hitNumberPadding = 10;

%numberHits = length(hits);

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