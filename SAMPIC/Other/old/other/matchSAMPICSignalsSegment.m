close all


%load all signal segments from folder
processingFolder = ['\\eosproject-smb\eos\project\p\picosec\testbeam\2022_May_h4\SAMPIC\Processing\Run' run.id]; 
signalsFolder = [processingFolder '\Signals']; 

infoFile = [processingFolder '\info.mat'];
load(infoFile);

refSignalsFile = [signalsFolder '\signalsRef.mat'];
load(refSignalsFile);

matchedEventsFolder = [processingFolder '\MatchedEvents'];
mkdir(matchedEventsFolder);


store_folder = ['\\eosproject-smb\eos\project\p\picosec\testbeam\2022_May_h4\Results\Run' run.id '-SAMPIC'];

channelsEnabled=[];

matchedEvents = [];

matchingDUTREFWindowWidth = 5; %max offset in ns between signal in DUT and REF for matching
matchingEventCounterWindowWidth = 10000000; %max offset in ns between signal in trigger and REF for matching

eventCounterLatency = 4; %ns delay between event counter and REF

maxNumberEventCounterForLatency = 1000;
maxNumberChannelForLatency = 1000;

%determine enabled channels, create arrays for cut
% for i = 1:length(signals)
%     if length(channelsEnabled)==0 || any(ismember(channelsEnabled,signals(i).ch)) == 0
%         channelsEnabled = [channelsEnabled;signals(i).ch];
%     end
%     signalChannels (i) = signals(i).ch;
% end
channelsEnabled = [0;1];


%sort signals into channels
%{
for i = 1:length(channelsEnabled)
    channelID = channelsEnabled(i);
    signalSelectionCut = signalChannels == channelID;
    eval(['signals_ch' num2str(channelID) ' = signals(signalSelectionCut);']);
    
    if  channelsEnabled(i) == referenceDetectorChannel
        %is reference channel
        signalsRef = signals(signalSelectionCut);
    end
    %    else
    %       %is dut channel
    %       eval(['ch_' num2str(channelsEnabled(i)) '=[];']);
    %    end
end
%}

for i = 1:length(signalsRef)
    refTime (i) = signalsRef(i).cell0Time;
end

figure
hold on
hist(refTime,10000);
xlabel('Timestamp (ns)');
ylabel('Events');
title(['SAMPIC - Reference Timestamp - Run ' run.id]);
grid on
saveas(gcf,[store_folder '\Run' run.id '_refTimestamps.png'])

pause(1);
medianLatencyErrorArray = [];

%%get latency for trigger and assemble trigger array
for i = 1:length(channelsEnabled)
    if  channelsEnabled(i) == referenceDetectorChannel
        %is reference channel
        %determine latency between reference and SRS counter
        
        scaledEventCounterTimes = eventCounterEntries(:,4);
        
        
        %eval(['signalsDUTchannel = signals_ch' num2str(channelsEnabled(i)) ';']);
        timestamps = [];
        timeDiffToRefArray = [];
        signalsRefTime = [];
        
        test = 1
        
        for pos = 1:length(signalsRef)
            signalsRefTime(pos) = signalsRef(pos).cell0Time;
        end
        
        % = signalsRefTime>9e10 & signalsRefTime<10e10;
        signalsRefCut = signalsRefTime;
        signalsRefCut = signalsRefTime;
        timeDiffToRefArray = [];
        
        
        refTimeSamplesForLatency = length(refTime);
        if refTimeSamplesForLatency>maxNumberEventCounterForLatency
            refTimeSamplesForLatency = maxNumberEventCounterForLatency;
        end
        
        for k = 1:refTimeSamplesForLatency
            %loop through all events in DUT to find time differences to
            %reference
            refSignal = refTime(k);
            timeDifferencesDUTRef = [];
            detectorHitTime = refSignal;
            %detectorHitTime = refSignal.cell0Time;
            timeDiffToRef = detectorHitTime-scaledEventCounterTimes;
            % hist(abs(timeDiffToRef));
            %pause(1);
            % close all
            
            %       figure
            % hold on
            %hist(timeDiffToRefArray,1000);
            %pause(1);
            %close all
            minTimeDiffToRef = min(abs(timeDiffToRef));
            %if minTimeDiffToRef<1000000000
            timeDiffToRefArray = [timeDiffToRefArray;minTimeDiffToRef];
            %end
        end
        
        %timeDiffCutLow = timeDiffToRefArray(timeDiffToRefArray<1e8);
        % timeDiffCutHigh = timeDiffToRefArray(timeDiffToRefArray>1e8);
        
        %numberLow = length(timeDiffCutLow)
        % numberHigh = length(timeDiffCutHigh)
        
        % ratio = numberLow/numberHigh
        
        timeDiffToRefArrayMask = timeDiffToRefArray<1e6;
        timeDiffToRefArrayCut = timeDiffToRefArray(timeDiffToRefArrayMask);
        
        %determine latency from min time diff histogram between ch and ref
        latency = median (timeDiffToRefArrayCut);
        eventCounterLatency = latency;
        % latency = -3; %override to set latency manually
        %latencyEntry = [scalingFactor latency numberLow numberHigh ratio];
        latencyEntry = [latency];
        medianLatencyErrorArray = [medianLatencyErrorArray;latencyEntry];
        figure
        hold on
        hist(timeDiffToRefArrayCut,1000);
        if isnan(latency)
            
        else
            xline(latency,'color','red','linewidth',2)
        end
        xlabel('Min time between SRS counter and REF (ns)');
        ylabel('Events');
        title(['SAMPIC - Minimum time between EventCounter REF - CH: ' num2str(channelsEnabled(i)) ' - Run ' run.id]);
        grid on
        %xlim([0 30000]);
        saveas(gcf,[store_folder '\Run' run.id '_latency_counter_scaling.png'])
        %close all
        pause(1);
        
        [refTimeCounts,refTimeCenters] = hist(signalsRefCut,10000);
        [eventCounterCounts, eventCounterCenters] = hist(scaledEventCounterTimes,10000);
        
        close all
        
        figure
        hold on
        plot(refTimeCenters,refTimeCounts);
        plot(eventCounterCenters,eventCounterCounts);
        if isnan(latency)
            
        else
            xline(latency,'color','red','linewidth',2)
        end
        xlabel('Time');
        ylabel('Events');
        title(['SAMPIC - Minimum time between EventCounter REF - CH: ' num2str(channelsEnabled(i)) ' - Run ' run.id]);
        legend('Ref hits time','Event Counter time')
        grid on
        
        saveas(gcf,[store_folder '\Run' run.id '_timeComparison_scaling.png'])
        %close all
        pause(1);
        
    end
end



        clear signals;



        signalSegmentPath = [signalsFolder '\segment1.mat'];
        load(signalSegmentPath);
        %signalChannels=zeros(length(signals),1);
        for channelPos = 1:length(signals)
            signalChannels (channelPos) = signals(channelPos).ch;
        end     


%% loop through DUT channels  - get latency
for i = 1:length(channelsEnabled)
    if  channelsEnabled(i) == referenceDetectorChannel
        
    else
        %is dut channel
        
        %load first segment of this channel for latency calculation
        signalsCut = signalChannels==channelsEnabled(i);
        signalsDUTchannel = signals(signalsCut);
        
       % eval(['signalsDUTchannel = signals_ch' num2str(channelsEnabled(i)) ';']);
        timestamps = [];
        timeDiffToRefArray = [];
        dutSamplesForLatency = length(signalsDUTchannel);
        if dutSamplesForLatency>maxNumberChannelForLatency
            dutSamplesForLatency = maxNumberChannelForLatency;
        end
        
        
        for k = 1:dutSamplesForLatency
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
            if minTimeDiffToRef<200
                timeDiffToRefArray = [timeDiffToRefArray;minTimeDiffToRef];
            end
        end
        
        %determine latency from min time diff histogram between ch and ref
        latency = median (timeDiffToRefArray);
        % latency = -3; %override to set latency manually
        
        figure
        hold on
        hist(timeDiffToRefArray,100);
        if isnan(latency)
            
        else
            xline(latency,'color','red','linewidth',2)
        end
        xlabel('Min time between DUT and REF (ns)');
        ylabel('Events');
        title(['SAMPIC - Minimum time between DUT and REF - CH: ' num2str(channelsEnabled(i)) ' - Run ' run.id]);
        grid on
        saveas(gcf,[store_folder '\Run' run.id '_latency_CH' num2str(channelsEnabled(i)) '.png'])
        
        pause(1);
        
        
        eval(['latency_ch' num2str(channelsEnabled(i)) ' = latency;']);
        
    end
    
end

channelsEnabled



%%loop through all segements
for segmentCounter = 1:numberSegments
clear signals
        signalSegmentPath = [signalsFolder '\segment' num2str(segmentCounter) '.mat'];
        load(signalSegmentPath);

           for channelPos = 1:length(signals)
            signalChannels (channelPos) = signals(channelPos).ch;
        end     
 
    

%% loop through DUT channels  - associate to reference
for i = 1:length(channelsEnabled)
    if  channelsEnabled(i) == referenceDetectorChannel
        
    else
        %is dut channel
                signalsCut = signalChannels==channelsEnabled(i);
        signalsDUTchannel = signals(signalsCut);

        %eval(['signalsDUTchannel = signals_ch' num2str(channelsEnabled(i)) ';']);
        timestamps = [];
        timeDiffToRefArray = [];
        
        %loop through all events and accept ones in window around latency
        
       % for k = 1:length(signalsDUTchannel)
%             dutSignal = signalsDUTchannel(k);
%             timeDifferencesDUTRef = [];
%             detectorHitTime = dutSignal.cell0Time;
%             detectorHitTimeUnix = dutSignal.unixTimestamp;
%             timestampEntry = [detectorHitTime detectorHitTimeUnix];
%             timestamps = [timestamps; timestampEntry];
%             timeDiffToRef = detectorHitTime-refTime;
            
            eval(['latency = latency_ch' num2str(channelsEnabled(i)) ';']);
            
            
            minTime = latency-matchingDUTREFWindowWidth;
            maxTime = latency+matchingDUTREFWindowWidth;
            
            matchedEventsTemp = [];
            
            %loop through events of this channel to find matches
            for k = 1:length(signalsDUTchannel)
                
                if mod(k,1000)==0
                    str_disp=sprintf('Processing DUT hit No. %d', k);
                    disp(str_disp);
                end
                
                dutSignal = signalsDUTchannel(k);
                timeDifferencesDUTRef = [];
                detectorHitTime = dutSignal.cell0Time;
                timeDiff = signalsRefTime-detectorHitTime;
                timeDiffTracker = scaledEventCounterTimes-detectorHitTime;
                
                timeMatchToRef = timeDiff>minTime & timeDiff<maxTime;
                
                matchIdx = find(timeMatchToRef==1);
                
                [trackerMinValue,trackerMinIdx] = min(abs(timeDiffTracker));
                trackerEntryMatched = eventCounterEntries(trackerMinIdx,:);
                
                if length(matchIdx==1) %unambiguous match
                    %found a single match in reference for this dut event -> unambiguois match
                    event.dut = signalsDUTchannel(k);
                    event.ref = signalsRef(matchIdx(1));
                    event.timeDiff = timeDiff;
                    
                    if trackerMinValue<50000
                        event.eventID = trackerEntryMatched(3);
                    else
                        event.eventID = 0;
                    end
                    event.trackerTimeDiff = trackerMinValue;
                    if length(matchedEventsTemp)==0
                        matchedEventsTemp = [event];
                    else
                        matchedEventsTemp(k) = event;
                    end
                   % matchedEventsTemp = [matchedEventsTemp;event];
                elseif  length(matchIdx>=1)
                    %found multiple matches
                    matchIdx
                end
                
            end
            
            %eval(['matchedEvents_ch' num2str(channelsEnabled(i)) ' = matchedEventsTemp;']);
            save([matchedEventsFolder '\ch' num2str(channelsEnabled(i)) '_segment' num2str(segmentCounter) '.mat'], 'matchedEventsTemp');

        end
        
    %end
end
channelsEnabled

end

