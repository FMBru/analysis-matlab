close all

%%settings for processing

addpath 'C:\Users\GDD\Documents\MATLAB\Picosec\CommonFunctions';
addpath 'C:\Users\GDD\Documents\MATLAB\Picosec\CommonFunctions\SCP';

% clear everything if not batch processing
%set run parameters
if exist('batchProcess','var') == 1
    run.id = runIDString;
    run.oscilloscope = oscilloscopeString;
    run.lecroy_name = '--Trace--'; %['Run' run.id];
    opts_MM.chID = 4;
    opts_MCP.chID = 0;
    opts_MM.ch_name = [DUTnameString run.lecroy_name]
    opts_MCP.ch_name = [MCPnameString run.lecroy_name]
    shouldSaveMAT = false;
    shouldUseEOSFolder = useEOSFolder;
    tracker.dutIndex = str2num(trackerDUTIndex); %1: miniCactus, 2: MM3 (VacChamberPos), 3: MM1(Multipad), 4: MM2 (on support plate), 5: MM4 (electron setup)
    
else
    %not batch processing
    %clear all
    run.lecroy_name = '--Trace--'; %['Run' run.id];
    
    %run.id = '362';
    run.oscilloscope = 'SAMPIC';
    opts_MM.chID = 22;
    opts_MCP.chID = 0;
    opts_MM.ch_name = ['C' num2str(opts_MM.chID)]; % LeCroy file name format MCP1 (used for timing)
    opts_MCP.ch_name = ['C' num2str(opts_MCP.chID)]; % MPC 2 (used as trigger)
    shouldSaveMAT = false;
    shouldUseEOSFolder = true;
    tracker.dutIndex = 3; %1: miniCactus, 2: MM3 (VacChamberPos), 3: MM1(Multipad), 4: MM2 (on support plate), 5: MM4 (electron setup)
    
end



% run file coordinates and number of files
minuit = 1;   % set minuit for optimizer if 1 or matlab 0
run.year = '2022/05 ';
run.name = ['BEAM ' run.year ' RUN ' run.id];
%run.path=['M:\CommonProject\PICOSEC\Testbeams\July2021\Oscilloscopes\POOL2\Run' run.id '\']; % use mounted Cernbox folder
% run.path=['D:\Matlab_timing\Run' run.id '\']; % use local folder with scope data
run.path=['C:\Users\GDD\Documents\Picosec\May22\' run.oscilloscope '\Run' run.id '\Run' run.id '\'];
run.path=['C:\Users\GDD\Documents\Picosec\May22\' run.oscilloscope '\Run' run.id '\'];
run.pathEOS=['\\eosproject-smb\eos\project\p\picosec\testbeam\2022_May_h4\' run.oscilloscope '\Run' run.id '\'];

if shouldUseEOSFolder
    tracker.path = ['\\eosproject-smb\eos\project\p\picosec\testbeam\2022_May_h4\tracker\reconstructed\asciiRun' run.id '.dat'];
else
    tracker.path = ['C:\Users\GDD\Documents\Picosec\May22\Tracker\reconstructed\asciiRun' run.id '.dat'];
end

%read tracker data from ASCII file
tracker.en = 1; %match eventIDs to tracking data and add XY to output

dutColArray = [1 4 5; 2 7 8; 3 10 11; 4 13 14; 5 16 17]; %[dID colXID colYID; ] -> May MM run

if tracker.en
    trackerFile = fopen(tracker.path,'rt');
    D = textscan(trackerFile, '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f', 'Delimiter','\t', 'HeaderLines',2, 'CollectOutput',1);
    tracker.data = cell2mat(D);
end

% options for processing micromegas channel (lot magic numbers from
% processing functionhave to be added here)
opts_MM.t_dead = 1;      % blanking time from start of the recording samples
opts_MM.t_prerms = 5;  % blanking time before global maximum samples
opts_MM.Ts = 1/20e9;     % sampling speed
opts_MM.Rin=50;          % input impedance
opts_MM.invert = 1;      % is inverted (1-inverted, 0-non-inverted)
opts_MM.type=0;          % detector type (0-DUT, 1-REF) for future use
%only for override - should set in top
%opts_MM.ch_name = ['C4' run.lecroy_name]; % LeCroy file name format MCP1 (used for timing)
opts_MM.en_plot = 0;     % enable debugging plots

% options for processing MCP channel (lot magic numbers from processing
% functionhave to be added here)
opts_MCP.t_dead = 1;
opts_MCP.t_prerms = 5;
opts_MCP.Ts = 1/20e9;
opts_MCP.Rin=50;
opts_MCP.invert = 1;
opts_MCP.type=1;
%only for override - should set in top
%opts_MCP.ch_name = ['C1' run.lecroy_name]; % MPC 2 (used as trigger)
opts_MCP.en_plot = 0;



%%settings for matching
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
channelsEnabled = [0;22];


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

eventPos = 1;
eventPosAnalysis = 1;

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
            
            matchedEventsDUT = [];
            
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
                    if length(matchedEventsDUT)==0
                        matchedEventsDUT = [event];
                    else
                        matchedEventsDUT(eventPos) = event;
                    end
                    % matchedEventsTemp = [matchedEventsTemp;event];
                    eventPos = eventPos+1;
                elseif  length(matchIdx>=1)
                    %found multiple matches
                    matchIdx
                end
                
            end
            matchedEventsNumber = length(matchedEventsDUT)
            %eval(['matchedEvents_ch' num2str(channelsEnabled(i)) ' = matchedEventsTemp;']);
            %save([matchedEventsFolder '\ch' num2str(channelsEnabled(i)) '_segment' num2str(segmentCounter) '.mat'], 'matchedEventsTemp');
            
            
            %%do processing for matched events in DUT
            %% DO PROCCEISNG
            tic;         % start time measurement
            k=1;         % valid data counter
            
            
            % go trough all of the events in matched events vector
            for m=1:length(matchedEventsDUT)
                
                if mod(m,1000)==0
                    str_disp=sprintf('Analyzing No. %d', m);
                    disp(str_disp);
                end
                
                %structure of device hit
                %cell0Time
                %unixTimestamp
                %ch
                %hit
                %waveform
                if length(matchedEventsDUT(m).dut)==0 || length(matchedEventsDUT(m).ref)==0
                    continue
                end
                dut = matchedEventsDUT(m).dut;
                dutWaveform = dut.waveform;
                
                ref = matchedEventsDUT(m).ref;
                refWaveform = ref.waveform;
                
                event_id = matchedEventsDUT(m).eventID;
                
                % generate virtual time vector (important to take care for trigger
                % offset in LeCroy scope)
                t_vec_mm=dutWaveform(:,1);
                t_vec_mcp=refWaveform(:,1);
                
                % subtract the earliest time
                etime=min([t_vec_mm; t_vec_mcp]);
                
                % process MM Picosec first to see if signal is valid
                if(minuit==1)
                    %MM_temp = process_signal_minuit(t_vec_mm-etime,dutWaveform(:,2),opts_MM)
                    MM_temp = process_signal_sampic(t_vec_mm-etime,dutWaveform(:,2),opts_MM,1);
                else
                    MM_temp = process_signal(t_vec_mm-etime,dutWaveform(:,2),opts_MM);
                end
                
                % if signal is valid process MCP
                if(MM_temp.fail==0)
                    if(minuit==1)
                        %MCP_temp = process_signal_minuit(t_vec_mcp-etime,refWaveform(:,2),opts_MCP);
                        MCP_temp = process_signal_sampic(t_vec_mcp-etime,refWaveform(:,2),opts_MCP,1);
                    else
                        MCP_temp = process_signal(t_vec_mcp-etime,refWaveform(:,2),opts_MCP);
                    end
                    
                    % if MCP is valid store data to structure array
                    if(MCP_temp.fail==0)
                        % process tracker ID channel
                        
                        MM_temp.event_id = event_id;
                        MCP_temp.event_id = event_id;
                        
                        
                        xPos = 0;
                        yPos = 0;
                        if tracker.en
                            %find corresponding tracker entry and save XY info
                            trackIdx = find(tracker.data(:,1) == event_id);
                            if trackIdx>0
                                %is valid trackIdx -> get XY
                                xPos = tracker.data(trackIdx,dutColArray(tracker.dutIndex,2));
                                yPos = tracker.data(trackIdx,dutColArray(tracker.dutIndex,3));
                                if size(xPos,1)>1
                                    xPos = xPos(1);
                                end
                                if size(yPos,1)>1
                                    yPos = yPos(1);
                                end
                            end
                            
                        end
                        MM_temp.x = xPos;
                        MM_temp.y = yPos;
                        trackerX(eventPosAnalysis) = MM_temp.x;
                        trackerY(eventPosAnalysis) = MM_temp.y;
                        eventIDArray(eventPosAnalysis) = MM_temp.event_id;
                        
                        
                        % store valid data into array of structures
                        MM_data(eventPosAnalysis)= MM_temp;
                        MCP_data(eventPosAnalysis)= MCP_temp;
                        time_diff(eventPosAnalysis) = MM_data(k).cfd.time-MCP_data(eventPosAnalysis).cfd.time;
                        time_diff_sigmoid(eventPosAnalysis) = MM_data(k).sigmoid.timepoint-MCP_data(eventPosAnalysis).sigmoid.timepoint;
                        MCP_maxy(eventPosAnalysis) = MCP_data(eventPosAnalysis).sig.max.y;
                        MM_maxy(eventPosAnalysis) = MM_data(eventPosAnalysis).sig.max.y;
                        trackerX(eventPosAnalysis) = MM_temp.x;
                        trackerY(eventPosAnalysis) = MM_temp.y;
                        
                        eventPosAnalysis=eventPosAnalysis+1;
                    end
                end
            end
            
            toc
            
            
        end
        
        
        %end
    end
    channelsEnabled;
    
end

