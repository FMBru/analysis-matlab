%clear all
close all
tic

shouldAnalyseWhenFinished = true;

eventPosProcessing = 1; %counter of processed events, comment when combining runs

%%when combining data from multiple runs, do not reset eventPosAnalysis and
%%do not clear all

addpath 'C:\Users\GDD\Documents\MATLAB\Picosec\CommonFunctions';
addpath 'C:\Users\GDD\Documents\MATLAB\Picosec\CommonFunctions\SCP';
addpath 'C:\Users\GDD\Documents\MATLAB\Picosec\SAMPIC\CrateSystem\Functions'

hitsPerSegment = 30000; %how many hits to save per segment

%reduced sampling for quick look at data - better to use ref SignalsStep
matchedEventsStep = 1; %step size when looping through events to reconstruct, set 1 for all, set higher for fast analysis of reduced data set
refSignalsStep = 1; %step size when looping through ref signals to reconstruct, set 1 for all, set higher for fast analysis of reduced data set


run.id = '081';
referenceDetectorChannel = 0; %channel number of detector used as timing reference - for latency relative to this detector
triggerSignalsChannel = 57; %channel number to which trigger signal is connected which is sent to SRS for sync with tracker
%channelToAnalyse = [1;2;3;4;5;6;7;8;9;10;11;12;13;14;15;16;17;18;19;20;21;22;23;24;25;26;27;28;29;30;31;32;33;34;35;36;37;38;39;40;41;42;43;44;45;46;47;48;49;50];
%channelToAnalyse = [1;2;3;4;5];
channelToAnalyse = 1:63;

run.path=['\\eosproject-smb\eos\project\p\picosec\testbeam\2022_October_h4\SAMPIC\Run' run.id '\'];
%run.path = '\\eosproject-smb\eos\project\p\picosec\testbeam\2021_October_h4\SAMPIC\Run297a\Run_SAMPIC_297a_Data_10_30_2021_9h_46min_Ascii\'
SAMPICFolderInfo = dir(run.path);

foundTriggerFile = false;
foundAsciiDataFile = false;
foundBinaryDataFile = false;
foundSegmentedBinaryDataFile = false;

dataFilePath = '';
triggerFilePath = '';

store_folder = ['\\eosproject-smb\eos\project\p\picosec\testbeam\2022_October_h4\Results\Run' run.id '-SAMPIC'];
mkdir(store_folder);

processingFolder = ['\\eosproject-smb\eos\project\p\picosec\testbeam\2022_October_h4\SAMPIC\Processing\Run' run.id ];
signalsFolder = [processingFolder '\Signals'];
mkdir(signalsFolder);

numberHitsForLatency = 100000;
processingSegmentSize = 50000;

numberHeaderLines = 7;

matchingDUTREFWindowWidth=3;

channelsEnabled = []; %ids of enabled channels


%load mapping
setMapping



%% config for processing


run.oscilloscope = 'SAMPIC';
opts_MM.ch_name = ['DUT CH']; 
opts_MCP.ch_name = ['C' num2str(referenceDetectorChannel)]; % MPC 2 (used as trigger)
shouldSaveMAT = false;
shouldUseEOSFolder = true;
tracker.dutIndex = 3; %1: miniCactus, 2: MM3 (VacChamberPos), 3: MM1(Multipad), 4: MM2 (on support plate), 5: MM4 (electron setup)

% run file coordinates and number of files
minuit = 1;   % set minuit for optimizer if 1 or matlab 0
run.year = '2022/10 ';
run.name = ['BEAM ' run.year ' RUN ' run.id];
run.pathEOS=['\\eosproject-smb\eos\project\p\picosec\testbeam\2022_October_h4\' run.oscilloscope '\Run' run.id '\'];

if shouldUseEOSFolder
    tracker.path = ['\\eosproject-smb\eos\project\p\picosec\testbeam\2022_October_h4\tracker\reconstructed\asciiRun' run.id '.dat'];
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

% options for processing micromegas channel
opts_MM.t_dead = 1;      % blanking time from start of the recording samples
opts_MM.t_prerms = 5;  % blanking time before global maximum samples
opts_MM.Ts = 1/20e9;     % sampling speed
opts_MM.Rin=50;          % input impedance
opts_MM.invert = 1;      % is inverted (1-inverted, 0-non-inverted)
opts_MM.type=0;          % detector type (0-DUT, 1-REF) for future use
%only for override - should set in top
%opts_MM.ch_name = ['C4' run.lecroy_name]; % LeCroy file name format MCP1 (used for timing)
opts_MM.en_plot = 0;     % enable debugging plots

% options for processing MCP channel
opts_MCP.t_dead = 1;
opts_MCP.t_prerms = 5;
opts_MCP.Ts = 1/20e9;
opts_MCP.Rin=50;
opts_MCP.invert = 1;
opts_MCP.type=1;
%only for override - should set in top
%opts_MCP.ch_name = ['C1' run.lecroy_name]; % MPC 2 (used as trigger)
opts_MCP.en_plot = 0;

%% start processing data
%check for segmented files
for i=1:length(SAMPICFolderInfo)
    fileInfo = SAMPICFolderInfo(i);
    if contains(fileInfo.name,'.bin_0001')                %check if ascii data file
        dataFilePath=[run.path fileInfo.name];
        if exist(dataFilePath,'file')==2
            foundSegmentedBinaryDataFile= true;
        end
    end
end



segmentedBinaryFilesArray = [];


for i=1:length(SAMPICFolderInfo)
    fileInfo = SAMPICFolderInfo(i);
    if contains(fileInfo.name,'SAMPIC_Trigger_Data')    %check if trigger file
        triggerFilePath=[run.path fileInfo.name];
        if exist(triggerFilePath,'file')==2
            foundTriggerFile= true;
            
        end
    elseif contains(fileInfo.name,'.dat')                %check if ascii data file
        dataFilePath=[run.path fileInfo.name];
        if exist(dataFilePath,'file')==2
            foundAsciiDataFile= true;
        end
    elseif contains(fileInfo.name,'.bin_')  %check if binary data file
        %is segmented file
        foundSegmentedBinaryDataFile = true;
        segmentedBinaryFilesArray = [segmentedBinaryFilesArray;fileInfo.name];
    elseif contains(fileInfo.name,'.bin')  %check if binary data file
        %is binary data file, check if single file or multiple
        dataFilePath=[run.path fileInfo.name];
        if exist(dataFilePath,'file')==2
            foundBinaryDataFile= true;
        end
    end
end

% segmentedBinaryFilesArray
%
% foundTriggerFile
% foundAsciiDataFile
% foundBinaryDataFile
% foundSegmentedBinaryDataFile
% dataFilePath

%read in trigger file
eventCounterEntries = readSAMPICTriggerBinary(triggerFilePath);
plot(eventCounterEntries(:,2),eventCounterEntries(:,1),'.');
xlabel('Time (ns)');
ylabel('SRS Event ID');
title(['SAMPIC - decoded SRS event ID over time - Run ' run.id]);
grid on
pause(1);
saveas(gcf,[store_folder '\Run' run.id '_SRSeventCounter.png'])

scaledEventCounterTimes = eventCounterEntries(:,2);

signals = [];

if foundAsciiDataFile
    %% read in ascii file
    fid = fopen(dataFilePath);
    
    %init array for file reading
    tline = fgetl(fid);
    lineNumber = 1;
    lineString = '';
    samplingFrequency = 1;
    timestep = 1;
    firstCell0Time = 0;
    cell0TimeArray = [];
    unixTimeArray = [];
    
    str_disp=sprintf('Reading in ASCII file');
    disp(str_disp);
    
    i=1;
    while ischar(tline) && lineNumber<100000000000
        
        if mod(lineNumber,1000)==0
            str_disp=sprintf('Reading line no. %d', lineNumber);
            disp(str_disp);
            
        end
        % disp(tline);
        lineString = tline;
        startString = lineString(1:3);
        if startString=='==='
            %is comment
            %disp('comment');
            if contains(lineString,'SAMPLING FREQUENCY')
                samplingFrequencyLineSplit = split(lineString);
                samplingFrequency = str2num(samplingFrequencyLineSplit{4})*1000000 ;
                timestep = 1/samplingFrequency;
            end
        else
            %disp('dataLine');
            %disp(startString);
            %is data line
            if startString=='Hit'
                hitLine = lineString;
                chLine = fgetl(fid);
                samplesLine = fgetl(fid);
                
                hitLineSplit = split(hitLine);
                hitNumber = str2num(hitLineSplit{2});
                hitUnixTime = str2num(hitLineSplit{5});
                
                chLineSplit = split(chLine);
                chNumber = str2num(chLineSplit{2}) ;
                cell0Time = str2num(chLineSplit{4}) ;
                
                if firstCell0Time == 0
                    firstCell0Time = cell0Time;
                end
                
                samplesLineSplit = split(samplesLine);
                numberDataSampels = length(samplesLineSplit)-2;
                %samplesVector = str2num(cell2mat(samplesLineSplit(2:(length(samplesLineSplit)-1))));
                
                sampleVector = [];
                for pos=1:(length(samplesLineSplit)-2)
                    pickIndex = pos+1;
                    sampleVector(pos) = str2num(samplesLineSplit{pickIndex});
                end
                
                %reducedCell0Time = cell0Time-firstCell0Time;
                reducedCell0Time = cell0Time;
                
                numberSamples = length(sampleVector);
                waveform = zeros(numberSamples,2);
                waveform(1:numberSamples,1) = (1:numberSamples)-1;
                waveform(1:numberSamples,1) = waveform(1:numberSamples,1)*timestep+cell0Time*0.000000001;
                waveform(1:numberSamples,2) = sampleVector;
                
                cell0TimeArray = [cell0TimeArray;reducedCell0Time];
                unixTimeArray = [unixTimeArray;hitUnixTime];
                
                hit.cell0Time = reducedCell0Time;
                hit.unixTimestamp = hitUnixTime; %or hitUnixTime
                hit.ch = chNumber;
                hit.hit = hitNumber;
                hit.waveform = waveform;
                
                
                if chNumber == referenceDetectorChannel
                    %is reference detector, store separately
                    
                    
                    if length(signalsRef)==0
                        signalsRef = [hit];
                        
                    else
                        signalsRef (length(signalsRef)+1) = hit;
                    end
                    
                    
                else
                    
                    if length(signals)==0
                        signals = [hit];
                        
                    else
                        signals (length(signals)+1) = hit;
                        % signals = [signals;hit];
                    end
                    i = i+1;
                    %                    length(signals);
                    %                     if length(signals)>hitsPerSegment
                    %                         %save signals as part and reset array
                    %                         segmentCounter
                    %                         save([signalsFolder '\segment' num2str(segmentCounter) '.mat'], 'signals');
                    %                         clear signals;
                    %                         signals = [];
                    %                         segmentCounter=segmentCounter+1;
                    %                     end
                end
                
                
                
                
                % signals = [signals;hit];
            end
        end
        
        
        tline = fgetl(fid);
        lineNumber = lineNumber+1;
    end
    fclose(fid);
    
elseif foundBinaryDataFile
    str_disp=sprintf('Reading in Binary data file');
    disp(str_disp);
    
    signals = [];
    chIDsArray = [];
    chTimesArray = [];
    
    
    if foundSegmentedBinaryDataFile
        dataFilePath=[run.path segmentedBinaryFilesArray(1,:)];
    end
    
    fid = fopen(dataFilePath);
    for i=1:numberHeaderLines
        tline = fgets(fid);
        startString = tline(1:3);
        
        if startString=='==='
            %is comment
            %disp('comment');
            if contains(tline,'SAMPLING FREQUENCY')
                samplingFrequencyLineSplit = split(tline);
                samplingFrequency = str2num(samplingFrequencyLineSplit{4})*1000000
                timestep = 1/samplingFrequency;
            end
        else
        end
    end
    pause(1)
    
    %going through some first lines to use for latency calculation
    i = 1;
    while ~feof(fid) && i<numberHitsForLatency
        if mod(i,10000)==0
            toc
            tic
            str_disp=sprintf('Reading in hit No. %d', i);
            disp(str_disp);
        end
        hitNumber = fread(fid,1,'int32',0);
        epochTime = fread(fid,1,'double',0);
        channel = fread(fid,1,'int32',0);
        TriggerCellTimeInstant = fread(fid,1,'double',0);
        CFDTimeInstant = fread(fid,1,'double',0);
        Baseline = fread(fid,1,'float64',0);
        PeakValue = fread(fid,1,'float64',0);
        Amplitude = fread(fid,1,'float32',0);
        dataSize = fread(fid,1,'int32',0);
        dataLength = dataSize;
        
        sampleVector = zeros(dataLength,1);
        for pos=1:dataLength
            sampleVector(pos) = fread(fid,1,'float32',0)  ;
        end
        
        numberSamples = length(sampleVector);
        waveform = zeros(numberSamples,2);
        waveform(1:numberSamples,1) = (1:numberSamples)-1;
        waveform(1:numberSamples,1) = waveform(1:numberSamples,1)*timestep+TriggerCellTimeInstant*0.000000001;
        waveform(1:numberSamples,2) = sampleVector;
        
        hit.cell0Time = TriggerCellTimeInstant;
        %hit.unixTimestamp = epochTime; %or hitUnixTime
        hit.ch = channel;
        %hit.hit = hitNumber;
        hit.waveform = waveform;
        
        if length(signals)==0
            signals = [hit];
            chIDsArray = [channel];
            chTimesArray = [TriggerCellTimeInstant];
            
        else
            if dataSize>10
                signals (length(signals)+1) = hit;
                chIDsArray (length(signals)) = channel;
                chTimesArray (length(signals)) = TriggerCellTimeInstant;
                
            end
        end
        i = i+1;
    end
    fclose(fid);
    
    
    
    
    
    
    %first numberFilesForLatency in signals array
    
    %% determine enabled channels from signals array
    
    for i = 1:length(signals)
        if length(channelsEnabled)==0 || any(ismember(channelsEnabled,signals(i).ch)) == 0
            channelsEnabled = [channelsEnabled;signals(i).ch];
        end
    end
    channelsEnabled = sort(channelsEnabled,'ascend');
    numberChannelsEnabled = length(channelsEnabled);
    chCountsArray = zeros(numberChannelsEnabled,1);
    chLatencyArray = zeros(numberChannelsEnabled,1);
    
    for chPos = 1:length(channelsEnabled)
        chID = channelsEnabled(chPos);
        chIDsMask = chIDsArray==chID;
        chCountsArray(chPos,1) = length(chIDsArray(chIDsMask));
    end
    chCountsArrayNumbers = chCountsArray;
    chCountsArray = chCountsArray/sum(chCountsArray);
    bar(channelsEnabled,chCountsArray);
    xlabel('Channel ID');
    ylabel('Hits (fraction of total) ');
    title(['Hits per SAMPIC channel - Run ' run.id]);
    grid on
    xlim([0 63]);
    pause(1);
    saveas(gcf,[store_folder '\Run' run.id '_CHHits.png'])
    %plot hitmap accross multipad
    VisualiseMultipad(channelsEnabled,chCountsArray,'Hitmap Multipad Picosec Relative',[store_folder '\Run' run.id '_HitmapRelative.png'],'Hits: %0.1f %');
    VisualiseMultipad(channelsEnabled,chCountsArrayNumbers,'Hitmap Multipad Picosec',[store_folder '\Run' run.id '_Hitmap.png'],'Hits: %0.1f %');
    
    %%determine latency for channels relative to reference detector
    %refIndexInEnabledChannels = find(channelsEnabled==referenceDetectorChannel);
    refChannelMask = chIDsArray==referenceDetectorChannel;
    refTime = chTimesArray(refChannelMask);
    
    %% loop through DUT channels  - get latency
    for chPos = 1:length(channelsEnabled)
        if  channelsEnabled(chPos) == referenceDetectorChannel
            %is ref, get latency for SRS to ref
            
            %eval(['signalsDUTchannel = signals_ch' num2str(channelsEnabled(i)) ';']);
            timestamps = [];
            timeDiffToRefArray = [];
            
            
            
            for k = 1:length(refTime)
                %loop through all events in DUT to find time differences to
                refSignal = refTime(k);
                timeDifferencesDUTRef = [];
                detectorHitTime = refSignal;
                timeDiffToRef = detectorHitTime-scaledEventCounterTimes;
                minTimeDiffToRef = min(abs(timeDiffToRef));
                timeDiffToRefArray = [timeDiffToRefArray;minTimeDiffToRef];
            end
            
            timeDiffToRefArrayMask = timeDiffToRefArray<1e6;
            timeDiffToRefArrayCut = timeDiffToRefArray(timeDiffToRefArrayMask);
            
            %determine latency from min time diff histogram between ch and ref
            latency = median (timeDiffToRefArrayCut);
            eventCounterLatency = latency;
            % latency = -3; %override to set latency manually
            %latencyEntry = [scalingFactor latency numberLow numberHigh ratio];
            %latencyEntry = [latency];
            %medianLatencyErrorArray = [medianLatencyErrorArray;latencyEntry];
            figure
            hold on
            hist(timeDiffToRefArrayCut,1000);
            if isnan(latency)
                
            else
                %    xline(latency,'color','red','linewidth',2)
            end
            xlabel('Min time between SRS counter and REF (ns)');
            ylabel('Events');
            title(['SAMPIC - Minimum time between EventCounter to RefChannel:  - Run ' run.id]);
            grid on
            xlim([0 30000]);
            saveas(gcf,[store_folder '\Run' run.id '_latency_counter_scaling.png'])
            %close all
            pause(1);
            
            [refTimeCounts,refTimeCenters] = hist(refTime,10000);
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
            title(['SAMPIC - Minimum time between EventCounter and Ref - Run ' run.id]);
            legend('Ref hits time','Event Counter time')
            grid on
            
            saveas(gcf,[store_folder '\Run' run.id '_timeComparison.png'])
            %close all
            pause(1);
        else
            %is dut channel
            
            %load first segment of this channel for latency calculation
            signalsCut = chIDsArray==channelsEnabled(chPos);
            signalsDUTchannel = signals(signalsCut);
            
            % eval(['signalsDUTchannel = signals_ch' num2str(channelsEnabled(i)) ';']);
            timestamps = [];
            timeDiffToRefArray = [];
            
            for k = 1:length(signalsDUTchannel)
                %loop through all events in DUT to find time differences to
                %reference
                dutSignal = signalsDUTchannel(k);
                timeDifferencesDUTRef = [];
                detectorHitTime = dutSignal.cell0Time;
                %timestampEntry = [detectorHitTime detectorHitTimeUnix];
                % timestamps = [timestamps; timestampEntry];
                timeDiffToRef = detectorHitTime-refTime;
                minTimeDiffToRef = min(abs(timeDiffToRef));
                if minTimeDiffToRef<200 && minTimeDiffToRef~0;
                    timeDiffToRefArray = [timeDiffToRefArray;minTimeDiffToRef];
                end
            end
            
            %determine latency from min time diff histogram between ch and ref
            latency = median (timeDiffToRefArray);
            % latency = -3; %override to set latency manually
            
            figure
            hold on
            hist(timeDiffToRefArray,1000);
            if isnan(latency)
                
            else
                xline(latency,'color','red','linewidth',1);
            end
            xlabel('Min time between DUT and REF (ns)');
            ylabel('Events');
            title(['SAMPIC - Minimum time between DUT and REF - CH: ' num2str(channelsEnabled(chPos)) ' - Run ' run.id]);
            grid on
            xlim([0 200]);
            saveas(gcf,[store_folder '\Run' run.id '_latency_CH' num2str(channelsEnabled(chPos)) '.png']);
            
            % pause(1);
            
            
            eval(['latency_ch' num2str(channelsEnabled(chPos)) ' = latency;']);
            chLatencyArray(chPos,1) = latency;
        end
        
    end
    figure
    plot(channelsEnabled,chLatencyArray,'.');
    xlabel('Channel ID');
    ylabel('Latency to ref (ns) ');
    title(['Latency to Ref SAMPIC channel - Run ' run.id]);
    grid on
    ylim([0 200]);
    xlim([0 63]);
    pause(1);
    saveas(gcf,[store_folder '\Run' run.id '_CHLatencyToRef.png'])
    
    %plot latency accross multipad
    VisualiseMultipad(channelsEnabled,abs(chLatencyArray),'Latency Multipad Picosec',[store_folder '\Run' run.id '_Latency.png'],'Latency: %0.1f ns');
    
    
    %% latencies determined - process binary file in segments
    
    %for timestamp overflows
    timestamp_prev = -1;
    timestamp_ov = 0;

    str_disp=sprintf('Reading in Binary data file for processing');
    disp(str_disp);
    
    if foundSegmentedBinaryDataFile == false
        segmentedBinaryFilesArray=[dataFilePath];
        currentFilePath = dataFilePath;
    end


    %loop through files
    for filePos=1:size(segmentedBinaryFilesArray,1)
        
                  str_disp=sprintf('Starting Binary File No. %d', filePos);
                disp(str_disp);

        
        currentFilePath = [run.path segmentedBinaryFilesArray(filePos,:)];
        if foundSegmentedBinaryDataFile == false
            currentFilePath = dataFilePath;
        end
        
        fid = fopen(currentFilePath);
        
        for i=1:numberHeaderLines
            tline = fgets(fid);
        end
        
        %going through all files for data
        i = 1;
        
        signals = [];
        chIDsArray = [];
        chTimesArray = [];
        chIndexArray = [];
        
        
        while ~feof(fid) %&& eventPosAnalysis<500000 %process full file
            if mod(i,10000)==0
                toc
                tic
                %str_disp=sprintf('Reading in hit No. %d', i);
                %disp(str_disp);
            end
            hitNumber = fread(fid,1,'int32',0);
            epochTime = fread(fid,1,'double',0);
            channel = fread(fid,1,'int32',0);
            TriggerCellTimeInstantRaw = fread(fid,1,'double',0);
            CFDTimeInstant = fread(fid,1,'double',0);
            Baseline = fread(fid,1,'float64',0);
            PeakValue = fread(fid,1,'float64',0);
            Amplitude = fread(fid,1,'float32',0);
            dataSize = fread(fid,1,'int32',0);
            dataLength = dataSize;
            
            
            % count overflows
            if(TriggerCellTimeInstantRaw < (timestamp_prev-10000))
                timestamp_ov = timestamp_ov + 1
            end
            timestamp_prev = TriggerCellTimeInstantRaw;

           % pause(1);
            %previous - is timestamp_prev-9511627776 needed?
%             if(TriggerCellTimeInstantRaw < (timestamp_prev-9511627776))
%                 timestamp_ov = timestamp_ov + 1
%             end
%             timestamp_prev = TriggerCellTimeInstantRaw;

            TriggerCellTimeInstant = timestamp_ov * 1099511627776 + TriggerCellTimeInstantRaw;
            
            
            
            sampleVector = zeros(dataLength,1);
            for pos=1:dataLength
                sampleVector(pos) = fread(fid,1,'float32',0)  ;
            end
            
            numberSamples = length(sampleVector);
            waveform = zeros(numberSamples,2);
            waveform(1:numberSamples,1) = (1:numberSamples)-1;
            waveform(1:numberSamples,1) = waveform(1:numberSamples,1)*timestep+TriggerCellTimeInstant*0.000000001;
            waveform(1:numberSamples,2) = sampleVector;
            
            hit.cell0Time = TriggerCellTimeInstant;
            %hit.unixTimestamp = epochTime; %or hitUnixTime
            hit.ch = channel;
            %hit.hit = hitNumber;
            hit.waveform = waveform;
            
            if length(signals)==0
                signals = [hit];
                chIDsArray = [channel];
                chTimesArray = [TriggerCellTimeInstant];
                chIndexArray = [1];
            else
                if dataSize>10
                    signals (length(signals)+1) = hit;
                    chIDsArray (length(signals)) = channel;
                    chTimesArray (length(signals)) = TriggerCellTimeInstant;
                    chIndexArray (length(signals)) = length(signals);
                end
            end
            
            eventCounterTempLatencyArray = [];
            
            if length(signals)>50000
                %% match and process signals
                str_disp=sprintf('Matching events', i);
                disp(str_disp);
                matchedEvents = [];
                
                %go through all ref signals
                refChannelMask = chIDsArray==referenceDetectorChannel;
                chTimesArrayRef = chTimesArray(refChannelMask);
                chIndexArrayRef = chIndexArray(refChannelMask);
                signalsRef = signals(refChannelMask);
                
                
                for refPos = 1:refSignalsStep:length(signalsRef)
                    refSignal = signalsRef(refPos);
                    
                    %go through all channels and check if there is a match
                    matchedChannels = [];
                    for chDUTPos = 1:length(channelsEnabled)
                        if  channelsEnabled(chDUTPos) == referenceDetectorChannel
                        else
                            %is dut channel
                            latency = chLatencyArray(chDUTPos);
                            dutChannelMask = chIDsArray==channelsEnabled(chDUTPos);
                            chTimesArrayDUT = chTimesArray(dutChannelMask);
                            chIndexArrayDUT = chIndexArray(dutChannelMask);
                            
                            minTime = latency-matchingDUTREFWindowWidth;
                            maxTime = latency+matchingDUTREFWindowWidth;
                            timeDiffArray = abs(chTimesArrayDUT-refSignal.cell0Time);
                            matchMask = timeDiffArray>minTime & timeDiffArray<maxTime;
                            
                            chIndexMatched = chIndexArrayDUT(matchMask);
                            if length(chIndexMatched)==1 %unambiguous hit
                                %accept as match
                                dutTemp = signals(chIndexMatched(1));
                                if length(matchedChannels)==0
                                    matchedChannels = [channelsEnabled(chDUTPos)];
                                else
                                    matchedChannels(length(matchedChannels)+1) = channelsEnabled(chDUTPos);
                                end
                                
                            else
                                dutTemp = [];
                            end
                            eval(['event.dut' num2str(channelsEnabled(chDUTPos)) ' = dutTemp;']);
                            
                        end
                        
                    end
                    event.matchedChannels = matchedChannels;
                    event.numberMatched = length(matchedChannels);
                    event.refSignal = refSignal;
                    nextIndex = length(matchedEvents)+1;
                    
                    if event.numberMatched>1
                        if nextIndex==1
                            matchedEvents = [event];
                        else
                            matchedEvents(nextIndex) = event;
                        end
                    end
                end
                
                
                %have matched events - go through and process events of
                %interest
                str_disp=sprintf('Processing events', i);
                disp(str_disp);
                
                for matchPos = 1:matchedEventsStep:length(matchedEvents)
                    event = matchedEvents(matchPos);
                    matchedIdx = find(event.matchedChannels==triggerSignalsChannel);
                    if 1==1 
                    %if length(matchedIdx)==1 %require trigger signals
                    %channel
                        %found trigger match signal
                                        str_disp=sprintf('Found trig channel', i);

                        %go through all matched channels and analyse for the
                        %ones under test
                        
                        for dutChannelPos = 1:length(event.matchedChannels)
                            matchedIdxToAnalyse = find(channelToAnalyse==event.matchedChannels(dutChannelPos)); %check if this should be analysed
                            if length(matchedIdxToAnalyse)==1
                                %should analyse
                                
                                %found match in DUT channel to analyse
                                eval(['dut = event.dut' num2str(event.matchedChannels(dutChannelPos)) ';']);
                                ref = event.refSignal;
                                %check match to SRS event counter
                                timeDiffTracker = scaledEventCounterTimes-ref.cell0Time;
                                [trackerMinValue,trackerMinIdx] = min(abs(timeDiffTracker));
                                trackerEntryMatched = eventCounterEntries(trackerMinIdx,:);
                                if size(trackerEntryMatched,1) %unambiguous match
                                    %found a single match in reference for this dut event -> unambiguois match
%                                     eventCounterTempLatencyArray = [eventCounterTempLatencyArray;trackerMinValue];
                                    if trackerMinValue<1e15
                                        event.eventID = trackerEntryMatched(1);
                                        %matched to tracker hit, process event
                                        
                                        dutWaveform = dut.waveform;
                                        refWaveform = ref.waveform;
                                        
                                        event_id = event.eventID;
                                        
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
                                                
                                                % store valid data into array of structures
                                                MM_data(eventPosProcessing)= MM_temp;
                                                MCP_data(eventPosProcessing)= MCP_temp;
                                                time_diff(eventPosProcessing) = MM_data(eventPosProcessing).cfd.time-MCP_data(eventPosProcessing).cfd.time;
                                                time_diff_sigmoid(eventPosProcessing) = MM_data(eventPosProcessing).sigmoid.timepoint-MCP_data(eventPosProcessing).sigmoid.timepoint;
                                                MCP_maxy(eventPosProcessing) = MCP_data(eventPosProcessing).sig.max.y;
                                                MM_maxy(eventPosProcessing) = MM_data(eventPosProcessing).sig.max.y;
                                                trackerX(eventPosProcessing) = MM_temp.x;
                                                trackerY(eventPosProcessing) = MM_temp.y;
                                                eventIDArray(eventPosProcessing) = MM_temp.event_id;
                                                dutChannelArray(eventPosProcessing) = event.matchedChannels(dutChannelPos);
                                                
                                                eventPosProcessing=eventPosProcessing+1;
                                            end
                                        end
                                        
                                    else
                                        
                                        
                                        event.eventID = 0;
                                    end
                                    %eventPos = eventPos+1;
                                else
                                    %found multiple matches
                                end
                                
                            end
                        end
                        
                       % matchedIdx = find(event.matchedChannels==channelToAnalyse);

                    end
                end
                
                %reset arrays
                signals = [];
                chIDsArray = [];
                chIDsArray = [];
                
                %finished batch of read and process files, then go to next
                str_disp=sprintf('Processed batch, total events processed: %i latestEventID: %i', eventPosProcessing,eventIDArray(eventPosProcessing-1));
                disp(str_disp);
                
                %stop processing to not crash
                %if eventPosAnalysis>450000
                %    break;
                %end
                
                %display hist of latency counter
%                 figure;
%                 hist(eventCounterTempLatencyArray,1000);
%                     xlabel('Latency SRS to signals');
%     ylabel('Counts ');
%     title(['Latency between SRS event counter and signal - Run ' run.id]);
%     grid on
% 
%                     saveas(gcf,[store_folder '\Run' run.id '_LatencySRS_batch' num2str(i) '.png'])
%                 pause(5);
% 
%                 close all
            end
            
            i = i+1;
        end
        fclose(fid);
    end
    
end

%save([store_folder '\Run' run.id '-' run.oscilloscope 'Processed.mat'], 'run', 'MM_data', 'MCP_data', 'time_diff', 'time_diff_sigmoid', 'MCP_maxy', 'MM_maxy', 'trackerX', 'trackerY', 'eventIDArray', 'dutChannelArray');

if shouldAnalyseWhenFinished
    AnalyseAllChannels
end
