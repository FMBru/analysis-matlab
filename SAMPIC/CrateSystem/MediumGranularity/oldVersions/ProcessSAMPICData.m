clearvars -except run
close all
tic

eventPosProcessing = 1; %counter of processed events

%%when combining data from multiple runs, do not reset eventPosAnalysis and
%%do not clear all

addpath '\\eosproject-smb\eos\project\p\picosec\testbeam\2025_July_h4\Matlab\CommonFunctions';
addpath 'C:\Users\GDD\Documents\MATLAB\Picosec\CommonFunctions\SCP';
addpath '\\eosproject-smb\eos\project\p\picosec\testbeam\2025_July_h4\Matlab\SAMPIC\CrateSystem\Functions';

run.id = '206';

numberHitsForLatency = 20000000;

%% reduced sampling for quick look at data - better to use refSignalsStep
matchedEventsStep = 1; %step size when looping through events to reconstruct, set 1 for all, set higher for fast analysis of reduced data set
refSignalsStep = 1; %step size when looping through ref signals to reconstruct, set 1 for all, set higher for fast analysis of reduced data set

numberFEBoards = 2;

referenceDetectorChannel = 0; %channel number of detector used as timing reference - for latency relative to this detector

%%analyse only channels listed here
channelToAnalyse = 1:110;

baseName = 'sampic_run1'; %basename used by SAMPIC for folders and files, static if using run folders
run.path=['\\eosproject-smb\eos\project\p\picosec\testbeam\2025_July_h4\SAMPIC\Run' run.id ];
storeFolderSignalRef = ['\\eosproject-smb\eos\project\p\picosec\testbeam\2025_July_h4\Results\Run' run.id '-SAMPIC\signals\REF'];
mkdir(storeFolderSignalRef);
storeFolderSignalDUT = ['\\eosproject-smb\eos\project\p\picosec\testbeam\2025_July_h4\Results\Run' run.id '-SAMPIC\signals\DUT'];
mkdir(storeFolderSignalDUT);
SAMPICFolderInfo = dir(run.path);

foundTriggerFile = false;
foundAsciiDataFile = false;
foundBinaryDataFile = false;
foundSegmentedBinaryDataFile = false;

dataFilePath = '';
triggerFilePath = '';

store_folder = ['\\eosproject-smb\eos\project\p\picosec\testbeam\2025_July_h4\Results\Run' run.id '-SAMPIC'];
mkdir(store_folder);

store_folderLatency = ['\\eosproject-smb\eos\project\p\picosec\testbeam\2025_July_h4\Results\Run' run.id '-SAMPIC\latencies'];
mkdir(store_folderLatency);

processingFolder = ['\\eosproject-smb\eos\project\p\picosec\testbeam\2025_July_h4\SAMPIC\Processing\Run' run.id ];
signalsFolder = [processingFolder '\Signals'];
mkdir(signalsFolder);


storeMatfileFolder = ['C:\Users\GDD\Documents\Picosec\Run' run.id '-SAMPIC\variables'];
storeMatfileFolder = ['F:\Processing\Picosec\Run' run.id '-SAMPIC\variables'];

processingBatchCounter = 1;
mkdir(storeMatfileFolder);


%create file to store variable

numberHeaderLines = 7; %for 64CH version
numberHeaderLines = 0; %for segemented files from crate version

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
run.year = '2025 July ';
run.name = ['BEAM ' run.year ' RUN ' run.id];
run.pathEOS=['\\eosproject-smb\eos\project\p\picosec\testbeam\2025_July_h4\' run.oscilloscope '\Run' run.id '\'];
tracker.path = ['\\eosproject-smb\eos\project\p\picosec\testbeam\2025_July_h4\tracker\reconstructed\asciiRun' run.id '.dat'];

%read tracker data from ASCII file
tracker.en = 1; %match eventIDs to tracking data and add XY to output
dutColArray = [1 4 5; 2 7 8; 3 10 11; 4 13 14; 5 16 17]; %[dID colXID colYID; ] -> May MM run

if tracker.en
    trackerFile = fopen(tracker.path,'rt');
    triggerFilePath=[run.path '\' baseName '\' baseName '_trigger_data.bin'];
    D = textscan(trackerFile, '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f', 'Delimiter','\t', 'HeaderLines',2, 'CollectOutput',1);
    tracker.data = cell2mat(D);
end


%% plot slices of tracker Data
numberTrackerEntries = length(tracker.data);
sliceCount = 1000;
numberFullSlices = floor(numberTrackerEntries/sliceCount);
close all;
for i=1:numberFullSlices
    startIndex = (i-1)*sliceCount+1;
    endIndex = i*sliceCount;

    xArray = tracker.data(startIndex:endIndex,4);
    yArray = tracker.data(startIndex:endIndex,5);
    scatter(xArray,yArray); hold on;

end
title(['Run ' run.id ' - trackerEvents: ' int2str(length(tracker.data))]);
saveas(gcf,[store_folder '\Run' run.id '_TrackerMap.png'])
pause(1);
close all;


% options for processing micromegas channel
opts_MM.t_dead = 1;      % blanking time from start of the recording samples
opts_MM.t_prerms = 5;  % blanking time before global maximum samples
opts_MM.Ts = 1/20e9;     % sampling speed
opts_MM.Rin=50;          % input impedance
opts_MM.invert = 1;      % is inverted (1-inverted, 0-non-inverted)
opts_MM.type=0;          % detector type (0-DUT, 1-REF) for future use
opts_MM.en_plot = 0;     % enable debugging plots

% options for processing MCP channel
opts_MCP.t_dead = 1;
opts_MCP.t_prerms = 5;
opts_MCP.Ts = 1/20e9;
opts_MCP.Rin=50;
opts_MCP.invert = 1;
opts_MCP.type=1;
opts_MCP.en_plot = 0;

run.savedSignals = 0;
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

triggerFilePath=[run.path '\' baseName '\' baseName '_trigger_data.bin'];

%read in run settings
runSettingsFilePath = [run.path '\' baseName '\Run_Settings.txt'];
fid = fopen(runSettingsFilePath);
while ~feof(fid)
    tline = fgets(fid);
    startString = tline(1:2);
    if startString=='=='
        %is comment
        %disp('comment');
        if contains(tline,'SamplingFrequency')
            samplingFrequencyLineSplit = split(tline);
            samplingFrequency = str2num(samplingFrequencyLineSplit{3})*1000000;
            timestep = 1/samplingFrequency;
        end
    end
end
periodFactor = (64*1000)/(samplingFrequency/1000000);

dataFilesBasePath = [run.path '\' baseName];
dataFilesList = [];

for febPos = 0:(numberFEBoards-1)
    febPath = [dataFilesBasePath '\feb' num2str(febPos)]
    febFilesList = dir(febPath);
    for i=1:length(febFilesList)
        fileInfo = febFilesList(i);
        if contains(fileInfo.name,'.bin')
            %is segemented binary data file, add to array of data files
            dataFileEntryTemp.path = [febPath '\' fileInfo.name];
            dataFileEntryTemp.feb = febPos;
            dataFilesList = [dataFilesList; dataFileEntryTemp];
        end
    end
end

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

if length(dataFilesList)>0
    str_disp=sprintf('Reading in Binary data file for latency determination');
    disp(str_disp);

    signals = [];
    chIDsArray = [];
    chTimesArray = [];

    dataFilesListInitialSamples = [];
    dataFilesListInitialSamplesFEBs = [];

    timestamp_prev = -1;
    timestamp_ov = 0;

    previousFEBID = -1;

    for posSample = 1:length(dataFilesList)
        fileEntry = dataFilesList(posSample);
        fileEntryFEB = fileEntry.feb;
        febExisting = find(dataFilesListInitialSamplesFEBs==fileEntryFEB);
        if length(febExisting)==0
            %not yet existing, add
            dataFilesListInitialSamples = [dataFilesListInitialSamples;fileEntry];
            dataFilesListInitialSamplesFEBs = [dataFilesListInitialSamplesFEBs;fileEntryFEB];
        end
    end

    for samplePos = 1:length(dataFilesListInitialSamples)
        %sample first files from each FEB
        dataFilePath=dataFilesListInitialSamples(samplePos).path;
        currentFEB = dataFilesListInitialSamples(samplePos).feb

        if previousFEBID < 0
            timestamp_prev = -1;
            timestamp_ov = 0;
        else
            if previousFEBID ~= currentFEB
                timestamp_prev = -1;
                timestamp_ov = 0;
            end
        end

        previousFEBID = currentFEB;

        fid = fopen(dataFilePath);
        i=1;
        %going through some first lines to use for latency calculationi = 1;
        while ~feof(fid) && i<numberHitsForLatency
            if mod(i,10000)==0
                toc
                tic
                str_disp=sprintf('Reading in hit No. %d', i);
                disp(str_disp);
            end
            channel = fread(fid,1,'int32',0);
            Cell0TimeStamp = fread(fid,1,'double',0);
            TOTValue = fread(fid,1,'float32',0);
            TimeInstant = fread(fid,1,'double',0);
            Baseline = fread(fid,1,'float32',0);
            PeakValue = fread(fid,1,'float32',0);
            Amplitude = fread(fid,1,'float32',0);
            dataSize = fread(fid,1,'int32',0);
            dataLength = dataSize;

            channel = channel+currentFEB*64;

            % count overflows of timestamp
            timestampRaw = Cell0TimeStamp;
            if(timestampRaw < (timestamp_prev-5000))
                timestamp_ov = timestamp_ov + 1;
            end
            timestamp_prev = timestampRaw;
            timestamp = timestamp_ov * 1099511627776*periodFactor + timestampRaw;

            sampleVector = zeros(dataLength,1);
            for pos=1:dataLength
                sampleVector(pos) = fread(fid,1,'float32',0)  ;
            end

            numberSamples = length(sampleVector);
            waveform = zeros(numberSamples,2);
            waveform(1:numberSamples,1) = (1:numberSamples)-1;
            waveform(1:numberSamples,1) = waveform(1:numberSamples,1)*timestep+timestamp*0.000000001;
            waveform(1:numberSamples,2) = sampleVector;

            hit.cell0Time = timestamp;
            hit.ch = channel;
            hit.waveform = waveform;
            hit.TOTvalue = TOTValue;
            hit.timeInstant = TimeInstant;
            hit.baseline = Baseline;
            hit.peakValue = PeakValue;
            hit.amplitude = Amplitude;

            if length(signals)==0
                signals = [hit];
                chIDsArray = [channel];
                chTimesArray = [timestamp];

            else
                if dataSize>10
                    signals (length(signals)+1) = hit;
                    chIDsArray (length(signals)) = channel;
                    chTimesArray (length(signals)) = timestamp;

                end
            end
            i = i+1;
        end
        fclose(fid);

    end

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
    xlim([0 128]);
    pause(1);
    saveas(gcf,[store_folder '\Run' run.id '_CHHits.png'])
    %plot hitmap accross multipad
    VisualiseMultipad(channelsEnabled,chCountsArray,'Hitmap Multipad Picosec Relative',[store_folder '\Run' run.id '_HitmapRelative.png'],'Hits','%0.1f %',0,0);
    VisualiseMultipad(channelsEnabled,chCountsArrayNumbers,'Hitmap Multipad Picosec',[store_folder '\Run' run.id '_Hitmap.png'],'Hits','%0.1f %',0,0);

    %%determine latency for channels relative to reference detector
    %refIndexInEnabledChannels = find(channelsEnabled==referenceDetectorChannel);
    refChannelMask = chIDsArray==referenceDetectorChannel;
    refTime = chTimesArray(refChannelMask);

    %% loop through DUT channels  - get latency
    for chPos = 1:length(channelsEnabled)
        if  channelsEnabled(chPos) == referenceDetectorChannel
            %is ref, get latency for SRS to ref

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

            figure
            hold on
            hist(timeDiffToRefArrayCut,1000);
            xlabel('Min time between SRS counter and REF (ns)');
            ylabel('Events');
            title(['SAMPIC - Minimum time between EventCounter to RefChannel:  - Run ' run.id]);
            grid on
            xlim([0 100000]);
            saveas(gcf,[store_folder '\Run' run.id '_latency_counter_scaling.png'])
            %close all
            pause(1);

            [refTimeCounts,refTimeCenters] = hist(refTime,10000);
            [eventCounterCounts, eventCounterCenters] = hist(scaledEventCounterTimes,10000);

            close all

            figure
            hold on
            plot(eventCounterCenters,eventCounterCounts);
            plot(refTimeCenters,refTimeCounts);

            if isnan(latency)
            else
                xline(latency,'color','red','linewidth',2)
            end
            xlabel('Time');
            ylabel('Events');
            title(['SAMPIC - Minimum time between EventCounter and Ref - Run ' run.id]);
            legend('Event Counter time','Ref hits time')
            grid on

            saveas(gcf,[store_folder '\Run' run.id '_timeComparison.png'])
            pause(1);
            close all
        else
            %is dut channel

            %check if channel included in analysis
            shouldAnalyseChannel = find(channelToAnalyse==channelsEnabled(chPos));
            if length(shouldAnalyseChannel)==1

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
                    timeDiffToRef = detectorHitTime-refTime;
                    minTimeDiffToRef = min(abs(timeDiffToRef));
                    if minTimeDiffToRef<200 && minTimeDiffToRef~0;
                        timeDiffToRefArray = [timeDiffToRefArray;minTimeDiffToRef];
                    end
                end

                %determine latency from min time diff histogram between ch and ref
                latency = median (timeDiffToRefArray);

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
                saveas(gcf,[store_folderLatency '\Run' run.id '_latency_CH' num2str(channelsEnabled(chPos)) '.png']);
                pause(1);
                close all;

                eval(['latency_ch' num2str(channelsEnabled(chPos)) ' = latency;']);
                chLatencyArray(chPos,1) = latency;
            end
        end

    end
    figure
    plot(channelsEnabled,chLatencyArray,'.');
    xlabel('Channel ID');
    ylabel('Latency to ref (ns) ');
    title(['Latency to Ref SAMPIC channel - Run ' run.id]);
    grid on
    ylim([0 50]);
    xlim([0 128]);
    pause(1);
    saveas(gcf,[store_folderLatency '\Run' run.id '_CHLatencyToRef.png'])

    %plot latency accross multipad
    VisualiseMultipad(channelsEnabled,abs(chLatencyArray),'Latency Multipad Picosec',[store_folder '\Run' run.id '_Latency.png'],'Latency','%0.1f ns',0,0);

    pause(1);
    close all;
    %% latencies determined - process binary file in segments

    %for timestamp overflows
    timestamp_prev = -1;
    timestamp_ov = 0;

    previousFEBID = -1;

    str_disp=sprintf('Reading in Binary data file for processing');
    disp(str_disp);


    %% loop through files to extract ref signals
    signalsRef = [];
    chIDsRefArray = [];
    chTimesRefArray = [];
    chIndexRefArray = [];

    for filePos=1:length(dataFilesList)

        str_disp=sprintf('Ref signals extraction - starting Binary File No. %d', filePos);
        disp(str_disp);
        currentFilePath = dataFilesList(filePos).path;
        currentFEB = dataFilesList(filePos).feb


        if previousFEBID < 0
            timestamp_prev = -1;
            timestamp_ov = 0;
        else
            if previousFEBID ~= currentFEB
                timestamp_prev = -1;
                timestamp_ov = 0;
            end
        end

        previousFEBID = currentFEB;

        fid = fopen(currentFilePath);

        for i=1:numberHeaderLines
            tline = fgets(fid);
        end
        i = 1;

        while ~feof(fid) %&& eventPosAnalysis<500000 %process full file
            if mod(i,10000)==0
                toc
                tic
            end
            channel = fread(fid,1,'int32',0);
            Cell0TimeStamp = fread(fid,1,'double',0);
            TOTValue = fread(fid,1,'float32',0);
            TimeInstant = fread(fid,1,'double',0);
            Baseline = fread(fid,1,'float32',0);
            PeakValue = fread(fid,1,'float32',0);
            Amplitude = fread(fid,1,'float32',0);
            dataSize = fread(fid,1,'int32',0);
            dataLength = dataSize;

            channel = channel+currentFEB*64;

            % count overflows of timestamp
            timestampRaw = Cell0TimeStamp;
            if(timestampRaw < (timestamp_prev-5000))
                timestamp_ov = timestamp_ov + 1;
            end
            timestamp_prev = timestampRaw;
            timestamp = timestamp_ov * 1099511627776*periodFactor + timestampRaw;


            sampleVector = zeros(dataLength,1);
            for pos=1:dataLength
                sampleVector(pos) = fread(fid,1,'float32',0)  ;
            end

            numberSamples = length(sampleVector);
            waveform = zeros(numberSamples,2);
            waveform(1:numberSamples,1) = (1:numberSamples)-1;
            waveform(1:numberSamples,1) = waveform(1:numberSamples,1)*timestep+timestamp*0.000000001;
            waveform(1:numberSamples,2) = sampleVector;

            hit.cell0Time = timestamp;
            %hit.unixTimestamp = epochTime; %or hitUnixTime
            hit.ch = channel;
            %hit.hit = hitNumber;
            hit.waveform = waveform;
            hit.TOTvalue = TOTValue;
            hit.timeInstant = TimeInstant;
            hit.baseline = Baseline;
            hit.peakValue = PeakValue;
            hit.amplitude = Amplitude;

            if length(channel)>0
                if channel == referenceDetectorChannel
                    if length(signalsRef)==0
                        signalsRef = [hit];
                        chIDsRefArray = [channel];
                        chTimesRefArray = [timestamp];
                        chIndexRefArray = [1];
                    else
                        if dataSize>10
                            signalsRef (length(signalsRef)+1) = hit;
                            chIDsRefArray (length(signalsRef)) = channel;
                            chTimesRefArray (length(signalsRef)) = timestamp;
                            chIndexRefArray (length(signalsRef)) = length(signalsRef);
                        end
                    end
                end
            end
            eventCounterTempLatencyArray = [];
            minTimeDiffs = [];
            i = i+1;
        end
        fclose(fid);
    end

    str_disp=sprintf('Finished extracing ref signals');
    disp(str_disp);

    timestamp_prev = -1;
    timestamp_ov = 0;

    previousFEBID = -1;

    %% loop through files for processing
    for filePos=1:length(dataFilesList)

        str_disp=sprintf('Starting Binary File No. %d', filePos);
        disp(str_disp);


        currentFilePath = dataFilesList(filePos).path;
        currentFEB = dataFilesList(filePos).feb

        if previousFEBID < 0
            timestamp_prev = -1;
            timestamp_ov = 0;
        else
            if previousFEBID ~= currentFEB
                timestamp_prev = -1;
                timestamp_ov = 0;
            end

        end

        previousFEBID = currentFEB;

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
            channel = fread(fid,1,'int32',0);
            Cell0TimeStamp = fread(fid,1,'double',0);
            TOTValue = fread(fid,1,'float32',0);
            TimeInstant = fread(fid,1,'double',0);
            Baseline = fread(fid,1,'float32',0);
            PeakValue = fread(fid,1,'float32',0);
            Amplitude = fread(fid,1,'float32',0);
            dataSize = fread(fid,1,'int32',0);
            dataLength = dataSize;

            channel = channel+currentFEB*64;

            % count overflows of timestamp
            timestampRaw = Cell0TimeStamp;
            if(timestampRaw < (timestamp_prev-5000))
                timestamp_ov = timestamp_ov + 1;
            end
            timestamp_prev = timestampRaw;
            timestamp = timestamp_ov * 1099511627776*periodFactor + timestampRaw;

            sampleVector = zeros(dataLength,1);
            for pos=1:dataLength
                sampleVector(pos) = fread(fid,1,'float32',0)  ;
            end

            numberSamples = length(sampleVector);
            waveform = zeros(numberSamples,2);
            waveform(1:numberSamples,1) = (1:numberSamples)-1;
            waveform(1:numberSamples,1) = waveform(1:numberSamples,1)*timestep+timestamp*0.000000001;
            waveform(1:numberSamples,2) = sampleVector;

            hit.cell0Time = timestamp;
            hit.ch = channel;
            hit.waveform = waveform;
            hit.TOTvalue = TOTValue;
            hit.timeInstant = TimeInstant;
            hit.baseline = Baseline;
            hit.peakValue = PeakValue;
            hit.amplitude = Amplitude;

            if length(channel)>0
                shouldAnalyseChannel = find(channel==channelToAnalyse);
                %do not accept ref signals - already saved
                if length(shouldAnalyseChannel)==1 && channel ~= referenceDetectorChannel

                    if length(signals)==0
                        signals = [hit];
                        chIDsArray = [channel];
                        chTimesArray = [timestamp];
                        chIndexArray = [1];
                    else
                        if dataSize>10
                            signals (length(signals)+1) = hit;
                            chIDsArray (length(signals)) = channel;
                            chTimesArray (length(signals)) = timestamp;
                            chIndexArray (length(signals)) = length(signals);
                        end
                    end
                end

            end
            eventCounterTempLatencyArray = [];
            minTimeDiffs = [];


            if length(signals)>5000
                %% match and process signals
                str_disp=sprintf('Matching events', i);
                disp(str_disp);
                matchedEvents = [];

                minSignalsTime = min(chTimesArray);
                maxSignalsTime = max(chTimesArray);

                chTimesArrayRef = chTimesRefArray;
                chIndexArrayRef = chIndexRefArray;
                refSignalsSubsetMask = chTimesArrayRef>minSignalsTime & chTimesArrayRef<maxSignalsTime;

                signalsRefSubset = signalsRef(refSignalsSubsetMask);

                %go through all ref signals - previously determined

                str_disp=sprintf('%d ref signals in subset', length(signalsRefSubset));
                disp(str_disp);

                for refPos = 1:refSignalsStep:length(signalsRefSubset)
                    refSignal = signalsRefSubset(refPos);
                    if mod(refPos,5000)==0
                        str_disp=sprintf('Matched ref signal %d', refPos);
                        disp(str_disp);
                    end

                    %go through all channels and check if there is a match
                    matchedChannels = [];

                    for chDUTPos = 1:length(channelsEnabled)
                        if  channelsEnabled(chDUTPos) ~= referenceDetectorChannel

                            shouldAnalyseChannel = find(channelsEnabled(chDUTPos)==channelToAnalyse);
                            % shouldAnalyseChannel=1;
                            if length(shouldAnalyseChannel)==1

                                %is dut channel
                                latency = chLatencyArray(chDUTPos);
                                dutChannelMask = chIDsArray==channelsEnabled(chDUTPos);
                                chTimesArrayDUT = chTimesArray(dutChannelMask);
                                chIndexArrayDUT = chIndexArray(dutChannelMask);

                                minTime = latency-matchingDUTREFWindowWidth;
                                maxTime = latency+matchingDUTREFWindowWidth;
                                timeDiffArray = abs(chTimesArrayDUT-refSignal.cell0Time);

                                minTimeDiff = min(timeDiffArray);
                                %minTimeDiffs = [minTimeDiffs;minTimeDiff];

                                matchMask = timeDiffArray>minTime & timeDiffArray<maxTime;

                                chIndexMatched = chIndexArrayDUT(matchMask);

                                if length(chIndexMatched)==1 %unambiguous hit
                                    %accept as match
                                    length(chIndexMatched);
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

                    end
                    event.matchedChannels = matchedChannels;
                    event.numberMatched = length(matchedChannels);
                    event.refSignal = refSignal;
                    nextIndex = length(matchedEvents)+1;

                    matchedEvents = [matchedEvents event];
                end

                %have matched events - go through and process events of
                %interest
                str_disp=sprintf('Processing events', i);
                disp(str_disp);

                for matchPos = 1:matchedEventsStep:length(matchedEvents)
                    event = matchedEvents(matchPos);
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

                                    run.savedSignals = run.savedSignals +1;
                                    if run.savedSignals<10
                                        shouldSave = 1;
                                    else
                                        shouldSave = 0;
                                    end

                                    % process MM Picosec first to see if signal is valid
                                    MM_temp = process_signal_sampic(t_vec_mm-etime,dutWaveform(:,2),opts_MM,1,run,shouldSave,storeFolderSignalDUT);

                                    % if signal is valid process MCP
                                    if(MM_temp.fail==0)
                                        MCP_temp = process_signal_sampic(t_vec_mcp-etime,refWaveform(:,2),opts_MCP,1,run,shouldSave,storeFolderSignalRef);

                                        % if MCP is valid store data to structure array
                                        if(MCP_temp.fail==0)
                                            % process tracker ID channel

                                            MM_temp.event_id = event_id;
                                            MCP_temp.event_id = event_id;

                                            event_id;
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

                                            %save data directly extracted sampic file

                                            MM_temp.sampic.tot = dut.TOTvalue;
                                            MM_temp.sampic.amp = dut.amplitude;

                                            MCP_temp.sampic.tot = ref.TOTvalue;
                                            MCP_temp.sampic.amp = ref.amplitude;

                                            % store valid data into array of structures
                                            MM_data(eventPosProcessing)= MM_temp;
                                            MCP_data(eventPosProcessing)= MCP_temp;
                                            time_diff(eventPosProcessing) = MM_data(eventPosProcessing).cfd.time-MCP_data(eventPosProcessing).cfd.time;
                                            time_diff_sigmoid(eventPosProcessing) = MM_data(eventPosProcessing).sigmoid.timepoint-MCP_data(eventPosProcessing).sigmoid.timepoint;
                                            MCP_maxy(eventPosProcessing) = MCP_data(eventPosProcessing).sig.max.y;
                                            MM_maxy(eventPosProcessing) = MM_data(eventPosProcessing).sig.max.y;
                                            MM_bgAvg(eventPosProcessing) = MM_data(eventPosProcessing).sig.blavg;
                                            MM_bgRMS(eventPosProcessing) = MM_data(eventPosProcessing).sig.blrms;
                                            trackerX(eventPosProcessing) = MM_temp.x;
                                            trackerY(eventPosProcessing) = MM_temp.y;
                                            eventIDArray(eventPosProcessing) = MM_temp.event_id;
                                            dutChannelArray(eventPosProcessing) = event.matchedChannels(dutChannelPos);
                                            currentChannel = event.matchedChannels(dutChannelPos);

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

                end
                signals;
                %% %reset arrays
                signals = [];
                chIDsArray = [];
                chTimesArray = [];
                chIndexArray = [];


                %finished batch of read and process files, then go to next
                str_disp=sprintf('Processed batch, total events analysed: %i latestEventID: %i', eventPosProcessing,eventIDArray(eventPosProcessing-1));
                disp(str_disp);

                str_disp=sprintf('Saving to Matfile');
                disp(str_disp);

                storeMatfilePath = [storeMatfileFolder '\processedBatch_' int2str(processingBatchCounter) '.mat'];
                m = matfile(storeMatfilePath,'Writable',true);

                save(storeMatfilePath,'MM_data');
                save(storeMatfilePath,'MCP_data','-append');
                save(storeMatfilePath,'time_diff','-append');
                save(storeMatfilePath,'time_diff_sigmoid','-append');
                save(storeMatfilePath,'MCP_maxy','-append');
                save(storeMatfilePath,'MM_maxy','-append');
                save(storeMatfilePath,'trackerX','-append');
                save(storeMatfilePath,'trackerY','-append');
                save(storeMatfilePath,'eventIDArray','-append');
                save(storeMatfilePath,'dutChannelArray','-append');
                save(storeMatfilePath,'currentChannel','-append');
                save(storeMatfilePath,'MM_bgAvg','-append');
                save(storeMatfilePath,'MM_bgRMS','-append');


                clear MM_data MCP_data time_diff time_diff_sigmoid MCP_maxy MM_maxy trackerX trackerY eventIDArray dutChannelArray currentChannel MM_bgAvg MM_bgRMS;

                eventPosProcessing = 1;
                processingBatchCounter = processingBatchCounter+1;

                str_disp=sprintf('Saved variables');
                disp(str_disp);

            end

            i = i+1;
        end
        fclose(fid);
    end

end

loadVariablesFromFiles

