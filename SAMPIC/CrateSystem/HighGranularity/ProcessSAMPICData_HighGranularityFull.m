clear all
close all
tic

eventPosProcessing = 1; %counter of processed events
eventIDArray = [];

%%when combining data from multiple runs, do not reset eventPosAnalysis and
%%do not clear all

addpath 'C:\Users\gdd.CERN\Documents\MATLAB\Picosec\CommonFunctions';
addpath 'C:\Users\gdd.CERN\Documents\MATLAB\Picosec\CommonFunctions\SCP';
addpath 'C:\Users\gdd.CERN\Documents\MATLAB\Picosec\SAMPIC\CrateSystem\Functions';
addpath '\\eosproject-smb\eos\project\p\picosec\testbeam\2024_September_h4\Matlab\SAMPIC\CrateSystem\HighGranularity\Functions';
addpath '\\eosproject-smb\eos\project\p\picosec\testbeam\2024_September_h4\Matlab\SAMPIC\CrateSystem\Functions';
addpath '\\eosproject-smb\eos\project\p\picosec\testbeam\2024_September_h4\Matlab\CommonFunctions';

hitsPerSegment = 30000; %how many hits to save per segment
numberHitsForLatency = 30000;
processingSegmentSize = 50000;

numberHitsForLatency = 20000000;

%% reduced sampling for quick look at data - better to use refSignalsStep
matchedEventsStep = 1; %step size when looping through events to reconstruct, set 1 for all, set higher for fast analysis of reduced data set
refSignalsStep = 1; %step size when looping through ref signals to reconstruct, set 1 for all, set higher for fast analysis of reduced data set

numberFEBoards = 2;

run.id = '315';
referenceDetectorChannel = 64; %channel number of detector used as timing reference - for latency relative to this detector
triggerSignalsChannel = 57; %not usinbg anymore, channel number to which trigger signal is connected which is sent to SRS for sync with tracker

%%analyse only channels listed here
channelToAnalyse = [64;65;66;67;68;69;70;71;72;73;74;75;76;77;78;79;80;81;82;83;84;85;86;87;88;89;90;91;92;93;94;95;96;97;98;99;100;101;0;1];

baseName = 'sampic_run1'; %basename used by SAMPIC for folders and files, static if using run folders
run.path=['\\eosproject-smb\eos\project\p\picosec\testbeam\2024_September_h4\SAMPIC\Run' run.id ];


storeFolderSignalRef = ['\\eosproject-smb\eos\project\p\picosec\testbeam\2024_September_h4\Results\Run' run.id '-SAMPIC\signals\REF'];
mkdir(storeFolderSignalRef);
storeFolderSignalDUT = ['\\eosproject-smb\eos\project\p\picosec\testbeam\2024_September_h4\Results\Run' run.id '-SAMPIC\signals\DUT'];
mkdir(storeFolderSignalDUT);

SAMPICFolderInfo = dir(run.path);

foundTriggerFile = false;
foundAsciiDataFile = false;
foundBinaryDataFile = false;
foundSegmentedBinaryDataFile = false;

dataFilePath = '';
triggerFilePath = '';

store_folder = ['\\eosproject-smb\eos\project\p\picosec\testbeam\2024_September_h4\Results\Run' run.id '-SAMPIC'];
mkdir(store_folder);

store_folderLatency = ['\\eosproject-smb\eos\project\p\picosec\testbeam\2024_September_h4\Results\Run' run.id '-SAMPIC\latencies'];
mkdir(store_folderLatency);

processingFolder = ['\\eosproject-smb\eos\project\p\picosec\testbeam\2024_September_h4\SAMPIC\Processing\Run' run.id ];
signalsFolder = [processingFolder '\Signals'];
mkdir(signalsFolder);

storeMatfileFolder = ['\\eosproject-smb\eos\project\p\picosec\testbeam\2024_September_h4\SAMPIC\Processing\Run' run.id '-SAMPIC\variables'];
storeMatfileFolder = ['C:\Users\gdd.CERN\Documents\Picosec\Run' run.id '-SAMPIC\variables'];
processingBatchCounter = 1;
mkdir(storeMatfileFolder);


%create file to store variable

numberHeaderLines = 7; %for 64CH version
numberHeaderLines = 0; %for segemented files from crate version

matchingDUTREFWindowWidth=3;
%matchingDUTREFWindowWidth=10;

channelsEnabled = []; %ids of enabled channels


%load mapping
setMappingHighGranularity



%% config for processing


run.oscilloscope = 'SAMPIC';
opts_MM.ch_name = ['DUT CH'];
opts_MCP.ch_name = ['C' num2str(referenceDetectorChannel)]; % MPC 2 (used as trigger)
shouldSaveMAT = false;
shouldUseEOSFolder = true;
tracker.dutIndex = 2; 

% run file coordinates and number of files
minuit = 1;   % set minuit for optimizer if 1 or matlab 0
run.year = '2024 September ';
run.name = ['BEAM ' run.year ' RUN ' run.id];
run.pathEOS=['\\eosproject-smb\eos\project\p\picosec\testbeam\2024_September_h4\' run.oscilloscope '\Run' run.id '\'];

if shouldUseEOSFolder
    %tracker.path = ['\\eosproject-smb\eos\project\p\picosec\testbeam\2022_October_h4\tracker\reconstructed\asciiRun' run.id '.dat'];
    tracker.path = ['\\eosproject-smb\eos\project\p\picosec\testbeam\2024_September_h4\tracker\reconstructed\asciiRun' run.id '.dat'];
else
    tracker.path = ['C:\Users\gdd.CERN\Documents\Picosec\May22\Tracker\reconstructed\asciiRun' run.id '.dat'];
end

%override tracker file path
%tracker.path = ['\\eosproject-smb\eos\project\p\picosec\testbeam\2024_June_h4\tracker\reconstructed\asciiRun232.dat'];

%read tracker data from ASCII file
tracker.en = 1; %match eventIDs to tracking data and add XY to output

dutColArray = [1 7 8; 2 10 11; 3 13 14; 4 16 17; 5 22 23; 6 25 26; 7 16 17]; %[dID colXID colYID; ] 

tracker.path

if tracker.en
    trackerFile = fopen(tracker.path,'rt');
    triggerFilePath=[run.path '\' baseName '\' baseName '_trigger_data.bin'];
    D = textscan(trackerFile, '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f', 'Delimiter','\t', 'HeaderLines',2, 'CollectOutput',1);
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
        
        if contains(tline,'SamplingFrequency')
            samplingFrequencyLineSplit = split(tline);
            samplingFrequency = str2num(samplingFrequencyLineSplit{3})*1000000;
            timestep = 1/samplingFrequency;
        end
    else
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



    pause(1)

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
        samplePos
        dataFilePath
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
        %going through some first lines to use for latency calculation i = 1;
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
                if dataSize>20
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

    padIDArray = [];
    for chPos = 1:length(channelsEnabled)
        chID = channelsEnabled(chPos);
        chIDsMask = chIDsArray==chID;
        chCountsArray(chPos,1) = length(chIDsArray(chIDsMask));
       padID = getPadForChannelNumberHighGranularity(chID);
       padIDArray = [padIDArray;padID];
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
    VisualiseHighGranularity(channelsEnabled,chCountsArray,'Hitmap Multipad Picosec Relative',[store_folder '\Run' run.id '_HitmapRelative.png'],'Hits','%0.1f %',0,0);
    VisualiseHighGranularity(channelsEnabled,chCountsArrayNumbers,'Hitmap Multipad Picosec',[store_folder '\Run' run.id '_Hitmap.png'],'Hits','%0.1f %',2000,6000);

hitArray = [channelsEnabled chCountsArrayNumbers padIDArray]

    %%determine latency for channels relative to reference detector
    refChannelMask = chIDsArray==referenceDetectorChannel;
    refTime = chTimesArray(refChannelMask);

    %% loop through DUT channels  - get latency
    for chPos = 1:length(channelsEnabled)
        if  channelsEnabled(chPos) == referenceDetectorChannel
            %is ref, get latency for SRS to ref
            timestamps = [];
            timeDiffToRefArray = [];

            for k = 1:length(refTime)
                % Loop through all events in DUT to find time differences to
                refSignal = refTime(k);
                timeDifferencesDUTRef = [];
                detectorHitTime = refSignal;
                timeDiffToRef = detectorHitTime-scaledEventCounterTimes;
                minTimeDiffToRef = min(abs(timeDiffToRef));
                timeDiffToRefArray = [timeDiffToRefArray;minTimeDiffToRef];
            end

            timeDiffToRefArrayMask = timeDiffToRefArray<1e6;
            timeDiffToRefArrayCut = timeDiffToRefArray(timeDiffToRefArrayMask);

            % Determine latency from min time diff histogram between ch and ref
            latency = median (timeDiffToRefArrayCut);
            eventCounterLatency = latency;

%             figure
%             hold on
%             hist(timeDiffToRefArrayCut,1000);
% %             if isnan(latency)
% % 
% %             else
% %                 %    xline(latency,'color','red','linewidth',2)
% %             end
%             xlabel('Min time between SRS counter and REF (ns)');
%             ylabel('Events');
%             title(['SAMPIC - Minimum time between EventCounter to RefChannel:  - Run ' run.id]);
%             grid on
%             xlim([0 100000]);
%             saveas(gcf,[store_folder '\Run' run.id '_latency_counter_scaling.png'])
%             %close all
%             pause(1);

            [refTimeCounts,refTimeCenters] = hist(refTime,10000);
            [eventCounterCounts, eventCounterCenters] = hist(scaledEventCounterTimes,10000);

            close all

%             figure
%             hold on
%             plot(eventCounterCenters,eventCounterCounts);
%             plot(refTimeCenters,refTimeCounts);
% 
%             if isnan(latency)
% 
%             else
%                 xline(latency,'color','red','linewidth',2)
%             end
%             xlabel('Time');
%             ylabel('Events');
%             title(['SAMPIC - Minimum time between EventCounter and Ref - Run ' run.id]);
%             legend('Event Counter time','Ref hits time')
%             grid on
% 
%             saveas(gcf,[store_folder '\Run' run.id '_timeComparison.png'])
%             %close all
%             pause(1);
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
                    dutWaveform = dutSignal.waveform;
                    dutBaseline = mean(dutWaveform(1:5,2));
                    dutWaveform(:,2) = -(dutWaveform(:,2)-dutBaseline);

                    timeDifferencesDUTRef = [];
                    detectorHitTime = dutSignal.cell0Time;
                    signalAmplitude = max(dutWaveform(:,2));
                    %timestampEntry = [detectorHitTime detectorHitTimeUnix];
                    % timestamps = [timestamps; timestampEntry];
                    timeDiffToRef = detectorHitTime-refTime;
                    minTimeDiffToRef = min(abs(timeDiffToRef));
                    if minTimeDiffToRef<200 && minTimeDiffToRef~0 && signalAmplitude>0.05;
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
                saveas(gcf,[store_folderLatency '\Run' run.id '_latency_CH' num2str(channelsEnabled(chPos)) '.png']);

                if channelsEnabled(chPos)==27
                    pause(5);
                end

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
    VisualiseHighGranularity(channelsEnabled,abs(chLatencyArray),'Latency Multipad Picosec',[store_folder '\Run' run.id '_Latency.png'],'Latency','%0.1f ns',0,0);

    pause(1);
    close all;
    %% latencies determined - process binary file in segments

    %for timestamp overflows
    timestamp_prev = -1;
    timestamp_ov = 0;

    previousFEBID = -1;

    str_disp=sprintf('Reading in Binary data file for processing');
    disp(str_disp);

    dataFilesList;


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
                shouldAnalyseChannel = find(channel==channelToAnalyse);
                %shouldAnalyseChannel
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

    timestamp_prev = -1;
    timestamp_ov = 0;

    previousFEBID = -1;

%%extracted ref times - check distribution of differences
refTimeDiffArray = [];

                for refPos = 2:length(chTimesRefArray)
                    timeDiff = chTimesRefArray(refPos)-chTimesRefArray(refPos-1);
                    if timeDiff<1e9
                        refTimeDiffArray = [refTimeDiffArray;timeDiff];
                    end
                end

        signals = [];
        chIDsArray = [];
        chTimesArray = [];
        chIndexArray = [];


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


        % currentFilePath = dataFilePath;

        fid = fopen(currentFilePath);

        for i=1:numberHeaderLines
            tline = fgets(fid);
        end

        % Going through all files for data
        i = 1;


        while ~feof(fid) %&& eventPosAnalysis<500000 %process full file
            if mod(i,10000)==0
                toc
                tic
                str_disp=sprintf('Reading in hit No. %d', i);
                disp(str_disp);
            end
            % hitNumber = fread(fid,1,'int32',0);
            % epochTime = fread(fid,1,'double',0);
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

            % Count overflows of timestamp
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
                shouldAnalyseChannel = find(channel==channelToAnalyse);
                %do not accept ref signals - already saved
                if length(shouldAnalyseChannel)==1 % | channel == referenceDetectorChannel

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

            i = i+1;
        end
        fclose(fid);
    end

end



%save signals
storeMatfilePath = [storeMatfileFolder '\signals_' int2str(processingBatchCounter) '.mat'];
 mkdir(storeMatfileFolder)
m = matfile(storeMatfilePath,'Writable',true);
%
save(storeMatfilePath,'signals');
save(storeMatfilePath,'chIDsArray','-append');
save(storeMatfilePath,'chTimesArray','-append');
save(storeMatfilePath,'chIndexArray','-append');



processSignals

loadVariablesFromFiles_HighGranularityMultiplePads

