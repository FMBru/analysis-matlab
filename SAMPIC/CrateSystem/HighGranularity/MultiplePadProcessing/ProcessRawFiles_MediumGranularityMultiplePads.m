% MATLAB/OCTAVE implementation of processing code for PICOSEC
% version date: 7 Nov 2021

close all
addpath 'C:\Users\GDD\Documents\MATLAB\Picosec\CommonFunctions';

enDebugPlots = true;
correctFilenumberOffset = false;
fileStepSize = 1; %1 analyse all, or larger step for quick monitoring
%read tracker data from ASCII file
tracker.en = 1; %match eventIDs to tracking data and add XY to output

numberEventsPerBatch = 1000;

%May22: %1: miniCactus, 2: MM3 (VacChamberPos), 3: MM1(Multipad), 4: MM2 (on support plate), 5: MM4 (electron setup)

    run.id = runIDString;
    run.pad = str2num(padIDString);
    run.oscilloscope = oscilloscopeString;
    run.lecroy_name = '--Trace--'; %['Run' run.id];
    opts_MM.chID = DUTnameString;
    opts_MCP.chID = MCPnameString;
    opts_MM.ch_name = ['C' DUTnameString run.lecroy_name];
    opts_MCP.ch_name = ['C' MCPnameString run.lecroy_name];
    shouldSaveMAT = false;
    shouldUseEOSFolder = useEOSFolder;
    tracker.dutIndex = str2num(trackerDUTIndex);  %DUT index: 1: MM1 (Multipad), 2: MM2(VacChamber), 3: MM3 (ElectronsetupMM), 4: MM4 (Picolarge main tracker), 5:MM5 (SiPM)
    numberFilesToAnalyse = noFileToAnalyze;


analysis.dutChannel = opts_MM.chID;

numberDebugPlots = 1;
saveSignalWaveformsNumber =2;
debugPlotCounter = 0;
run.savedSignals = 0;

%correct for mistake in GDD scope where C2 files are saved starting at
%00001 instead of 00000
shouldCorrectFileCounter_GDD_C2 = false;

% run file coordinates and number of files
minuit = 1;   % set minuit for optimizer if 1 or matlab 0

run.year = '2023 Apr ';
run.name = ['BEAM ' run.year ' RUN ' run.id];
run.pathEOS=['\\eosproject-smb\eos\project\p\picosec\testbeam\2023_April_h4\' run.oscilloscope '\Run' run.id '\'];
disp(run.pathEOS)
run.nfiles = find_fileNo(run.pathEOS);

tracker.path = ['\\eosproject-smb\eos\project\p\picosec\testbeam\2023_April_h4\tracker\reconstructed\asciiRun' run.id '.dat'];

if(numberFilesToAnalyse>0)
    run.nfiles =numberFilesToAnalyse;
end
    
run.lecroy_name = '--Trace--'; %['Run' run.id];

run.debugPlotsPath=['\\eosproject-smb\eos\project\p\picosec\testbeam\2023_April_h4\Results\DebugPlots\Run' run.id '-'  run.oscilloscope '\'];
mkdir(run.debugPlotsPath);

storeFolderSignalRef = ['\\eosproject-smb\eos\project\p\picosec\testbeam\2023_April_h4\Results\Run' run.id '-' run.oscilloscope '-dut' opts_MM.chID '-ref' opts_MCP.chID '\signals\REF'];
mkdir(storeFolderSignalRef);
storeFolderSignalDUT = ['\\eosproject-smb\eos\project\p\picosec\testbeam\2023_April_h4\Results\Run' run.id '-' run.oscilloscope '-dut' opts_MM.chID '-ref' opts_MCP.chID '\signals\DUT'];
mkdir(storeFolderSignalDUT);

storeMatfileFolder = ['\\eosproject-smb\eos\project\p\picosec\testbeam\2023_April_h4\Results\Run' run.id '-Picolarge\variables'];
processingBatchCounter = 1;
mkdir(storeMatfileFolder);



dutColArray = [1 10 11; 2 16 17; 3 4 5; 4 13 14; 5 7 8]; %[dID colXID colYID; ] -> July 2022 MM run
%DUT index: 1: MM1 (Multipad), 2: MM2(VacChamber), 3: MM3
%(ElectronsetupMM), 4: MM4 (Picolarge main tracker), 5:MM5 (SiPM) July 2022

if tracker.en
    trackerFile = fopen(tracker.path,'rt');
    D = textscan(trackerFile, '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f', 'Delimiter','\t', 'HeaderLines',2, 'CollectOutput',1);
    tracker.data = cell2mat(D);
end

% options for processing micromegas channel (lot magic numbers from
% processing functionhave to be added here)
opts_MM.t_dead = 1;      % blanking time from start of the recording samples
opts_MM.t_prerms = 20;  % blanking time before global maximum samples 
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
opts_MCP.t_prerms = 20;
opts_MCP.Ts = 1/20e9;
opts_MCP.Rin=50;
opts_MCP.invert = 1;
opts_MCP.type=1;
%only for override - should set in top
%opts_MCP.ch_name = ['C1' run.lecroy_name]; % MPC 2 (used as trigger)
opts_MCP.en_plot = 0;

% options for extracting tracker ID from the bitstream
opts_TR.ch_name = ['C3' run.lecroy_name];
opts_TR.baud_rate = 40e6; % 40Mbps baudrate
opts_TR.n_bits = 16;      % number of bits after start bit

%% DO PROCCEISNG
tic;         % start time measurement
k=1;         % valid data counter
j=1;         % all data counter
event_id_prev = -1;
event_id_ov = 0;
monitorPlotCounter = 0;

for (ff=1:fileStepSize:run.nfiles)
    % display current file number
    str_disp=sprintf('Loading file set No. %d', ff);
    disp(str_disp);
    
    % generate file name strings
    
    if correctFilenumberOffset
        %file numbers are shifted
        ch_mm_str=sprintf('%s%s%05d.trc', run.pathEOS, opts_MM.ch_name,ff);
        ch_mcp_str=sprintf('%s%s%05d.trc', run.pathEOS, opts_MCP.ch_name,ff+1); %add ff+1 if file numbers are shifted
        ch_tr_str=sprintf('%s%s%05d.trc', run.pathEOS, opts_TR.ch_name,ff+1); %add ff+1 if file numbers are shifted
    else
        %no problem with filenumbers
        ch_mm_str=sprintf('%s%s%05d.trc', run.pathEOS, opts_MM.ch_name,ff);
        ch_mcp_str=sprintf('%s%s%05d.trc', run.pathEOS, opts_MCP.ch_name,ff); 
        ch_tr_str=sprintf('%s%s%05d.trc', run.pathEOS, opts_TR.ch_name,ff);
    end
    
    if shouldCorrectFileCounter_GDD_C2
        if strcmp(run.oscilloscope,'GDD')
            %shift counter for C2 by one
            if strcmp(opts_MM.ch_name,['C2' run.lecroy_name])
                ch_mm_str=sprintf('%s%s%05d.trc', run.path, opts_MM.ch_name,(ff+1));
                if shouldUseEOSFolder
                    ch_mm_str=sprintf('%s%s%05d.trc', run.pathEOS, opts_MM.ch_name,(ff+1));
                end
            end
            
            if strcmp(opts_MCP.ch_name,['C2' run.lecroy_name])
                ch_mcp_str=sprintf('%s%s%05d.trc', run.path, opts_MCP.ch_name,(ff+1));
                
            end
        end
    end
    
    filesExists = 0;
    
    if exist(ch_mm_str,'file')==2 && exist(ch_mcp_str,'file')==2 && exist(ch_tr_str,'file')==2
        filesExists=1;
        
        % decode binary files
        ch_mm = ReadLeCroyBinaryWaveform(ch_mm_str);
        ch_mcp = ReadLeCroyBinaryWaveform(ch_mcp_str);
        ch_tr = ReadLeCroyBinaryWaveform(ch_tr_str);
        
        str_disp=sprintf('Processing file set No. %d', ff);
        disp(str_disp);
    else
        str_disp=sprintf('Skipping file set No. %d', ff);
        disp(str_disp);
    end
    
    if filesExists==1
        % get number of segments (events) in one file
        nTRCseg = ch_mcp.info.nbSegments;
        % get segment length
        lTRCseg = size(ch_mcp.y,1);
        % calculate sampling time
        Ts = ch_mcp.x(2,1) - ch_mcp.x(1,1);
        
        % go trough all of the events in the file
        for m=1:nTRCseg
            % generate virtual time vector (important to take care for trigger
            % offset in LeCroy scope)
            t_vec_mm=(0:lTRCseg-1)'*Ts + ch_mm.trigger_offset(m);
            t_vec_mcp=(0:lTRCseg-1)'*Ts + ch_mcp.trigger_offset(m);
            
            % subtract the earliest time
            etime=min([t_vec_mm; t_vec_mcp]);
            
            dutSignalValid = false;
            refSignalValid = false;
            trackerSignalValid = false;
            
            MM_temp.x = 0;
            MM_temp.y = 0;
            
            %% check if DUT, REF and tracker contain valid signals
            if run.savedSignals<saveSignalWaveformsNumber
            shouldSave = 1;
            else 
            shouldSave = 0;
            end
            if(minuit==1)
                MM_temp = process_signal_minuit_picolarge(t_vec_mm-etime,ch_mm.y(:,m),opts_MM,run,shouldSave,storeFolderSignalDUT);
            else
                MM_temp = process_signal(t_vec_mm-etime,ch_mm.y(:,m),opts_MM);
            end
            if(MM_temp.fail==0)
                dutSignalValid = true;
            end
            
            if(minuit==1)
                MCP_temp = process_signal_minuit(t_vec_mcp-etime,ch_mcp.y(:,m),opts_MCP,run,shouldSave,storeFolderSignalRef);
            else
                MCP_temp = process_signal(t_vec_mcp-etime,ch_mcp.y(:,m),opts_MCP);
            end
            
            
            run.savedSignals = run.savedSignals +1;
            
            if(MCP_temp.fail==0)
                refSignalValid = true;
            end
            
                xPos = 0;
                yPos = 0;
                currentEventID = 0;

         if tracker.en
                % process tracker ID channel
                event_id = process_tr_bitstream(t_vec_mcp, ch_tr.y(:,m), opts_TR);
                
                % count overflows
                if(event_id < event_id_prev)
                    event_id_ov = event_id_ov + 1;
                end
                event_id_prev = event_id;
                
                currentEventID = event_id_ov * 65536 + event_id;
                MM_temp.event_id = currentEventID;
                
                MCP_temp.event_id = MM_temp.event_id;
                
                %find corresponding tracker entry and save XY info
                trackIdx = find(tracker.data(:,1) == MM_temp.event_id);
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
                    trackerSignalValid = true;
                end
            end
            eventsValidREF(j) = refSignalValid;
            eventsValidDUT(j) = dutSignalValid;
            eventsValidTracker(j) = trackerSignalValid;
            eventsTrackerX(j) = xPos;
            eventsTrackerY(j) = yPos;
            
            %refSignalValid
            dutSignalValid;
            if dutSignalValid
            else
                MM_temp.fail;
            end
            %trackerSignalValid
            
%                     subplot(1,2,1);
%                     plot(t_vec_mm-etime,ch_mm.y(:,m));
%                     hold on;
%                     plot(t_vec_mcp-etime,ch_mcp.y(:,m));
%                     grid on
%                     xlim([200e-9 300e-9])
%                     title('DUT / REF');
%                     ylabel('Voltage, V');
%                     xlabel('Time, ns');
%                     legend('DUT','REF');
%                     grid on;
%                                         
%                     subplot(1,2,2);
%                     plot(t_vec_mcp-etime,ch_tr.y(:,m));
%                     grid on
%                     title('Event counter');
%                     ylabel('Voltage, V');
%                     xlabel('Time, ns');
%                     legend('DUT','REF');
%                     grid on;
            
            %pause(1);
           % close all
            
            if refSignalValid & dutSignalValid & trackerSignalValid
                eventsTimeDiffSigmoid(j) = MM_temp.sigmoid.timepoint-MCP_temp.sigmoid.timepoint;
                eventsREFmaxY(j) = MCP_temp.sig.max.y;
                eventsDUTmaxY(j) = MM_temp.sig.max.y;
            else 
                eventsTimeDiffSigmoid(j) = 0;
                eventsREFmaxY(j) = 0;
                eventsDUTmaxY(j) = 0;
            end
            
            %% checked if signals valid
            
            %if valid, save
            if refSignalValid & dutSignalValid & trackerSignalValid
                MM_temp.x = xPos;
                MM_temp.y = yPos;
                trackerX(k) = MM_temp.x;
                trackerY(k) = MM_temp.y;
                if tracker.en == 1
                eventIDArray(k) = MM_temp.event_id;
                else 
                    eventIDArray(k) = 0;
                end
                
                if enDebugPlots == true & debugPlotCounter<numberDebugPlots
                    subplot(1,2,1);
                    plot(t_vec_mm-etime,ch_mm.y(:,m));
                    hold on;
                    plot(t_vec_mcp-etime,ch_mcp.y(:,m));
                    grid on
                    xlim([200e-9 300e-9])
                    title('DUT / REF');
                    ylabel('Voltage, V');
                    xlabel('Time, ns');
                    legend('DUT','REF');
                    grid on;
                                        
                    subplot(1,2,2);
                    plot(t_vec_mcp-etime,ch_tr.y(:,m));
                    grid on
                    title('Event counter');
                    ylabel('Voltage, V');
                    xlabel('Time, ns');
                    legend('DUT','REF');
                    grid on;
                    
                    %saveDebugPlot
                    saveas(gcf,[run.debugPlotsPath '\Run' run.id '_' int2str(debugPlotCounter) '_signals.png'])
                    pause(1);
                    close all;
                    debugPlotCounter = debugPlotCounter+1;
                    
                 
                end
                
                % store valid data into array of structures
                MM_data(k)= MM_temp;
                MCP_data(k)= MCP_temp;
               % time_diff(k) = MM_data(k).cfd.time-MCP_data(k).cfd.time;
                time_diff(k) = MM_data(k).sigmoid.timepoint-MCP_data(k).sigmoid.timepoint; %cfd often fails, not used, overwrite
                time_diff_sigmoid(k) = MM_data(k).sigmoid.timepoint-MCP_data(k).sigmoid.timepoint;
                MCP_maxy(k) = MCP_data(k).sig.max.y;
                MM_maxy(k) = MM_data(k).sig.max.y;
                trackerX(k) = MM_temp.x;    
                trackerY(k) = MM_temp.y;
                eventIDArray(k) = MM_temp.event_id;
                padIDArray(k) = run.pad;
                acceptedSignalArray(k) = 1;
                
                k=k+1;
            else
                %not accepted fit, still save waveform

    
    xMM = ch_mm.x(:,m);
    yMM = ch_mm.y(:,m);
        if(opts_MM.invert)
        yMM =-yMM;        % invert signal from MM 
    end
                MM_temp.sig.blrms = std(yMM(10:100));
                MM_temp.sig.blavg = mean(yMM(10:100));

                yMM = yMM - MM_temp.sig.blavg;
                MM_temp.sig.length = length(yMM);
                [MM_temp.sig.max.y, MM_temp.sig.max.idx] = max(yMM);


    xRef = ch_mcp.x(:,m);
    yRef = ch_mcp.y(:,m);
        if(opts_MCP.invert)
        yRef =-yRef;        % invert signal from MM 
    end
                MCP_temp.sig.blrms = std(yRef(10:100));
                MCP_temp.sig.blavg = mean(yRef(10:100));

                yRef = yRef - MCP_temp.sig.blavg;
                MCP_temp.sig.length = length(yRef);
                [MCP_temp.sig.max.y, MCP_temp.sig.max.idx] = max(yRef);


numberSamplesExtent = 20;
    signal.start = MM_temp.sig.max.idx - numberSamplesExtent;  % 5 samples in original code
    if (signal.start < 1)
        signal.start = 1;
    end
    signal.end = MM_temp.sig.max.idx + numberSamplesExtent;
    if (signal.end > MM_temp.sig.length)   
        signal.end = MM_temp.sig.length;
    end
     xfit = xMM(signal.start:signal.end);
     yfit = yMM(signal.start:signal.end);


     numberSamplesExtent = 20;
    signal.start = MCP_temp.sig.max.idx - numberSamplesExtent;  % 5 samples in original code
    if (signal.start < 1)
        signal.start = 1;
    end
    signal.end = MCP_temp.sig.max.idx + numberSamplesExtent;
    if (signal.end > MCP_temp.sig.length)   
        signal.end = MCP_temp.sig.length;
    end
     xfitMCP = xRef(signal.start:signal.end);
     yfitMCP = yRef(signal.start:signal.end);



                MM_temp.cfd = 0;
                MM_temp.sig.blrms = std(yMM(10:100));
                MM_temp.sig.blavg = mean(yMM(10:100));
                MM_temp.sig.xfit = xfit;
                MM_temp.sig.yfit = yfit;
                MM_temp.sigmoid = 0;
                MM_temp.x = xPos;
                MM_temp.y = yPos;

                MCP_temp.cfd = 0;
                MCP_temp.sig.blrms = std(ch_mcp.y(10:100,m));
                MCP_temp.sig.blavg = mean(ch_mcp.y(10:100,m));
                MCP_temp.sig.xfit = xfitMCP;
                MCP_temp.sig.yfit = yfitMCP;
                MCP_temp.sigmoid = 0;

                MM_temp.event_id = currentEventID;
                MM_data(k)= MM_temp;
                MCP_data(k)= MCP_temp;
                time_diff(k) = 0;
                time_diff_sigmoid(k) = 0;
                MCP_maxy(k) = MCP_data(k).sig.max.y;
                MM_maxy(k) = MM_data(k).sig.max.y;
                trackerX(k) = MM_temp.x;    
                trackerY(k) = MM_temp.y;
                eventIDArray(k) = MM_temp.event_id;
                padIDArray(k) = run.pad;
                acceptedSignalArray(k) = 0;
                
                k=k+1;

            
            end
            
            j = j+1;
        end
        toc
    end


    if k>numberEventsPerBatch || ff==run.nfiles
    
    %finished batch of read and process files, then go to next
                str_disp=sprintf('Processed batch, total events analysed: %i latestEventID: %i', k,eventIDArray(k-1));
                disp(str_disp);
                
                str_disp=sprintf('Saving to Matfile');
                disp(str_disp);
                
                storeMatfilePath = [storeMatfileFolder '\processedBatch_' int2str(processingBatchCounter) '_pad' int2str(run.pad) '.mat'];
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
                save(storeMatfilePath,'padIDArray','-append');
                save(storeMatfilePath,'acceptedSignalArray','-append');
                
                
                clear MM_data MCP_data time_diff time_diff_sigmoid MCP_maxy MM_maxy trackerX trackerY eventIDArray padIDArray acceptedSignalArray;
                k = 1;
                processingBatchCounter = processingBatchCounter+1;
                str_disp=sprintf('Saved variables');
                disp(str_disp);

    end
end

toc