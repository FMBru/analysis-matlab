% MATLAB/OCTAVE implementation of processing code for PICOSEC
% version date: 7 Nov 2021

close all

lastEvents = true;
middleEvents = false;
eventCutter = 200;

% clear everything if not batch processing
%set run parameters
if exist('batchProcess','var') == 1
    run.id = runIDString;
    run.oscilloscope = oscilloscopeString;
    run.lecroy_name = '--Trace--'; %['Run' run.id];
    opts_MM.ch_name = [DUTnameString run.lecroy_name]
    opts_MCP.ch_name = [MCPnameString run.lecroy_name]
    shouldSaveMAT = true;
    shouldUseEOSFolder = useEOSFolder;
    tracker.dutIndex = str2num(trackerDUTIndex); %1: miniCactus, 2: MM3 (VacChamberPos), 3: MM1(Multipad), 4: MM2 (on support plate), 5: MM4 (electron setup)
    
else
    %not batch processing
    %clear all
    run.lecroy_name = '--Trace--'; %['Run' run.id];
    
    run.id = '255';
    run.oscilloscope = 'Pool5';
    opts_MM.ch_name = ['C4' run.lecroy_name]; % LeCroy file name format MCP1 (used for timing)
    opts_MCP.ch_name = ['C2' run.lecroy_name]; % MPC 2 (used as trigger)
    shouldSaveMAT = false;
    shouldUseEOSFolder = true;
    tracker.dutIndex = 2; %1: miniCactus, 2: MM3 (VacChamberPos), 3: MM1(Multipad), 4: MM2 (on support plate), 5: MM4 (electron setup)
    
end

%correct for mistake in GDD scope where C2 files are saved starting at
%00001 instead of 00000
shouldCorrectFileCounter_GDD_C2 = true;


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
    run.nfiles = find_fileNo(run.pathEOS)
    tracker.path = ['C:\Users\GDD\Documents\Picosec\May22\Tracker\reconstructed\asciiRun' run.id '.dat'];
    
    %tracker.path = ['\\eosproject-smb\eos\project\p\picosec\testbeam\2022_May_h4\tracker\reconstructed\asciiRun' run.id '.dat'];
else
    run.nfiles = find_fileNo(run.path);
    tracker.path = ['C:\Users\GDD\Documents\Picosec\May22\Tracker\reconstructed\asciiRun' run.id '.dat'];
end
%run.nfiles = 50; % override number of file sets to process (see how many trc files are in the folder)
%run.lecroy_name = '--Trace--'; %['Run' run.id];
run.lecroy_name = '--Trace--'; %['Run' run.id];


%read tracker data from ASCII file
%tracker.path = ['E:\MCPtiming\Run103GDD\run0103ascii.dat'];
tracker.en = 1; %match eventIDs to tracking data and add XY to output
%MCP run

%dutColArray = [1 4 5; 2 7 8; 3 10 11; 4 13 14; 5 16 17]; %[dID colXID colYID; ] -> from MCP run
%dutColArray = [1 4 5; 2 7 8; 3 10 11; 4 13 14; 5 16 17]; %[dID colXID colYID; ] -> from MCP run
dutColArray = [1 4 5; 2 7 8; 3 10 11; 4 13 14; 5 16 17]; %[dID colXID colYID; ] -> May MM run

if tracker.en
    trackerFile = fopen(tracker.path,'rt');
    D = textscan(trackerFile, '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f', 'Delimiter','\t', 'HeaderLines',2, 'CollectOutput',1);
    tracker.data = cell2mat(D);
end

% options for processing micromegas channel (lot magic numbers from
% processing functionhave to be added here)
opts_MM.t_dead = 1;      % blanking time from start of the recording samples
opts_MM.t_prerms = 200;  % blanking time before global maximum samples
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
opts_MCP.t_prerms = 200;
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
event_id_prev = -1;
event_id_ov = 0;

for (ff=0:run.nfiles)
    % display current file number
    str_disp=sprintf('Loading file set No. %d', ff);
    disp(str_disp);
    
    % generate file name strings
    
    
    ch_mm_str=sprintf('%s%s%05d.trc', run.path, opts_MM.ch_name,ff);
    ch_mcp_str=sprintf('%s%s%05d.trc', run.path, opts_MCP.ch_name,ff);
    ch_tr_str=sprintf('%s%s%05d.trc', run.path, opts_TR.ch_name,ff);
    
    
    if shouldUseEOSFolder
        ch_mm_str=sprintf('%s%s%05d.trc', run.pathEOS, opts_MM.ch_name,ff);
        ch_mcp_str=sprintf('%s%s%05d.trc', run.pathEOS, opts_MCP.ch_name,ff);
        ch_tr_str=sprintf('%s%s%05d.trc', run.pathEOS, opts_TR.ch_name,ff);
    end
    
    if shouldCorrectFileCounter_GDD_C2
        if strcmp(run.oscilloscope,'GDD')
            %shift counter for C2 by one
            if strcmp(opts_MM.ch_name,['C2' run.lecroy_name])
                ch_mm_str=sprintf('%s%s%05d.trc', run.path, opts_MM.ch_name,(ff+1));
                
            end
            if strcmp(opts_MCP.ch_name,['C2' run.lecroy_name])
                ch_mcp_str=sprintf('%s%s%05d.trc', run.path, opts_MCP.ch_name,(ff+1));
                
            end
        end
    end
    
    filesExists = 0;
    
    if exist(ch_mm_str,'file')==2 && exist(ch_mcp_str,'file')==2 && exist(ch_tr_str,'file')==2
        filesExists=1;
        
        
        % read files using third party function
        ch_mm = ReadLeCroyBinaryWaveform(ch_mm_str);
        ch_mcp = ReadLeCroyBinaryWaveform(ch_mcp_str);
        ch_tr = ReadLeCroyBinaryWaveform(ch_tr_str);
        
%                 figure
%                 histogram(ch_mm.trigger_time,100)
%                 title("Time distribution before cut");
        
        if(lastEvents)
            ch_mm.trigger_time = ch_mm.trigger_time(end-eventCutter:end);
            ch_mm.trigger_offset = ch_mm.trigger_offset(end-eventCutter:end);
            ch_mm.x = ch_mm.x(:,end-eventCutter:end);
            ch_mm.y = ch_mm.y(:,end-eventCutter:end);
            
            ch_mcp.trigger_time = ch_mcp.trigger_time(end-eventCutter:end);
            ch_mcp.trigger_offset = ch_mcp.trigger_offset(end-eventCutter:end);
            ch_mcp.x = ch_mcp.x(:,end-eventCutter:end);
            ch_mcp.y = ch_mcp.y(:,end-eventCutter:end);
            
            ch_tr.trigger_time = ch_tr.trigger_time(end-eventCutter:end);
            ch_tr.trigger_offset = ch_tr.trigger_offset(end-eventCutter:end);
            ch_tr.x = ch_tr.x(:,end-eventCutter:end);
            ch_tr.y = ch_tr.y(:,1:end- eventCutter);
            
        elseif(middleEvents)
            timeDiff = ch_mm.trigger_time(2:size(ch_mm.trigger_time,1)) - ...
                ch_mm.trigger_time(1:size(ch_mm.trigger_time,1)-1);
            
            timeDiffThr = 0.2;
            
            densityEvents = timeDiff < timeDiffThr;
            
            indexNoise = find(densityEvents == 0);
            
            if(isempty(indexNoise))
                indexNoise = 0;
            end
            
            indexMiddle = (length( ch_mm.trigger_time) + indexNoise(end)+1)/2;
            indexBegin = indexMiddle - eventCutter/2;
            indexEnd = indexMiddle + eventCutter/2;
            
            if(~isempty(indexNoise) && (indexBegin <= indexNoise(end) || indexEnd > length(ch_mm.trigger_time)))
                disp("cutted")
                disp(string(indexNoise(end)))
                continue;
            end
            
            ch_mm.trigger_time = ch_mm.trigger_time(indexBegin:indexEnd);
            ch_mm.trigger_offset = ch_mm.trigger_offset(indexBegin:indexEnd);
            ch_mm.x = ch_mm.x(:,indexBegin:indexEnd);
            ch_mm.y = ch_mm.y(:,indexBegin:indexEnd);
            
            ch_mcp.trigger_time = ch_mcp.trigger_time(indexBegin:indexEnd);
            ch_mcp.trigger_offset = ch_mcp.trigger_offset(indexBegin:indexEnd);
            ch_mcp.x = ch_mcp.x(:,indexBegin:indexEnd);
            ch_mcp.y = ch_mcp.y(:,indexBegin:indexEnd);
            
            ch_tr.trigger_time = ch_tr.trigger_time(indexBegin:indexEnd);
            ch_tr.trigger_offset = ch_tr.trigger_offset(indexBegin:indexEnd);
            ch_tr.x = ch_tr.x(:,indexBegin:indexEnd);
            ch_tr.y = ch_tr.y(:,indexBegin:indexEnd);
            
        else
            
            timeDiff = ch_mm.trigger_time(2:size(ch_mm.trigger_time,1)) - ...
                ch_mm.trigger_time(1:size(ch_mm.trigger_time,1)-1);
            
            %             figure
            %             histogram(timeDiff,100)
            %             title("Time difference in file");
            
            timeDiffThr = 0.2;
            
            densityEvents = timeDiff < timeDiffThr;
            
            indexNoise = find(densityEvents == 0);
            
            if(~isempty(indexNoise) && eventCutter+indexNoise(end)+1>length(ch_mm.trigger_time))
                disp("cutted")
                disp(string(indexNoise(end)))
                continue;
            end
            
            if(isempty(indexNoise))
                indexNoise = 0;
            end
            
            ch_mm.trigger_time = ch_mm.trigger_time(indexNoise(end)+1:eventCutter+indexNoise(end)+1);
            ch_mm.trigger_offset = ch_mm.trigger_offset(indexNoise(end)+1:eventCutter+indexNoise(end)+1);
            ch_mm.x = ch_mm.x(:,indexNoise(end)+1:eventCutter+indexNoise(end)+1);
            ch_mm.y = ch_mm.y(:,indexNoise(end)+1:eventCutter+indexNoise(end)+1);
            
            ch_mcp.trigger_time = ch_mcp.trigger_time(indexNoise(end)+1:eventCutter+indexNoise(end)+1);
            ch_mcp.trigger_offset = ch_mcp.trigger_offset(indexNoise(end)+1:eventCutter+indexNoise(end)+1);
            ch_mcp.x = ch_mcp.x(:,indexNoise(end)+1:eventCutter+indexNoise(end)+1);
            ch_mcp.y = ch_mcp.y(:,indexNoise(end)+1:eventCutter+indexNoise(end)+1);
            
            ch_tr.trigger_time = ch_tr.trigger_time(indexNoise(end)+1:eventCutter+indexNoise(end)+1);
            ch_tr.trigger_offset = ch_tr.trigger_offset(indexNoise(end)+1:eventCutter+indexNoise(end)+1);
            ch_tr.x = ch_tr.x(:,indexNoise(end)+1:eventCutter+indexNoise(end)+1);
            ch_tr.y = ch_tr.y(:,indexNoise(end)+1:eventCutter+indexNoise(end)+1);
            
        end
%         
%                 figure
%                 histogram(ch_mm.trigger_time,100)
%                 title("Time distribution after cut");
        
        str_disp=sprintf('Processing file set No. %d', ff);
        disp(str_disp);
        
    else
        str_disp=sprintf('Skipping file set No. %d', ff);
        disp(str_disp);
        
    end
    
    if filesExists==1
        
        % get number of segments (events) in one file
        nTRCseg = length(ch_mcp.trigger_time);
        % get segment length
        lTRCseg = length(ch_mcp.y);
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
            
            % process MM Picosec first to see if signal is valid
            if(minuit==1)
                MM_temp = process_signal_minuit(t_vec_mm-etime,ch_mm.y(:,m),opts_MM);
            else
                MM_temp = process_signal(t_vec_mm-etime,ch_mm.y(:,m),opts_MM);
            end
            
            % if signal is valid process MCP
            if(MM_temp.fail==0)
                if(minuit==1)
                    MCP_temp = process_signal_minuit(t_vec_mcp-etime,ch_mcp.y(:,m),opts_MCP);
                else
                    MCP_temp = process_signal(t_vec_mcp-etime,ch_mcp.y(:,m),opts_MCP);
                end
                
                % if MCP is valid store data to structure array
                if(MCP_temp.fail==0)
                    % process tracker ID channel
                    event_id = process_tr_bitstream(t_vec_mcp, ch_tr.y(:,m), opts_TR);
                    
                    % count overflows
                    if(event_id < event_id_prev)
                        event_id_ov = event_id_ov + 1;
                    end
                    event_id_prev = event_id;
                    
                    MM_temp.event_id = event_id_ov * 65536 + event_id;
                    MCP_temp.event_id = MM_temp.event_id;
                    
                    
                    xPos = 0;
                    yPos = 0;
                    if tracker.en
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
                        end
                        
                    end
                    MM_temp.x = xPos;
                    MM_temp.y = yPos;
                    trackerX(k) = MM_temp.x;
                    trackerY(k) = MM_temp.y;
                    eventIDArray(k) = MM_temp.event_id;
                    
                    
                    % store valid data into array of structures
                    MM_data(k)= MM_temp;
                    MCP_data(k)= MCP_temp;
                    time_diff(k) = MM_data(k).cfd.time-MCP_data(k).cfd.time;
                    time_diff_sigmoid(k) = MM_data(k).sigmoid.timepoint-MCP_data(k).sigmoid.timepoint;
                    MCP_maxy(k) = MCP_data(k).sig.max.y;
                    MM_maxy(k) = MM_data(k).sig.max.y;
                    trackerX(k) = MM_temp.x;
                    trackerY(k) = MM_temp.y;
                    
                    k=k+1;
                end
            end
        end
    end
    toc
end
%% save scope settings
run.scope_set_ch_mcp = ch_mcp.info;
run.scope_set_ch_mm = ch_mm.info;
run.scope_set_ch_tr = ch_tr.info;

%s1Signal = MM_maxy;
%s2Signal = MCP_maxy;
%eventIDScintillators = eventIDArray;
%save(['C:\Users\GDD\Documents\Picosec\May22\Analysed\Run' run.id '-' run.oscilloscope '-S1+S2.mat'], 'run', 's1Signal', 's2Signal', 'eventIDScintillators');


%% save data to mat file
if shouldSaveMAT
    save(['C:\Users\GDD\Documents\Picosec\May22\Analysed\Marta\Run' run.id '-' run.oscilloscope '.mat'], 'run', 'MM_data', 'MCP_data', 'time_diff', 'time_diff_sigmoid', 'MCP_maxy', 'MM_maxy', 'trackerX', 'trackerY');
end
toc