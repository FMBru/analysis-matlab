% MATLAB/OCTAVE implementation of processing code for PICOSEC
% version date: 7 Nov 2021

% clear everything
clear all
close all
% run file coordinates and number of files
minuit = 1;   % set minuit for optimizer if 1 or matlab 0
run.id = '103';
run.year = '2021/10 MCP';
run.name = ['BEAM ' run.year ' RUN ' run.id];
%run.path=['M:\CommonProject\PICOSEC\Testbeams\July2021\Oscilloscopes\POOL2\Run' run.id '\']; % use mounted Cernbox folder
% run.path=['D:\Matlab_timing\Run' run.id '\']; % use local folder with scope data
run.path=['E:\MCPtiming\Run103GDD\Run' run.id '\'];

% run.path=['C:\Users\Marinko\Desktop\PICOSEC_tim\runs\beam2021\Run' run.id '\'];
run.nfiles = find_fileNo(run.path);
%run.nfiles = 50; % override number of file sets to process (see how many trc files are in the folder)
run.lecroy_name = '--Trace--'; %['Run' run.id];


%read tracker data from ASCII file
tracker.path = ['E:\MCPtiming\Run103GDD\run0103ascii.dat'];
tracker.en = 1; %match eventIDs to tracking data and add XY to output
tracker.dutIndex = 1; %1: DUT1 (Multipad), 2: DUT2 (MultipadRes), 3: DUT3 (ThinGap), 4:DUT4 (VacChamber), 5:DUT0(Cactus)
dutColArray = [1 7 8; 2 10 11; 3 13 14; 4 16 17; 5 4 5];

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
opts_MM.ch_name = ['C2' run.lecroy_name]; % LeCroy file name format
opts_MM.en_plot = 0;     % enable debugging plots

% options for processing MCP channel (lot magic numbers from processing
% functionhave to be added here)
opts_MCP.t_dead = 1;
opts_MCP.t_prerms = 200;
opts_MCP.Ts = 1/20e9;
opts_MCP.Rin=50;
opts_MCP.invert = 1;
opts_MCP.type=1;
opts_MCP.ch_name = ['C1' run.lecroy_name];
opts_MCP.en_plot = 1;

% options for extracting tracker ID from the bitstream
opts_TR.ch_name = ['C3' run.lecroy_name];
opts_TR.baud_rate = 40e6; % 40Mbps baudrate
opts_TR.n_bits = 16;      % number of bits after start bit

%% DO PROCCEISNG
tic;         % start time measurement
k=1;         % valid data counter

for (ff=0:run.nfiles)
    % display current file number
    str_disp=sprintf('Loading file set No. %d', ff);
    disp(str_disp);
    
    % generate file name strings
    ch_mm_str=sprintf('%s%s%05d.trc', run.path, opts_MM.ch_name,ff);
    ch_mcp_str=sprintf('%s%s%05d.trc', run.path, opts_MCP.ch_name,ff);
    ch_tr_str=sprintf('%s%s%05d.trc', run.path, opts_TR.ch_name,ff);
    
    % read files using third party function
    ch_mm = ReadLeCroyBinaryWaveform(ch_mm_str);
    ch_mcp = ReadLeCroyBinaryWaveform(ch_mcp_str);
    ch_tr = ReadLeCroyBinaryWaveform(ch_tr_str);
    
    str_disp=sprintf('Processing file set No. %d', ff);
    disp(str_disp);
    
    % get number of segments (events) in one file
    nTRCseg = ch_mcp.info.nbSegments;
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
                MM_temp.event_id = process_tr_bitstream(t_vec_mcp, ch_tr.y(:,m), opts_TR);
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
    toc
end

%% save data to mat file
save(['Run' run.id '.mat'], 'run', 'MM_data', 'MCP_data', 'time_diff', 'time_diff_sigmoid', 'MCP_maxy', 'MM_maxy', 'trackerX', 'trackerY');



toc