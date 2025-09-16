% MATLAB/OCTAVE implementation of processing code for PICOSEC
% version date: 7 Nov 2021

% clear everything
clear all
close all
% run file coordinates and number of files
minuit = 1;   % set minuit for optimizer if 1 or matlab 0
run.id = '558';
run.year = '2022/07';
run.name = ['PICOSEC - single PE ' run.year ' - RUN ' run.id];
run.path=['C:\Users\GDD\Documents\Marta\TestBeamJuly22\Run' run.id '\'];
%run.path=['M:\CommonProject\PICOSEC\Testbeams\July2021\Oscilloscopes\POOL2\Run' run.id '\']; % use mounted Cernbox folder
% run.path=['D:\Matlab_timing\Run' run.id '\']; % use local folder with scope data
%run.path=['C:\Users\PICOSEC\Documents\MATLAB\Analysis\Timingv4 - Laser\LaserTest\Run' run.id '\'];

% run.path=['C:\Users\Marinko\Desktop\PICOSEC_tim\runs\beam2021\Run' run.id '\'];
run.nfiles = find_fileNo(run.path);
%run.nfiles = 9; % override number of file sets to process (see how many trc files are in the folder)
run.lecroy_name = '--Trace--'; %['Run' run.id];

% options for processing micromegas channel (lot magic numbers from
% processing functionhave to be added here)
opts_MM.t_dead = 1;      % blanking time from start of the recording samples
opts_MM.t_prerms = 20;  % blanking time before global maximum samples
opts_MM.Ts = 1/10e9;     % sampling speed
opts_MM.Rin=50;          % input impedance
opts_MM.invert = 0;      % is inverted (1-inverted, 0-non-inverted)
opts_MM.type=0;          % detector type (0-DUT, 1-REF) for future use
opts_MM.ch_name = ['C4' run.lecroy_name]; % LeCroy file name format
opts_MM.en_plot = 1;     % enable debugging plots

% options for processing MCP channel (lot magic numbers from processing
% functionhave to be added here)
opts_MCP.t_dead = 1;
opts_MCP.t_prerms = 20;
opts_MCP.Ts = 1/10e9;
opts_MCP.Rin=50;
opts_MCP.invert = 0;
opts_MCP.type=1;
opts_MCP.ch_name = ['C4' run.lecroy_name];
opts_MCP.en_plot = 0;


%% DO PROCCEISNG
tic;         % start time measurement
k=1;         % valid data counter

maxAmpMM = [];
chiVector = [];

for (ff=0:run.nfiles)
    % display current file number
    str_disp=sprintf('Loading file set No. %d', ff);
    disp(str_disp);
    
    % generate file name strings
    ch_mm_str=sprintf('%s%s%05d.trc', run.path, opts_MM.ch_name,ff);
    ch_mcp_str=sprintf('%s%s%05d.trc', run.path, opts_MCP.ch_name,ff);
    
    % read files using third party function
    ch_mm = ReadLeCroyBinaryWaveform(ch_mm_str);
    ch_mcp = ReadLeCroyBinaryWaveform(ch_mcp_str);
    
    str_disp=sprintf('Processing file set No. %d', ff);
    disp(str_disp);
    
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
        
        % process MM Picosec first to see if signal is valid
        if(minuit==1)
            MM_temp = process_signal_minuit(t_vec_mm-etime,ch_mm.y(:,m),opts_MM);
            %pause(1);
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
               % MM_temp.event_id = process_tr_bitstream(t_vec_mcp, ch_tr.y(:,m), opts_TR);
               % MCP_temp.event_id =MM_temp.event_id;
                % store valid data into array of structures
                MM_data(k)= MM_temp;
                MCP_data(k)= MCP_temp;
                time_diff(k) = MM_data(k).cfd.time-MCP_data(k).cfd.time;
                time_diff_sigmoid(k) = MM_data(k).sigmoid.timepoint-MCP_data(k).sigmoid.timepoint;
                MCP_maxy(k) = MCP_data(k).sig.max.y;
                MM_maxy(k) = MM_data(k).sig.max.y;
                
                k=k+1;
            end
       end
    end
    toc
end

%% save data to mat file
save(['Run' run.id '.mat'], 'run', 'MM_data', 'MCP_data', 'time_diff', 'time_diff_sigmoid', 'MCP_maxy', 'MM_maxy');



toc