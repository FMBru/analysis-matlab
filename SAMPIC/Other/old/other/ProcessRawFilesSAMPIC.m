% MATLAB/OCTAVE implementation of processing code for PICOSEC
% version date: 7 Nov 2021

close all

addpath 'C:\Users\GDD\Documents\MATLAB\Picosec\CommonFunctions';
addpath 'C:\Users\GDD\Documents\MATLAB\Picosec\CommonFunctions\SCP';

processingFolder = ['\\eosproject-smb\eos\project\p\picosec\testbeam\2022_May_h4\SAMPIC\Processing\Run' run.id]; 
signalsFolder = [processingFolder '\Signals']; 
infoFile = [processingFolder '\info.mat'];
load(infoFile);

matchedEventsFolder = [processingFolder '\MatchedEvents'];


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
    opts_MM.chID = 4;
    opts_MCP.chID = 0;
    opts_MM.ch_name = ['C4' ]; % LeCroy file name format MCP1 (used for timing)
    opts_MCP.ch_name = ['C2' ]; % MPC 2 (used as trigger)
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



for segmentCounter = 1:numberSegments
    %eval(['matchEventsDUT = matchedEvents_ch' num2str(opts_MM.chID) ';']);
    
        matchedEventsSegmentPath = [matchedEventsFolder '\ch' num2str(opts_MM.chID) '_segment' num2str(segmentCounter) '.mat'];
        load(matchedEventsSegmentPath);
        
    matchEventsDUT = matchedEventsTemp;

%% DO PROCCEISNG
tic;         % start time measurement
k=1;         % valid data counter


        % go trough all of the events in matched events vector
        for m=1:length(matchEventsDUT)
            
            if mod(m,100)==0
                str_disp=sprintf('Processing event No. %d', m);
         disp(str_disp);
            end

            %structure of device hit
            %cell0Time
            %unixTimestamp
            %ch
            %hit
            %waveform

            dut = matchEventsDUT(m).dut;
            dutWaveform = dut.waveform;
            
            ref = matchEventsDUT(m).ref;
            refWaveform = ref.waveform;
            
            event_id = matchEventsDUT(m).eventID;
            
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
    
    toc
end
    %% save data to mat file
if shouldSaveMAT
    save(['C:\Users\GDD\Documents\Picosec\May22\Analysed\Run' run.id '-' run.oscilloscope '.mat'], 'run', 'MM_data', 'MCP_data', 'time_diff', 'time_diff_sigmoid', 'MCP_maxy', 'MM_maxy', 'trackerX', 'trackerY');
end
toc