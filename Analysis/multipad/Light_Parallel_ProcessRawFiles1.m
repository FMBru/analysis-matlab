% MATLAB/OCTAVE implementation of processing code for PICOSEC
% version date: 7 Nov 2021
close all
addpath '\\eosproject-smb\eos\project\p\picosec\testbeam\2025_July_h4\Matlab\CommonFunctions';
addpath '..\..\Analysis\Functions';

%May22: %1: miniCactus, 2: MM3 (VacChamberPos), 3: MM1(Multipad), 4: MM2 (on support plate), 5: MM4 (electron setup)

% clear everything if not batch processing
%set run parameters
if exist('batchProcess','var') == 1
    run.id = runIDString;
    run.oscilloscope = oscilloscopeString;
    run.lecroy_name = '--Trace--'; %['Run' run.id];
    opts_MM.chID = DUTnameString;
    opts_MCP.chID = MCPnameString;
    opts_MM.ch_name = ['C' DUTnameString run.lecroy_name];
    opts_MCP.ch_name = ['C' MCPnameString run.lecroy_name];
    shouldSaveMAT = true;
    shouldUseEOSFolder = useEOSFolder;
    tracker.dutIndex = str2num(trackerDUTIndex);  %DUT index: 1: MM1 (Multipad), 2: MM2(VacChamber), 3: MM3 (ElectronsetupMM), 4: MM4 (Picolarge main tracker), 5:MM5 (SiPM)
    numberFilesToAnalyse = noFileToAnalyze;
    runInfoString = runInfoDesc;
else
    %not batch processing - List here Channels and DUT if not using batch
    %processing
    clear all

    run.lecroy_name = '--Trace--'; %['Run' run.id];
    
    run.id = '366';
    run.oscilloscope = 'Pool3';
    opts_MM.chID = '4';
    opts_MCP.chID = '1';
    opts_MM.ch_name = ['C' opts_MM.chID run.lecroy_name]; % LeCroy file name format MCP1 (used for timing)
    opts_MCP.ch_name = ['C' opts_MCP.chID run.lecroy_name]; % MPC 2 (used as trigger)
    shouldSaveMAT = false;
    shouldUseEOSFolder = true;
    tracker.dutIndex = 1; %DUT index: 1: MM1 (Multipad), 2: MM2(VacChamber), 3: MM3 (ElectronsetupMM), 4: MM4 (Picolarge main tracker), 5:MM5 (SiPM)
    numberFilesToAnalyse = 0; %max number of files to analyse, 0 -> analyse all in folder
    runInfoString = 'Multipad double DLC pad 34 - 180 um - PC Ti 2p4 nm - Flushing - A275V C550V (1 file)';
end

enDebugPlots = false;
correctFilenumberOffset = false;

% if strcmp(run.oscilloscope, 'Pool2') && strcmp(opts_MM.chID, '2')
%     correctFilenumberOffset = true;
%     disp('Corrected for file offset');
% end

shouldPlotOutliers = false;
shouldPlotTrackerWaveform = false;

%read tracker data from ASCII file
tracker.en = 1; %match eventIDs to tracking data and add XY to output

numberDebugPlots = 4;
saveSignalWaveformsNumber = 8;

debugPlotCounter = 0;
run.savedSignals = 0;

%correct for mistake in GDD scope where C2 files are saved starting at
%00001 instead of 00000
shouldCorrectFileCounter_GDD_C2 = false;

% run file coordinates and number of files
minuit = 1;   % set minuit for optimizer if 1 or matlab 0

run.year = '2025 July ';
run.name = ['BEAM ' run.year ' RUN ' run.id];
run.pathEOS=['\\eosproject-smb\eos\project\p\picosec\testbeam\2025_July_h4\' run.oscilloscope '\Run' run.id '\'];
%run.pathEOS=['\\eosproject-smb\eos\project\p\picosec\testbeam\2022_July_h4\' run.oscilloscope '\Run' run.id '\'];
disp(run.pathEOS)
run.nfiles = find_fileNo(run.pathEOS);

tracker.path = ['\\eosproject-smb\eos\project\p\picosec\testbeam\2025_July_h4\tracker\reconstructed\asciiRun' run.id '.dat'];
%tracker.path = ['\\eosproject-smb\eos\project\p\picosec\testbeam\2022_July_h4\tracker\reconstructed\asciiRun' run.id '.dat'];

if(numberFilesToAnalyse>0)
    run.nfiles =numberFilesToAnalyse;
end
    
run.lecroy_name = '--Trace--'; %['Run' run.id];

run.debugPlotsPath=['\\eosproject-smb\eos\project\p\picosec\testbeam\2025_July_h4\Results\DebugPlots\Run' run.id '-'  run.oscilloscope '\'];
mkdir(run.debugPlotsPath);

storeFolderSignalRef = ['\\eosproject-smb\eos\project\p\picosec\testbeam\2025_July_h4\Results\Run' run.id '-' run.oscilloscope '-dut' opts_MM.chID '-ref' opts_MCP.chID '-' runInfoString '\signals\REF'];
mkdir(storeFolderSignalRef);
storeFolderSignalDUT = ['\\eosproject-smb\eos\project\p\picosec\testbeam\2025_July_h4\Results\Run' run.id '-' run.oscilloscope '-dut' opts_MM.chID '-ref' opts_MCP.chID '-' runInfoString '\signals\DUT'];
mkdir(storeFolderSignalDUT);


%Detector Mapping September 2024 
%DUT index: 1:MM1, 2:MM2, 3:MM3, 4:MM4, 5:MM5, 6:MM6 - USTC 7:MM7-PIMENT96
dutColArray = [1 4 5; 2 10 11; 3 13 14; 4 16 17; 5 22 23; 6 25 26; 7 19 20]; %[dID colXID colYID; ] 

if tracker.en
    trackerFile = fopen(tracker.path,'rt');
    D = textscan(trackerFile, '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f', 'Delimiter','\t', 'HeaderLines',2, 'CollectOutput',1);
    tracker.data = cell2mat(D);

    if length(tracker.data)<1000
           error('Tracker file asciiRunXXX.dat required but not found or empty. PLease check tracker reconstruction and ROOT->ASCII conversion.')
    end
end

% options for processing micromegas channel (lot magic numbers from
% processing functionhave to be added here)
opts_MM.t_dead = 1;      % blanking time from start of the recording samples
opts_MM.t_prerms = 20;  % blanking time before global maximum samples
opts_MM.t_prerms = 200;  % blanking time before global maximum samples
opts_MM.Ts = 1/10e9;     % sampling speed
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
opts_MCP.t_prerms = 200;
opts_MCP.Ts = 1/10e9;
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
% tic;         % start time measurement
k=1;         % valid data counter
j=1;         % all data counter
event_id_prev = -1;
event_id_ov = 0;
monitorPlotCounter = 0;

% 
if isempty(gcp('nocreate'))
    parpool('local', 6); %activate a parallelization pool
end

tic;         % start time measurement
% Pre-allocation
eventsValidDUT_all = cell(run.nfiles, 1);
eventsValidREF_all = cell(run.nfiles, 1);
eventsValidTracker_all = cell(run.nfiles, 1);

%MM_data_onlyValid = cell(run.nfiles, 1);
MM_maxy_all = cell(run.nfiles, 1);
MM_maxy_onlyValid = cell(run.nfiles, 1);

e_peak_MM_onlyValid = cell(run.nfiles, 1);
e_lead_MM_onlyValid = cell(run.nfiles, 1);
riseTime_onlyValid = cell(run.nfiles, 1);

%MCP_data_onlyValid = cell(run.nfiles, 1);
MCP_maxy_all = cell(run.nfiles, 1);
MCP_maxy_onlyValid = cell(run.nfiles, 1);

e_peak_MCP_onlyValid = cell(run.nfiles, 1);

eventsTrackerX_all = cell(run.nfiles, 1);
eventsTrackerY_all = cell(run.nfiles, 1);
eventsTrackerX_onlyValid = cell(run.nfiles, 1);
eventsTrackerY_onlyValid = cell(run.nfiles, 1);

time_diff_onlyValid = cell(run.nfiles, 1);

time_diff_sigmoid_onlyValid = cell(run.nfiles, 1);
time_diff_sigmoid_all = cell(run.nfiles, 1);

eventIDArray_all = cell(run.nfiles, 1);

for ff=1:run.nfiles
    % display current file number
    str_disp=sprintf('Loading file set No. %d', ff);
    disp(str_disp);
    
    % generate file name strings
    
    if correctFilenumberOffset
        %file numbers are shifted
        ch_mm_str=sprintf('%s%s%05d.trc', run.pathEOS, opts_MM.ch_name,ff+1);
        ch_mcp_str=sprintf('%s%s%05d.trc', run.pathEOS, opts_MCP.ch_name,ff); %add ff+1 if file numbers are shifted
        ch_tr_str=sprintf('%s%s%05d.trc', run.pathEOS, opts_TR.ch_name,ff); %add ff+1 if file numbers are shifted
    else
        %no problem with filenumbers
        ch_mm_str=sprintf('%s%s%05d.trc', run.pathEOS, opts_MM.ch_name,ff);
        ch_mcp_str=sprintf('%s%s%05d.trc', run.pathEOS, opts_MCP.ch_name,ff); 
        ch_tr_str=sprintf('%s%s%05d.trc', run.pathEOS, opts_TR.ch_name,ff);
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

        if isfield(ch_mm, 'trigger_offset')==0
            filesExists = 0;
            disp('Something went wrong with the oscilloscope');
        end
    else
        str_disp=sprintf('Skipping file set No. %d', ff);
        disp(str_disp);
    end


%     if ff == 2
%         filesExists = 0;
%     end
    
    if filesExists==1
        % get number of segments (events) in one file
        nTRCseg = ch_mcp.info.nbSegments;
        % get segment length
        lTRCseg = size(ch_mcp.y,1);
        % calculate sampling time
        Ts = ch_mcp.x(2,1) - ch_mcp.x(1,1);

        %temporary storage for each segment
        %localMMData = [];
        local_eventIDArray = zeros(1,nTRCseg);
        %localMCPData = [];
        local_time_diff = zeros(1,nTRCseg);
        local_time_diff_sigmoid = zeros(1,nTRCseg);
        local_e_peak_MM = zeros(1,nTRCseg);
        local_e_lead_MM = zeros(1,nTRCseg);
        local_riseTime = zeros(1,nTRCseg);
        local_e_peak_MCP = zeros(1,nTRCseg);
        localMaxYMCP = zeros(1,nTRCseg);
        localMaxY = zeros(1,nTRCseg);
        localValidDUT = zeros(1,nTRCseg);
        localValidREF = zeros(1,nTRCseg);
        localValidTracker = zeros(1,nTRCseg);
        localTrackerX = zeros(1,nTRCseg);
        localTrackerY = zeros(1,nTRCseg);


	    %-------- BROADCASTING --------
	    MM_triggeroffset = ch_mm.trigger_offset;
        MM_Y = ch_mm.y;

	    MCP_triggeroffset = ch_mcp.trigger_offset;
        MCP_Y = ch_mcp.y;

        enableTracker = tracker.en;
        if enableTracker
            TR_Y = ch_tr.y;
            TR_DATA = tracker.data;
            TR_DUTINDEX = tracker.dutIndex;

            %further broadcasting
            TR_DATA_eventIDColumn = TR_DATA(:,1);
            TR_DATA_xPosColumn = TR_DATA(:,dutColArray(TR_DUTINDEX,2));
            TR_DATA_yPosColumn = TR_DATA(:,dutColArray(TR_DUTINDEX,3));
        end


    
	    %only way that I've found to handle overflow problem (any better idea is welcome)
    
	    t_vec_mcp=(0:lTRCseg-1)'*Ts + MCP_triggeroffset(1);	
	    event_id_init = process_tr_bitstream(t_vec_mcp, ch_tr.y(:,1), opts_TR);
	    t_vec_mcp=(0:lTRCseg-1)'*Ts + MCP_triggeroffset(nTRCseg);
	    event_id_fin = process_tr_bitstream(t_vec_mcp, ch_tr.y(:,nTRCseg), opts_TR);
    
	    %if the final event has a counter lower than the initial one overflow has occured	
	    if event_id_init > event_id_fin 
	        overflowOccurs = true;
	    else
	        overflowOccurs = false;
	    end
        
        %if overflow occured during data saving
        if event_id_init < event_id_prev
	        event_id_ov = event_id_ov + 1;
        end

	    if overflowOccurs
	        event_id_ov = event_id_ov + 1;
	    end
        
        %update event_id_prev for the next file
        event_id_prev = event_id_fin;

        % go trough all of the events in the file
        parfor m=1:nTRCseg
        %for m=1:nTRCseg
            
            %disp(['Processing seg. ' int2str(m)])
            MM_temp = struct();

            % generate virtual time vector (important to take care for trigger
            % offset in LeCroy scope)
            t_vec_mm=(0:lTRCseg-1)'*Ts + MM_triggeroffset(m);
            t_vec_mcp=(0:lTRCseg-1)'*Ts + MCP_triggeroffset(m);
            
            % subtract the earliest time
            etime=min([t_vec_mm; t_vec_mcp]);
            
            dutSignalValid = false;
            refSignalValid = false;
            trackerSignalValid = false;
            
            MM_temp.x = 0;
            MM_temp.y = 0;
            
            eventIDTracker = 0;
            xPos = 0;
            yPos = 0;

       	    if enableTracker

                % process tracker ID channel
                event_id = process_tr_bitstream(t_vec_mcp, TR_Y(:,m), opts_TR);
                
                if shouldPlotTrackerWaveform
                    close all
                    plot(t_vec_mcp, TR_Y(:,m));
                    pause(0.5);
                end
                
                if overflowOccurs
	               if event_id < 65536 && event_id > event_id_init
	                  eventIDTracker = (event_id_ov-1) * 65536 + event_id;
	               else
	                  eventIDTracker = event_id_ov * 65536 + event_id;
	               end
	            else
	               eventIDTracker = event_id_ov * 65536 + event_id;
	            end
    
                %find corresponding tracker entry and save XY info
                trackIdx = find(TR_DATA_eventIDColumn == eventIDTracker);
                if trackIdx>0
                    %is valid trackIdx -> get XY
                    xPos = TR_DATA_xPosColumn(trackIdx); %TR_DATA(trackIdx,dutColArray(TR_DUTINDEX,2));
                    yPos = TR_DATA_yPosColumn(trackIdx); %TR_DATA(trackIdx,dutColArray(TR_DUTINDEX,3));
                    if size(xPos,1)>1
                        xPos = xPos(1);
                    end
                    if size(yPos,1)>1
                        yPos = yPos(1);
                    end
                    trackerSignalValid = true;
                else
                    str_disp=sprintf('No tracking data found for event %d', eventIDTracker);
                    %disp(str_disp);
                end
            end
            
            %% check if DUT, REF and tracker contain valid signals
            if ff ==1 && m<saveSignalWaveformsNumber
            shouldSave = 1;
            else 
            shouldSave = 0;
            end


            if(minuit==1)
                MM_temp = process_signal_minuit(t_vec_mm-etime,MM_Y(:,m),opts_MM,shouldSave,storeFolderSignalDUT,eventIDTracker);
            else
                MM_temp = process_signal(t_vec_mm-etime,MM_Y(:,m),opts_MM);
            end
            if(MM_temp.fail==0)
                dutSignalValid = true;
                localMaxY(m) = MM_temp.sig.max.y;
            else
                localMaxY(m) = 0;
            end
            
            if(minuit==1)
                MCP_temp = process_signal_minuit(t_vec_mcp-etime,MCP_Y(:,m),opts_MCP,shouldSave,storeFolderSignalRef,eventIDTracker);
            else
                MCP_temp = process_signal(t_vec_mcp-etime,MCP_Y(:,m),opts_MCP);
            end

            MCP_temp.event_id = eventIDTracker;
            MM_temp.event_id = eventIDTracker;
                        
            if(MCP_temp.fail==0)
                refSignalValid = true;
                localMaxYMCP(m) = MCP_temp.sig.max.y;
            else
                localMaxYMCP(m) = 0;
            end
            
         
            localValidREF(m) = refSignalValid;
            localValidDUT(m) = dutSignalValid;
            localValidTracker(m) = trackerSignalValid;
            localTrackerX(m) = xPos;
            localTrackerY(m) = yPos;

            if enableTracker == 1
                local_eventIDArray(m) = eventIDTracker;
            else 
                local_eventIDArray(m) = 0;
            end
	    
            
            
            %% checked if signals valid


            %if valid, save
            if refSignalValid && dutSignalValid && trackerSignalValid
                
                % store valid data into array of structures
                local_time_diff(m) = MM_temp.cfd.time-MCP_temp.cfd.time;
		        local_time_diff_sigmoid(m) = MM_temp.sigmoid.timepoint-MCP_temp.sigmoid.timepoint;

                %localMMData = [localMMData; MM_temp]; 
                %localMCPData = [localMCPData; MCP_temp];
                local_e_peak_MM(m) = MM_temp.sig.charge.e_peak;
                local_e_peak_MCP(m) =  MCP_temp.sig.charge.e_peak;
                local_e_lead_MM(m) =  MM_temp.sig.charge.lead_edge;
                local_riseTime(m) = MM_temp.sigmoid.timepoint90-MM_temp.sigmoid.timepoint10;

            else
    		    local_time_diff(m) = 0;
		        local_time_diff_sigmoid(m) = 0;
                local_e_peak_MM(m) = 0;
                local_e_peak_MCP(m) =  0;
                local_e_lead_MM(m) =  0;
                local_riseTime(m) = 0;
		        localMaxY(m) = 0;
                localMaxYMCP(m) = 0;
            end
            
            
        end
        %toc

	    %make them logical arrays
	    localValidDUT = logical(localValidDUT);
	    localValidREF = logical(localValidREF);
	    localValidTracker = logical(localValidTracker);
    
	    eventsValidDUT_all{ff} = localValidDUT;
        eventsValidREF_all{ff} = localValidREF;
        eventsValidTracker_all{ff} = localValidTracker;
    
	    valid = logical(localValidDUT & localValidREF & localValidTracker);

        %saving data from each file in this cell array
        %MM_data_onlyValid{ff} = localMMData; % --->this will go in MM_data
	    MM_maxy_all{ff} = localMaxY; % --->this will go in eventsDUTmaxY
        MM_maxy_onlyValid{ff} = localMaxY(valid); % --->this will go in MM_maxy
        
        e_peak_MM_onlyValid{ff} = local_e_peak_MM(valid); % ---> this will go in e_peak_MM
        e_lead_MM_onlyValid{ff} = local_e_lead_MM(valid); % ---> this will go in e_lead_MM
        riseTime_onlyValid{ff} = local_riseTime(valid); % ---> this will go in riseTime

	    %MCP_data_onlyValid{ff} = localMCPData;	% --->this will go in MCP_data
	    MCP_maxy_all{ff} = localMaxYMCP;
        MCP_maxy_onlyValid{ff} = localMaxYMCP(valid); % --->this will go in MCP_maxy

        e_peak_MCP_onlyValid{ff} = local_e_peak_MCP(valid); % ---> this will go in e_peak_MCP

        eventsTrackerX_all{ff} = localTrackerX;
        eventsTrackerY_all{ff} = localTrackerY;
	    eventsTrackerX_onlyValid{ff} = localTrackerX(valid); %---> this will go in trackerX
	    eventsTrackerY_onlyValid{ff} = localTrackerY(valid); %---> this will go in trackerY

	    time_diff_onlyValid{ff} = local_time_diff(valid); %---> this will go in time_diff
    
	    time_diff_sigmoid_onlyValid{ff} = local_time_diff_sigmoid(valid); %---> this will go in time_diff_sigmoid
	    time_diff_sigmoid_all{ff} = local_time_diff_sigmoid; %---> this will go in eventsTimeDiffSigmoid
	    
	    eventIDArray_all{ff} = local_eventIDArray(valid);
        
    end



end

%reunite all the data from each run file
%MM_data = vertcat(MM_data_onlyValid{:});
MM_maxy = horzcat(MM_maxy_onlyValid{:});
eventsDUTmaxY = horzcat(MM_maxy_all{:});

e_peak_MM = horzcat(e_peak_MM_onlyValid{:});
e_lead_MM = horzcat(e_lead_MM_onlyValid{:});
riseTime = horzcat(riseTime_onlyValid{:});

%MCP_data = vertcat(MCP_data_onlyValid{:});
MCP_maxy = horzcat(MCP_maxy_onlyValid{:});
eventsREFmaxY = horzcat(MCP_maxy_all{:});

e_peak_MCP = horzcat(e_peak_MCP_onlyValid{:});

eventsValidDUT = horzcat(eventsValidDUT_all{:});
eventsValidREF = horzcat(eventsValidREF_all{:});
eventsValidTracker = horzcat(eventsValidTracker_all{:});

eventsTrackerX = horzcat(eventsTrackerX_all{:});
eventsTrackerY = horzcat(eventsTrackerY_all{:});
trackerX = horzcat(eventsTrackerX_onlyValid{:});
trackerY = horzcat(eventsTrackerY_onlyValid{:});

time_diff = horzcat(time_diff_onlyValid{:});
time_diff_sigmoid = horzcat(time_diff_sigmoid_onlyValid{:});
eventsTimeDiffSigmoid = horzcat(time_diff_sigmoid_all{:});

eventIDArray = horzcat(eventIDArray_all{:});


%% save scope settings
run.scope_set_ch_mcp = ch_mcp.info;
run.scope_set_ch_mm = ch_mm.info;
run.scope_set_ch_tr = ch_tr.info;

%% save data to mat file
if shouldSaveMAT
    save(['\\eosproject-smb\eos\project\p\picosec\testbeam\2025_July_h4\Matlab\Run' run.id '-' run.oscilloscope '-' runInfoString '.mat'], 'time_diff_sigmoid', 'MCP_maxy', 'MM_maxy', 'trackerX', 'trackerY', 'eventsDUTmaxY', 'e_peak_MM', 'e_lead_MM', 'riseTime', 'eventsREFmaxY', 'e_peak_MCP', 'eventsValidDUT', 'eventsValidREF', 'eventsValidTracker', 'eventsTrackerX', 'eventsTrackerY', 'eventsTimeDiffSigmoid', 'eventIDArray');
end
toc

%AnalyseRun