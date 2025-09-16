%VERSION DATE: 30/06/2025 - FRANCESCO (summer student)
%code used to analyze single photoelectron response of PICOSEC using UV LED
% there's no timing info, used only to check if it possible to extract the
% number of p.e. emitted knowing the single p.e. response when using the same DUT 

close all
addpath '\\eosproject-smb\eos\project\p\picosec\testbeam\2024_September_h4\Matlab\CommonFunctions';
addpath '\\eosproject-smb\eos\project\p\picosec\testbeam\2024_September_h4\Matlab\Analysis\Functions';

%not used now --
enDebugPlots = true;
correctFilenumberOffset = false;

shouldPlotOutliers = false;
shouldPlotTrackerWaveform = false;
% --

% clear everything if not batch processing
%set run parameters
if exist('batchProcess','var') == 1
    run.id = runIDString;
    run.oscilloscope = oscilloscopeString;
    run.lecroy_name = '----Trace----'; %['Run' run.id];
    opts_MM.chID = DUTnameString;
    opts_MM.ch_name = ['C' DUTnameString run.lecroy_name];
    opts_MM.en_noiseRejection = enableNoiseRej; %enable strong noise rejection
    opts_MM.en_filter= enableFilter;    %enable filtering 
    shouldSaveMAT = false;
    shouldUseEOSFolder = useEOSFolder;
    tracker.dutIndex = str2num(trackerDUTIndex);  %DUT index: 1: MM1 (Multipad), 2: MM2(VacChamber), 3: MM3 (ElectronsetupMM), 4: MM4 (Picolarge main tracker), 5:MM5 (SiPM)
    numberFilesToAnalyse = noFileToAnalyze;
    runInfoString = runInfoDesc;
else
    %not batch processing - List here Channels and DUT if not using batch
    %processing
    clear all
    run.lecroy_name = '----Trace----'; %['Run' run.id];
    opts_MM.en_noiseRejection = 0; %enable strong noise rejection
    opts_MM.en_filter= 0;    %enable filtering 
    run.id = ['0807_Ti_C100_A275_C470_newCSA_spectrum'];
    run.oscilloscope = 'Pool1';
    opts_MM.chID = '3';
    analysis.dutChannel = opts_MM.chID; 
    opts_MM.ch_name = ['C' opts_MM.chID run.lecroy_name]; 
    shouldSaveMAT = false;
    shouldUseEOSFolder = true;
    tracker.dutIndex = 1; %DUT index: 1: MM1 (Multipad), 2: MM2(VacChamber), 3: MM3 (ElectronsetupMM), 4: MM4 (Picolarge main tracker), 5:MM5 (SiPM)
    numberFilesToAnalyse = 0; %max number of files to analyse, 0 -> analyse all in folder
    runInfoString =  run.id;
end

numberDebugPlots = 10;
saveSignalWaveformsNumber = 0;%numberDebugPlots;

debugPlotCounter = 0;
run.savedSignals = 0;


% run file coordinates and number of files
minuit = 1;   % set minuit for optimizer if 1 or matlab 0

run.year = '2025 June';
run.name = ['LAB ' run.year ' RUN ' run.id];
run.pathEOS=['\\eosproject-smb\eos\project\p\picosec\lab\Francesco\UV_LED_labtest\metallic_singlechannel\sealedMode\newCSA\' run.id '\']; %path where files to analyze are stored
disp(run.pathEOS)
run.nfiles = find_fileNo(run.pathEOS);

if(numberFilesToAnalyse>0)
    run.nfiles =numberFilesToAnalyse;
end


run.debugPlotsPath=['\\eosproject-smb\eos\project\p\picosec\lab\Francesco\UV_LED_labtest\metallic_singlechannel\DebugPlots\sealedMode\' run.id];
if exist(run.debugPlotsPath, 'dir') == 7
    rmdir(run.debugPlotsPath, 's');
end
mkdir(run.debugPlotsPath);

% add_str = '';
% if opts_MM.en_noiseRejection
%     add_str = [add_str '_rej'];
% end
% if opts_MM.en_filter
%     add_str = [add_str '_filt'];
% end

storeFolderSignalDUT = ['\\eosproject-smb\eos\project\p\picosec\lab\Francesco\UV_LED_labtest\metallic_singlechannel\SignalPlots\' run.id];% add_str];

if exist(storeFolderSignalDUT, 'dir') == 7
    rmdir(storeFolderSignalDUT, 's');
end
mkdir(storeFolderSignalDUT);

% options for processing micromegas channel (lot magic numbers from
% processing function have to be added here)
opts_MM.t_dead = 1;      % blanking time from start of the recording samples
opts_MM.t_prerms = 20;   % blanking time before global maximum samples
opts_MM.Ts = 1/10e9;     % sampling speed
opts_MM.Rin=50;          % input impedance
opts_MM.invert = 1;      % is inverted (1-inverted, 0-non-inverted)
opts_MM.type=0;          % detector type (0-DUT, 1-REF) for future use
opts_MM.en_plot = 0;     % enable debugging plots
% opts_MM.en_noiseRejection = 0; %enable strong noise rejection
% opts_MM.en_filter= 0;    %enable filtering 
opts_MM.shouldExcludeWidePeaks = 0; %when a signal is too wide (width = time between end of the start and the end of e-peak) could be a multipeak 
%opts_MM.saturation_value = 0.075; %maximum value that can be registered with the oscilloscope in certain acquisition settings

%% DO PROCCEISNG
tic;         % start time measurement
k=1;         % valid data counter
j=1;         % all data counter

eventIDTracker = 0;
for (ff=1:run.nfiles)      
    % display current file number
    str_disp=sprintf('Loading file set No. %d', ff);
    disp(str_disp);
    
    ch_mm_str=sprintf('%s%s%05d.trc', run.pathEOS, opts_MM.ch_name,ff);
    disp(ch_mm_str);

    filesExists = 0;
    
    if exist(ch_mm_str,'file')==2 
        filesExists=1;
        
        % decode binary files
        ch_mm = ReadLeCroyBinaryWaveform(ch_mm_str);
        
        str_disp=sprintf('Processing file set No. %d', ff);
        disp(str_disp);
    else
        str_disp=sprintf('Skipping file set No. %d', ff);
        disp(str_disp);
    end
    
    if filesExists==1
        % get number of segments (events) in one file
        nTRCseg = ch_mm.info.nbSegments;
        % get segment length
        lTRCseg = size(ch_mm.y,1);
        % calculate sampling time
        Ts = ch_mm.x(2,1) - ch_mm.x(1,1);
        
        % go trough all of the events in the file
        for m=1:nTRCseg
            % generate virtual time vector (important to take care for trigger
            % offset in LeCroy scope)
            t_vec_mm=(0:lTRCseg-1)'*Ts; % + ch_mm.trigger_offset(m);
            
            % subtract the earliest time
            etime=min(t_vec_mm);
            
            dutSignalValid = false;

            MM_temp.x = 0;
            MM_temp.y = 0;

            %% check if DUT, REF and tracker contain valid signals
            if run.savedSignals<saveSignalWaveformsNumber
                shouldSave = 1;
                run.savedSignals = run.savedSignals +1;
            else
                shouldSave = 0;
            end
            
            
            MM_temp = process_signal_minuit(ch_mm.x(:,m),ch_mm.y(:,m),opts_MM,run,shouldSave,storeFolderSignalDUT,eventIDTracker);
            

            eventIDTracker = eventIDTracker +1 ;
            
            if(MM_temp.fail==0)
                dutSignalValid = true;
            end
            
            eventsValidDUT(j) = dutSignalValid;
            
            
            
            
            %% checked if signals valid
            
            %if valid, save
            if dutSignalValid  %refSignalValid

                if debugPlotCounter<numberDebugPlots
                    figure(2)
                    plot(t_vec_mm-etime,ch_mm.y(:,m));
                    title('DUT / REF');
                    ylabel('Voltage, V');
                    xlabel('Time, ns');
                    legend('DUT');               
                 
                    %saveDebugPlot
                    
                    saveas(gcf,[run.debugPlotsPath '\Run' run.id '_' int2str(debugPlotCounter) '_signals.png'])
                    pause(1);
                    close all;
                    debugPlotCounter = debugPlotCounter+1;
                   
                 
                end
                
                MM_temp.waveformX = t_vec_mm;
                MM_temp.waveformY = ch_mm.y(:,m);
                
                % store valid data into array of structures
                MM_data(k)= MM_temp;    % has the virtual time (s) vector in waveform and real time info (ns) in sigmoid
                time_diff(k) = MM_data(k).cfd.time;%-MCP_data(k).cfd.time;
                time_diff_sigmoid(k) = MM_data(k).sigmoid.timepoint;%-MCP_data(k).sigmoid.timepoint;
                MM_maxy(k) = MM_data(k).sig.max.y;
                
                k=k+1;
            end
            
            j = j+1;
        end
        %j = j+1;
    end
    toc
    close all;
end

% efficiency = (k-1)/eventIDTracker;

%% save scope settings
%run.scope_set_ch_mcp = ch_mcp.info;
%run.scope_set_ch_mm = ch_mm.info;

%% save data to mat file
% if shouldSaveMAT
%     save(['C:\Users\GDD\Documents\Picosec\Apr23\Analysed\Run' run.id '-' run.oscilloscope '.mat'], 'run', 'MM_data', 'MCP_data', 'time_diff', 'time_diff_sigmoid', 'MCP_maxy', 'MM_maxy', 'trackerX', 'trackerY');
% end
% toc

AnalyseRun_LED