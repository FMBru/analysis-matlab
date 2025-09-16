%VERSION DATE: 30/06/2025 - FRANCESCO (summer student)
%code used to analyze single photoelectron response of PICOSEC using UV LED
% there's no timing info, used only to check if it possible to extract the
% number of p.e. emitted knowing the single p.e. response when using the same DUT 

close all
addpath '\\eosproject-smb\eos\project\p\picosec\testbeam\2025_July_h4\Matlab\CommonFunctions';
addpath '\\eosproject-smb\eos\project\p\picosec\testbeam\2025_July_h4\Matlab\Analysis\Functions';


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
    numberFilesToAnalyse = noFileToAnalyze;
    opts_MM.wholeIonTail = enableIonTail;
    runInfoString = runInfoDesc;
else
    %not batch processing - List here Channels and DUT if not using batch
    %processing
    clear all
    run.lecroy_name = '----Trace----'; %['Run' run.id];
    opts_MM.en_noiseRejection = 1; %enable strong noise rejection
    opts_MM.en_filter= 1;    %enable filtering 
    run.id = ['013'];
    run.oscilloscope = 'Pool2';
    opts_MM.chID = '4';
    opts_MM.ch_name = ['C' opts_MM.chID run.lecroy_name]; 
    runInfoString = 'MM2-sealed-275-420';
    numberFilesToAnalyse = 0; %max number of files to analyse, 0 -> analyse all in folder
    opts_MM.wholeIonTail = 1; %enable ion tail processing
end

numberDebugPlots = 0;
saveSignalWaveformsNumber = 10;

debugPlotCounter = 0;
run.savedSignals = 0;

% run file coordinates and number of files
minuit = 1;   % set minuit for optimizer if 1 or matlab 0

%run.pathEOS=['\\eosproject-smb\eos\project\p\picosec\testbeam\2025_July_h4\' run.oscilloscope '\Run' run.id '\']; %path where files to analyze are stored
run.pathEOS=['\\eosproject-smb\eos\project\p\picosec\lab\Marta\PhotoDetectorResTi\PhotoDetector_A275V_C' run.id 'V\'];
disp(run.pathEOS)
run.nfiles = find_fileNo(run.pathEOS);



if(numberFilesToAnalyse>0)
    run.nfiles =numberFilesToAnalyse;
end

% 
run.debugPlotsPath=['\\eosproject-smb\eos\project\p\picosec\lab\Marta\PhotoDetectorResTi\Results\PhotoDetector_A275V_C' run.id 'V\DebugPlots\'];%'-'  run.oscilloscope '\'];
% if exist(run.debugPlotsPath, 'dir') == 7
%     rmdir(run.debugPlotsPath, 's');
% end
mkdir(run.debugPlotsPath);

% add_str = '';
% if opts_MM.en_noiseRejection
%     add_str = [add_str '_rej'];
% end
% if opts_MM.en_filter
%     add_str = [add_str '_filt'];
% end
% 
storeFolderSignalDUT = ['\\eosproject-smb\eos\project\p\picosec\lab\Marta\PhotoDetectorResTi\Results\PhotoDetector_A275V_C' run.id 'V\SignalPlots\'];
% 
% if exist(storeFolderSignalDUT, 'dir') == 7
%     rmdir(storeFolderSignalDUT, 's');
% end
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
opts_MM.shouldExcludeSaturatedSignal = 1; %maximum value that can be registered with the oscilloscope in certain acquisition settings

%% DO PROCESSING

if isempty(gcp('nocreate'))
    parpool('local', 4); %activate a parallelization pool
end

%tic;         % start time measurement
k=1;         % valid data counter
j=1;         % all data counter

% Pre-allocation
MM_data_all = cell(run.nfiles, 1); 
MM_maxy_all = cell(run.nfiles, 1);
eventsValidDUT_all = cell(run.nfiles, 1);

%broadcasting variabless for the parfor loop
options = opts_MM;
chName = opts_MM.ch_name;
pathEOS = run.pathEOS;

titleStr = {'Saturated'; 'tPrerms'; 'sigBelow3RMS'; 'newNoiseRej'; 'ePeakFinderOld'; 'ePeakWidthNeg'; 'SmoothedRangeTooHigh'; 'newFiltRej'; 'NoSmoothedPeak'; 'ePeakEndOutside'};

eventIDTracker = 0;
shouldSave = 0;

for ff=1:run.nfiles 
    
    %temporary storage for each file
%     localData = []; 
%     localMaxY = [];
%     localValid = [];

    
    %dummy struct
    ch_mm = struct();
    
    % display current file number
    str_disp=sprintf('Loading file set No. %d', ff);
    disp(str_disp);
    
    ch_mm_str=sprintf('%s%s%05d.trc', pathEOS, chName, ff);
    disp(ch_mm_str);

    filesExists = 0;
    
    if exist(ch_mm_str,'file')==2 
        
        try
        filesExists = 1;
        
        % decode binary files
        ch_mm = ReadLeCroyBinaryWaveform(ch_mm_str);
        
        str_disp = sprintf('Processing file set No. %d', ff);
        %disp(str_disp);
        catch ME
            disp(['Errore nella lettura del file ' ch_mm_str]);
            disp(['Errore: ' ME.message]);
            filesExists = 0;
        end
%         filesExists=1;
%         
%         % decode binary files
%         ch_mm = ReadLeCroyBinaryWaveform(ch_mm_str);
%         
%         str_disp=sprintf('Processing file set No. %d', ff);
%         disp(str_disp);
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

        %temporary storage for each file
        localData = []; 
        localMaxY = zeros(1,nTRCseg);
        localValid = zeros(1,nTRCseg);
    
        %tic
        % go trough all of the events in the file
        parfor m=1:nTRCseg
%         for m=1:nTRCseg
            
            %disp(['Segment' int2str(m)]);
            %unique ID for the event
            eventIDTracker = (ff-1)*nTRCseg + m - 1 ;
             
            % generate virtual time vector (important to take care for trigger
            % offset in LeCroy scope)
            t_vec_mm=(0:lTRCseg-1)'*Ts; % + ch_mm.trigger_offset(m);
            dummyVector = 1:lTRCseg;
            
            % subtract the earliest time
            etime=min(t_vec_mm);

            %% check if DUT, REF and tracker contain valid signals           
            if ff==1 && m<10
                shouldSave = 1;
            else
                shouldSave = 0;
            end
            
            MM_temp = process_signal_minuitFra(ch_mm.x(:,m),ch_mm.y(:,m),options,shouldSave,storeFolderSignalDUT,eventIDTracker);            
            dutSignalValid = (MM_temp.fail == 0);
            %localValid(end+1) = dutSignalValid;
            localValid(m) = dutSignalValid;
            
            %run.savedSignals = run.savedSignals+1;
            
            %% checked if signals valid
            
            %if valid, save
            %if dutSignalValid          
            MM_temp.waveformX = t_vec_mm;
            MM_temp.waveformY = ch_mm.y(:,m);

            % store valid data into array of structures
            localData = [localData; MM_temp];                
            %localMaxY(end+1) = MM_temp.sig.max.y;
            localMaxY(m) = MM_temp.sig.max.y;

%             else
%                 figure
%                 plot(ch_mm.y(:,m));
%                 if options.en_filter && MM_temp.fail > 6
%                     title(titleStr{MM_temp.fail});
%                     saveas(gcf,[storeFolderSignalDUT '\Noise' int2str(eventIDTracker) '.png']);
%                 end
%                 if options.en_noiseRejection &&  ~options.en_filter && MM_temp.fail == 4
%                     title(titleStr{MM_temp.fail});
%                     saveas(gcf,[storeFolderSignalDUT '\Noise' int2str(eventIDTracker) '.png']);
%                 end
%                 if ~options.en_noiseRejection && ~options.en_filter
%                     title(titleStr{MM_temp.fail});
%                     if (MM_temp.fail == 1)
%                         hold on
%                         %scatter( dummyVector(MM_temp.saturatedIdx),ch_mm.y(MM_temp.saturatedIdx,m), 'r');
%                         title([titleStr{MM_temp.fail} ' ' int2str(sum(MM_temp.saturatedIdx))]);
%                     end
%                     saveas(gcf,[storeFolderSignalDUT '\Noise' int2str(eventIDTracker) '.png']);
%                     hold off
%                 end
%                 
%                 close all;
            %end
        end
        
        %saving data from each file in this cell array
        MM_data_all{ff} = localData;
        MM_maxy_all{ff} = localMaxY;
        eventsValidDUT_all{ff} = localValid;
    end
    %toc
    close all;
end


%reunite all the data from each run file
MM_data = vertcat(MM_data_all{:});
MM_maxy = horzcat(MM_maxy_all{:});
eventsValidDUT = horzcat(eventsValidDUT_all{:});

%toc
% efficiency = (k-1)/eventIDTracker;

AnalyseRunFrancesco