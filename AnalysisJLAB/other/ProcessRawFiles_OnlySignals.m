              % MATLAB/OCTAVE implementation of processing code for PICOSEC
% version date: 7 Nov 2021
    clear all
addpath 'C:\Users\GDD\Documents\MATLAB\Picosec\CommonFunctions';

    
    runIDsToAnalyse = ["500C-350A","485C-360A","475C-360A","480C-360A","490C-360A","500C-360A-bis","480C-350A","490C-350A","450C-350A","460C-350A","470C-350A","470C-360A","520C-350A","515C-350A","510C-350A","475C-350A","500C-360A","450C-385A"];

for runPos = 1:length(runIDsToAnalyse)

runString = convertStringsToChars(runIDsToAnalyse(runPos));

clear MM_data MM_maxy MM_riseTime

close all


enDebugPlots = true;
numberEventsToPlot = 10;

%May22: %1: miniCactus, 2: MM3 (VacChamberPos), 3: MM1(Multipad), 4: MM2 (on support plate), 5: MM4 (electron setup)

    run.lecroy_name = '-Trace-'; %['Run' run.id];
    
    run.id = 'uRWELL';
    run.oscilloscope = 'LabTest';
    opts_MM.chID = '2';
    analysis.dutChannel = opts_MM.chID;
    shouldSaveMAT = false;
    shouldUseEOSFolder = true;
    tracker.dutIndex = 2; %DUT index: 1: MM1 (Multipad), 2: MM2(VacChamber), 3: MM3 (ElectronsetupMM), 4: MM4 (Picolarge main tracker), 5:MM5 (SiPM)
    numberFilesToAnalyse = 0; %max number of files to analyse, 0 -> analyse all in folder

    opts_MM.ch_name = ['C1' run.lecroy_name];


% run file coordinates and number of files
minuit = 1;   % set minuit for optimizer if 1 or matlab 0
run.year = '2022/10 ';
run.name = ['BEAM ' run.year ' RUN ' run.id];
run.path=['C:\Users\GDD\Documents\Picosec\October22\' run.oscilloscope '\Run' run.id '\Run' run.id '\'];

%run.pathEOS=['\\eosproject-smb\eos\project\p\picosec\lab\uRWELL\Data\450C-350A\'];
%store_folder=['\\eosproject-smb\eos\project\p\picosec\lab\uRWELL\Results\450C-350A'];

run.pathEOS=['\\eosproject-smb\eos\project\p\picosec\lab\uRWELL\Data\' runString '\'];
store_folder=['\\eosproject-smb\eos\project\p\picosec\lab\uRWELL\Results\' runString];

store_folderWaveforms=[store_folder '\Waveforms'];

mkdir(store_folder);
mkdir(store_folderWaveforms);

if shouldUseEOSFolder
    run.nfiles = find_fileNo(run.pathEOS)
else
    run.nfiles = find_fileNo(run.path);
end

%override run path
%run.pathEOS = ['\\eosproject-smb\eos\project\p\picosec\testbeam\2022_May_h4\LaserTest\LaserTestJune23\A275C500\'];


if(numberFilesToAnalyse>0)
    run.nfiles =numberFilesToAnalyse;
end

run.lecroy_name = '-Trace-'; %['Run' run.id];


% options for processing micromegas channel (lot magic numbers from
% processing functionhave to be added here)
opts_MM.t_dead = 1;      % blanking time from start of the recording samples
opts_MM.t_prerms = 20;  % blanking time before global maximum samples
opts_MM.t_prerms = 200;  % blanking time before global maximum samples
opts_MM.Ts = 1/20e9;     % sampling speed
opts_MM.Rin=50;          % input impedance
opts_MM.invert = 1;      % is inverted (1-inverted, 0-non-inverted)
opts_MM.type=0;          % detector type (0-DUT, 1-REF) for future use
%only for override - should set in top
%opts_MM.ch_name = ['C4' run.lecroy_name]; % LeCroy file name format MCP1 (used for timing)
opts_MM.en_plot = 0;     % enable debugging plots


%% DO PROCCEISNG
tic;         % start time measurement
k=1;         % valid data counter
event_id_prev = -1;
event_id_ov = 0;


for (ff=1:run.nfiles)
    % display current file number
    str_disp=sprintf('Loading file set No. %d', ff);
    disp(str_disp);
    
    % generate file name strings
    
    
    ch_mm_str=sprintf('%s%s%05d.trc', run.path, opts_MM.ch_name,ff);
    
    
    if shouldUseEOSFolder
        ch_mm_str=sprintf('%s%s%05d.trc', run.pathEOS, opts_MM.ch_name,ff);
    end
    
    
   
    filesExists = 0;
    
    if exist(ch_mm_str,'file')==2 
        filesExists=1;
        
        
        % read files using third party function
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
            t_vec_mm=(0:lTRCseg-1)'*Ts + ch_mm.trigger_offset(m);
            
            % subtract the earliest time
            etime=min(t_vec_mm);
            
            % process MM Picosec first to see if signal is valid
            if(minuit==1)
                MM_temp = process_signal_minuit(t_vec_mm-etime,ch_mm.y(:,m),opts_MM);
            else
                MM_temp = process_signal(t_vec_mm-etime,ch_mm.y(:,m),opts_MM);
            end
            
           if MM_temp.fail==0
               
               
               if k<numberEventsToPlot && enDebugPlots == true
                   figure
                   plot((t_vec_mm-etime)*1e9,-ch_mm.y(:,m));
                   grid on
                   xlim([0 150])

                   grid on
                  xlabel('Time (ns)')
                ylabel('Amplitude (V)');

                   pause(2);
                   saveas(gcf,[store_folderWaveforms '\Event' num2str(k) '.png'])

                   close all;
                   
               end
               
                    
                    
                 % store valid data into array of structures
                    MM_data(k)= MM_temp;
                    MM_maxy(k) = MM_data(k).sig.max.y;
                    MM_riseTime(k)= MM_data(k).sigmoid.timepoint90-MM_data(k).sigmoid.timepoint10;% extract electrom peak charge
                    
                    k=k+1;
           end
                
                  
                end
            end
        end
    
    toc


%% save data to mat file
if shouldSaveMAT
    save(['C:\Users\GDD\Documents\Picosec\May22\Analysed\Run' run.id '-' run.oscilloscope '.mat'], 'run', 'MM_data', 'MCP_data', 'time_diff', 'time_diff_sigmoid', 'MCP_maxy', 'MM_maxy', 'trackerX', 'trackerY');
end
toc

 %% plot rise time distribution
            figure
            
            riseTimeMask = abs(MM_riseTime)<10;
            
            h=histogram(MM_riseTime(riseTimeMask),100);
            hold on
            xbins = h.BinEdges(1:end-1)+h.BinWidth/2;
            
            xlabel('Rise time (ns)')
            ylabel('Events');
            grid on
            xlim([0 10]);
            %legend('MM','MM fit');
            title_str = sprintf('Signal rise time \\mu = %4.4f ns', mean(MM_riseTime(riseTimeMask)));
            title(title_str)
            saveas(gcf,[store_folder '\Run' run.id '_riseTimeHist.png'])
            
    %% plot amplitude distribution
            figure
            h=histogram(MM_maxy,100);
            hold on
            xbins = h.BinEdges(1:end-1)+h.BinWidth/2;
            
            xlabel('Amplitude (V)')
            ylabel('Events');
            grid on
                        xlim([0 1.2]);

            %legend('MM','MM fit');
            title_str = sprintf('Amplitude \\mu = %4.4f V', mean(MM_maxy));
            title(title_str)
            saveas(gcf,[store_folder '\Run' run.id '_ampHist.png'])
            
end
                     