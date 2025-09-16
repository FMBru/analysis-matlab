close all

addpath 'C:\Users\GDD\Documents\MATLAB\Picosec\CommonFunctions';
addpath '..\Functions'

shouldPlotWaveform = false;
shouldSaveWaveformSamples = 5;

numberBGSamples = 100;

%% check if offsets for file numbers exist from oscilloscope saving
offsetTrackerFileNumber = 0;

%% Tracker setup
trackerExist = 0 %set 0 for SPE runs

trigIntervalTimes = [];
    

    
    %%beam test
    run.id = '500V';
    run.year = '2024Lab';

    channel.id = 'C3';
       
    run.name = ['PICOSEC - RUN ' run.id ' - Channel ' channel.id];

%     run.path=['\\eosproject-smb\eos\project\p\picosec\lab\Chiara\amplitudeVsTime\C103\' run.id '\'];
  %  run.path=['\\eosproject-smb\eos\project\p\picosec\testbeam\2023_August_h4\Pool2\Run' run.id '\'];
   % run.path=['D:\2023\picosecSPE_B4C9nm_AntonijaDet\231103\SlowPreamp\' run.id '\'];
         %run.path=['\\eosproject-smb\eos\project\p\picosec\lab\ResistiveMM-SinglePad\ResistiveMM-SPE-SlowPreamp\275A-540C-\signal\'];
         run.path=['\\eosproject-smb\eos\project\p\picosec\lab\ResistiveMM-SinglePad\RMM-SPE-SlowPream-060624\10mVdiv\325A-500C\signal\'];


%     run.nfiles = find_fileNo(run.path);
     run.lecroy_name = '--Trace--'; %['Run' run.id];

    run.nfiles =10; % override number of file sets to process (see how many trc files are in the folder)
%     run.lecroy_name = '_'; %['Run' run.id];
   


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %store_folderWaveforms = ['\\eosproject-smb\eos\project\p\picosec\testbeam\2023_July_h4\Results\SPE\Gauss+Polya\Run' run.id '\Waveforms\'];
            %mkdir(store_folderWaveforms);

    
    %%%%%%%%%
    
    str_disp=sprintf('- - - - - ANALYSING RUN%s - - - - - ', run.id);
    disp(str_disp);

    %option for subtracting noise from signal
    %noiseFile = 305;


    shiftX = 0;
    shiftY = 0;

    if trackerExist == 1
        tracker.path = ['\\eosproject-smb\eos\project\p\picosec\testbeam\2023_April_h4\tracker\reconstructed\asciiRun' run.id '.dat'];
        tracker.en = 1; %match eventIDs to tracking data and add XY to output
        
        tracker.dutIndex = 2; %DUT index: 1: MM1 (Multipad), 2: MM2(VacChamber), 3: MM3 (ElectronsetupMM), 4: MM4 (Picolarge main tracker), 5:MM5 (SiPM)

        %dutColArray = [1 4 5; 2 7 8; 3 10 11; 4 13 14; 5 16 17]; %[dID colXID colYID; ] -> May MM run
        %dutColArray = [1 10 11; 2 16 17; 3 4 5; 4 13 14; 5 7 8]; %[dID colXID colYID; ] -> July 2022 MM run
        dutColArray = [1 10 11; 2 13 14; 3 19 20; 4 22 23; 5 4 5; 6 25 26; 7 16 17] %July 2023
        %DUT index: 1: MM1 (Multipad), 2: MM2(VacChamber), 3: MM3
        %(ElectronsetupMM), 4: MM4 (Picolarge main tracker), 5:MM5 (SiPM) July 2022

        if tracker.en
            trackerFile = fopen(tracker.path,'rt');
            D = textscan(trackerFile, '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f', 'Delimiter','\t', 'HeaderLines',2, 'CollectOutput',1);
            tracker.data = cell2mat(D);
        end

        % options for extracting tracker ID from the bitstream
        opts_TR.ch_name = ['C3' run.lecroy_name];
        opts_TR.baud_rate = 40e6; % 40Mbps baudrate
        opts_TR.n_bits = 16;      % number of bits after start bit

    end
    %%
    opts_MM.ch_name = [channel.id run.lecroy_name]; % LeCroy file name format

    opts_MM.en_plot = 0;     % enable debugging plots

    event_id_ov = 0;

    mm_max_y = [];
    mm_bgLevel = [];
    mm_bgRMS = [];
    noise_mm_max_y = [];
    int_vector = [];

    k=1; 
    event_id_prev = -1;
    
  %  run.nfiles = 10; %run only part of files
    for (ff=0:run.nfiles)
        % generate file name strings
        ch_mm_str=sprintf('%s%s%05d.trc', run.path, opts_MM.ch_name,ff);

        if trackerExist == 1
            ch_tr_str=sprintf('%s%s%05d.trc', run.path, opts_TR.ch_name,ff+offsetTrackerFileNumber); 
        end

        % check if the file number exists, if not, go to next number
        if not(isfile(ch_mm_str))
            continue
        end

        % display current file number
        str_disp=sprintf('Loading file set No. %d', ff);
        disp(str_disp);


        %filesExists = 0;
        % if exist(ch_mm_str,'file')==2 && exist(ch_tr_str,'file')==2
        %      filesExists=1;

        % read files using third party function
        ch_mm = ReadLeCroyBinaryWaveform(ch_mm_str);
            
        if trackerExist == 1
            ch_tr = ReadLeCroyBinaryWaveform(ch_tr_str);
        end

        str_disp=sprintf('Processing file set No. %d', ff);
        disp(str_disp);

        %else
        %    str_disp=sprintf('Skipping file set No. %d', ff);
        %    disp(str_disp);

        %end

        % get number of segments (events) in one file
        nTRCseg = ch_mm.info.nbSegments;
        % get segment length
        lTRCseg = size(ch_mm.y,1);
        % calculate sampling time
        Ts = ch_mm.x(2,1) - ch_mm.x(1,1);
        lastTrigTime = 0;
     
        for (i=1:nTRCseg)
            trigTime =  ch_mm.x(1,i);
            if i>1
                trigInterval = lastTrigTime-trigTime;
                trigIntervalTimes = [trigIntervalTimes;trigInterval];

            end
            

            lastTrigTime = trigTime;

            if trackerExist == 1
            t_vec_mm=(0:lTRCseg-1)'*Ts + ch_mm.trigger_offset(i);
       maxTrackerAmp = max(-ch_tr.y(:,i));
       sumTrackerAmp = sum(-ch_tr.y(:,i));
                if tracker.en %&& maxTrackerAmp>0.8 && sumTrackerAmp>500
                    
                    % process tracker ID channel
                    event_id = process_tr_bitstream(t_vec_mm, ch_tr.y(:,i), opts_TR);
%                      plot(t_vec_mm, ch_tr.y(:,i));
%                     pause(2);
%                      close all
                    % count overflows
                    if(event_id < event_id_prev)
                        event_id_ov = event_id_ov + 1
                                             plot(t_vec_mm, ch_tr.y(:,i));
                    pause(5);
                     close all

                    end
                    event_id_prev = event_id;

                    MM_temp.event_id = event_id_ov * 65536 + event_id;

                    xPos = 0;
                    yPos = 0;

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


                    MM_temp.x = xPos;
                    MM_temp.y = yPos;
                    trackerX(k) = MM_temp.x;
                    trackerY(k) = MM_temp.y;
                    eventIDArray(k) = MM_temp.event_id;
                end
            
            
            % store valid data into array of structures
            MM_data(k)= MM_temp;
            trackerX(k) = MM_temp.x;
            trackerY(k) = MM_temp.y;
        
            k=k+1;
            end
            
            
            % entire event signal
            eventXvec = ch_mm.x(1:end-1,i);
            eventYvec = -ch_mm.y(1:end-1,i);
            
            %extract only around peak
%             lengthSample = 250; %number points to use before and after peak
%             [max_y,iMaxY] = max(eventYvec);
%             eventXvec = ch_mm.x((iMaxY-lengthSample:iMaxY+lengthSample),i);
%             eventYvec = -ch_mm.y((iMaxY-lengthSample:iMaxY+lengthSample),i);

            
            
            bg_level = mean(eventYvec(1:numberBGSamples));
            eventYvec = eventYvec - bg_level ;             
            bg_levelCorrected = mean(eventYvec(1:numberBGSamples));
            bg_rms = std(eventYvec(1:numberBGSamples));


            [max_y,iMaxY] = max(eventYvec);
            max_position = eventXvec(iMaxY);
            
            if shouldPlotWaveform
               plot(eventXvec, eventYvec,'k'); hold on
               plot(eventXvec(1:numberBGSamples), eventYvec(1:numberBGSamples),'r'); 
               plot(eventXvec(iMaxY), max_y,'.r');  hold off
               pause(1); 
               close all;
            end
            

            if shouldSaveWaveformSamples >0 && length(mm_bgLevel)<shouldSaveWaveformSamples
              
               plot(eventXvec-min(eventXvec), eventYvec); 
               ylim([-0.005 0.1] );
               pause(1);
               %saveas(gcf,[store_folderWaveforms 'RUN' run.id ' - Channel' channel.id '- waveform' int2str(length(mm_bgLevel)) '.png'])
               close all
            end

            
%                         %Fourier transformation for the first events
%                         if shouldSaveWaveformSamples >0 && length(mm_bgLevel)<shouldSaveWaveformSamples-3
%               
%                             y_ft = fft(eventYvec);
%                             fs = 1/Ts;
%                             f = (0:length(y_ft)-1)*fs/length(y_ft);
% 
% %                             n_len=length(eventXvec-min(eventXvec));
% %                             fshift = (-n_len/2:n_len/2-1)*(fs/n_len);
% %                             yshift=fftshift(y_ft);
% %                             plot(fshift,abs(yshift))
%                          
%                             plot(f, abs(y_ft));
%                             %xlim([-0.05e9 5.5e9] );
%                             set(gca, 'YScale', 'log') %plot on log scale
%                             ylim([0 1]);
%                             xlabel('Frequency (Hz)')
%                             ylabel('Magnitude')
%                             title('Fourier Transformation')
%                pause(5);
% close all
%                         end
            

            % peak is 200 events before max to 200 events after max
            iStart = iMaxY - 200;
            iEnd = iMaxY + 200;

            % remove noise for curve fit
            %if (max_y<3E-3)
            %    continue
            %end

            % deal with edge cases where whole peak isn't in the array
            if (iStart < 1)
                iStart = 1;
            end

            if (iEnd > length(eventXvec))
                iEnd = length(eventXvec);
            end

            % this is the part of the signal that makes the peak
            iEventXvec = eventXvec(iStart:iEnd);
            iEventYvec = eventYvec(iStart:iEnd);


            %% saving to data vectors
            int_iEventYvec = sum(iEventYvec); % sum of the amps in the peak. Histogram of the integral of the peak?
            int_vector= [int_vector; int_iEventYvec];
            mm_max_y = [mm_max_y; max_y]; %concatenate onto mm_max_y for each file
            mm_bgLevel = [mm_bgLevel; bg_level]; 
            mm_bgRMS = [mm_bgRMS; bg_rms]; 

            close all

%                    if (max_y > 3e-2 && max_y < 3.6e-2)
%                         hold on
%                         plot(eventXvec, eventYvec)
%                         max_y
%                         plot (iEventXvec,iEventYvec)
%                         pause(2);
%                   end
            %ret = exist(noiseFile);
            %print('ret: %d', ret);

            if ( exist('noiseFile') && (ff >= noiseFile ))
                debug = sprintf('Noise File No. %d', ff);
                disp(debug);

                % entire event signal
                noise_eventXvec = ch_mm.x(1:end-1,i);
                noise_eventYvec = -ch_mm.y(1:end-1,i);

                bg_level = mean(noise_eventYvec(1:20));
                noise_eventYvec = noise_eventYvec - bg_level ;

                %don't need max value
                [noise_max_y,noise_iMaxY] = max(noise_eventYvec);

                max_position = noise_eventXvec(iMaxY);

                % peak is 200 events before max to 200 events after max
                %iStart = iMaxY - 200;
                %iEnd = iMaxY + 200;
                % remove noise for curve fitnoiseFile
                %if (max_y<3E-3)
                %    continue
                %end
                %if (iStart < 1)
                %    iStart = 1;
                %end
                %if (iEnd > length(eventXvec))
                %    iEnd = length(eventXvec);
                %end
                %iEventXvec = eventXvec(iStart:iEnd);
                %iEventYvec = eventYvec(iStart:iEnd);


                %int_iEventYvec = sum(iEventYvec);
                %int_vector= [int_vector; int_iEventYvec];
                noise_mm_max_y = [noise_mm_max_y; noise_max_y];

                close all
            end
        end
    end

    meanTimeInterval = mean(trigIntervalTimes);
    meanRate = 1/meanTimeInterval

   % AnalyseRun_PEAnalysis_noise_plus_sig
