clear all
close all


run.id = '780';
run.year = '2022/07';
run.name = ['PICOSEC - single PE ' run.year ' - RUN ' run.id];
run.path=['C:\Users\Michaela\Documents\4th Year\CERN\MATLAB\MM2 DLC\Run' run.id '\'];
run.nfiles = find_fileNo(run.path);
%run.nfiles =200; % override number of file sets to process (see how many trc files are in the folder)
run.lecroy_name = '--Trace--'; %['Run' run.id];

% option for subtracting noise from signal
noiseFile = 305;

opts_MM.ch_name = ['C4' run.lecroy_name]; % LeCroy file name format

opts_MM.en_plot = 1;     % enable debugging plots


    mm_max_y = [];
    noise_mm_max_y = [];
    int_vector = [];

for (ff=250:run.nfiles)
    % display current file number
    str_disp=sprintf('Loading file set No. %d', ff);
    disp(str_disp);
    
    % generate file name strings
    ch_mm_str=sprintf('%s%s%05d.trc', run.path, opts_MM.ch_name,ff);
    
    % read files using third party function
    ch_mm = ReadLeCroyBinaryWaveform(ch_mm_str);
    
    str_disp=sprintf('Processing file set No. %d', ff);
    disp(str_disp);
    
    % get number of segments (events) in one file
    nTRCseg = ch_mm.info.nbSegments;
    % get segment length
    lTRCseg = size(ch_mm.y,1);    
    % calculate sampling time
    Ts = ch_mm.x(2,1) - ch_mm.x(1,1);
    
    
    for (i=1:nTRCseg)
       % entire event signal
       eventXvec = ch_mm.x(1:end-1,i);
       eventYvec = -ch_mm.y(1:end-1,i);
       
       bg_level = mean(eventYvec(1:20));
       eventYvec = eventYvec - bg_level ;
       
       [max_y,iMaxY] = max(eventYvec);
       
       max_position = eventXvec(iMaxY);
       
       % peak is 200 events before max to 200 events after max
       iStart = iMaxY - 200;
       iEnd = iMaxY + 200;
       
       % remove noise for curve fit
       %if (max_y<3E-3)
       %    continue
       %end
       
       if (iStart < 1)
           iStart = 1;
       end    
       
       if (iEnd > length(eventXvec))
           iEnd = length(eventXvec);
       end
       
       % this is the part of the signal that makes the peak
       iEventXvec = eventXvec(iStart:iEnd);
       iEventYvec = eventYvec(iStart:iEnd);
       
       
       int_iEventYvec = sum(iEventYvec); % sum of the amps in the peak. Histogram of the integral of the peak?
       int_vector= [int_vector; int_iEventYvec];
       mm_max_y = [mm_max_y; max_y]; %concatenate onto mm_max_y for each file
       
       close all
      
%        if (max_y > 3e-3)
%             hold on
%             plot(eventXvec, eventYvec)
%             max_y
%             plot (iEventXvec,iEventYvec)
%             pause(2);
%       end
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
       
       % remove noise for curve fit
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


Amplitude2_histfix