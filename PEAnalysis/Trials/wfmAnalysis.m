clear all
close all

addpath 'C:\Users\GDD\Documents\MATLAB\Picosec\CommonFunctions';
addpath '..\Functions'

shouldPlotWaveform = false;
shouldSaveWaveformSamples = 0;

numberBGSamples = 10;

% Check if offsets for file numbers exist from oscilloscope saving
offsetTrackerFileNumber = 0;

%% Tracker setup
trackerExist = 0; %set 0 for SPE runs, and lab
% Define run values
cathodeV = 500;
anodeV = 295;
mVdiv = 10;

run.id = sprintf('%dV', cathodeV);
run.year = '2024Lab';   
channel.id = 'C3';
run.name = ['PICOSEC - RUN ' run.id ' - Channel ' channel.id];
run.path = sprintf('\\\\eosproject-smb\\eos\\project\\p\\picosec\\lab\\ResistiveMM-SinglePad\\RMM-SPE-SlowPream-Elec-100624\\%dmVdiv\\%dA-%dC\\noise\\', mVdiv, anodeV, cathodeV);
run.lecroy_name = '--Trace--'; %['Run' run.id];
run.nfiles =15; % override number of file sets to process (see how many trc files are in the folder)
store_folder = sprintf('\\\\eosproject-smb\\eos\\project\\p\\picosec\\lab\\ResistiveMM-SinglePad\\RMM-SPE-SlowPream-Elec-100624\\%dmVdiv\\%dA-%dC\\Result\\', mVdiv, anodeV, cathodeV);
mkdir(store_folder);

str_disp=sprintf('- - - - - ANALYSING RUN %s - - - - - ', run.id);
disp(str_disp);

shiftX = 0;
shiftY = 0;

if trackerExist == 1
    tracker.path = ['\\eosproject-smb\eos\project\p\picosec\testbeam\2023_July_h4\tracker\reconstructed\asciiRun' run.id '.dat'];
    tracker.en = 1; %match eventIDs to tracking data and add XY to output
    
    tracker.dutIndex = 2; %DUT index: 1: MM1 (Multipad), 2: MM2(VacChamber), 3: MM3 (ElectronsetupMM), 4: MM4 (Picolarge main tracker), 5:MM5 (SiPM)

    %dutColArray = [1 4 5; 2 7 8; 3 10 11; 4 13 14; 5 16 17]; %[dID colXID colYID; ] -> May MM run
    dutColArray = [1 10 11; 2 16 17; 3 4 5; 4 13 14; 5 7 8]; %[dID colXID colYID; ] -> July 2022 MM run

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

opts_MM.ch_name = [channel.id run.lecroy_name]; % LeCroy file name format
opts_MM.en_plot = 1;     % enable debugging plots
event_id_ov = 0;

mm_max_y = [];
mm_bgLevel = [];
mm_bgRMS = [];
noise_mm_max_y = [];
int_vector = [];

k=1; 
event_id_prev = -1;

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

    % read files using third party function
    ch_mm = ReadLeCroyBinaryWaveform(ch_mm_str);
        
    if trackerExist == 1
        ch_tr = ReadLeCroyBinaryWaveform(ch_tr_str);
    end

    str_disp=sprintf('Processing file set No. %d', ff);
    disp(str_disp);

    % get number of segments (events) in one file
    nTRCseg = ch_mm.info.nbSegments;
    % get segment length
    lTRCseg = size(ch_mm.y,1);
    % calculate sampling time
    Ts = ch_mm.x(2,1) - ch_mm.x(1,1);

 
    for (i=1:nTRCseg)
        
        if trackerExist == 1
        t_vec_mm=(0:lTRCseg-1)'*Ts + ch_mm.trigger_offset(i);
        
        maxTrackerAmp = max(-ch_tr.y(:,i));
        sumTrackerAmp = sum(-ch_tr.y(:,i));
        if tracker.en %&& maxTrackerAmp>0.8 && sumTrackerAmp>500
                
                % process tracker ID channel
                event_id = process_tr_bitstream(t_vec_mm, ch_tr.y(:,i), opts_TR);
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
        
        
        % Pretrigger region for noise analysis
        %eventXvec = ch_mm.x(1:end-end*0.8,i);
        %eventYvec = -ch_mm.y(1:end-end*0.8,i);
        
        % Full range events
         eventXvec = ch_mm.x(1:end-1,i);
         eventYvec = -ch_mm.y(1:end-1,i);

        
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
            close all
        end
        

        iStart = iMaxY - 200;
        iEnd = iMaxY + 200;

        if (iStart < 1)
            iStart = 1;
        end

        if (iEnd > length(eventXvec))
            iEnd = length(eventXvec);
        end

        % this is the part of the signal that makes the peak
        iEventXvec = eventXvec(iStart:iEnd);
        iEventYvec = eventYvec(iStart:iEnd);


        % Save to data vectors
        int_iEventYvec = sum(iEventYvec); % sum of the amps in the peak. Histogram of the integral of the peak?
        int_vector= [int_vector; int_iEventYvec];
        mm_max_y = [mm_max_y; max_y]; % concatenate onto mm_max_y for each file
        mm_bgLevel = [mm_bgLevel; bg_level]; 
        mm_bgRMS = [mm_bgRMS; bg_rms]; 

        close all

       

            % entire event signal
            noise_eventXvec = ch_mm.x(1:end-1,i);
            noise_eventYvec = -ch_mm.y(1:end-1,i);

            bg_level = mean(noise_eventYvec(1:20));
            noise_eventYvec = noise_eventYvec - bg_level ;

            %don't need max value
            [noise_max_y,noise_iMaxY] = max(noise_eventYvec);

            max_position = noise_eventXvec(iMaxY);

        
            noise_mm_max_y = [noise_mm_max_y; noise_max_y];

            close all

        %end

    end

    end


%AnalyseRun_PEAnalysis_noise
close all

%% ANALYSIS
% histogram and polya fit
format long; % print 15 decimals intead of 4delat

spacial_cut = 1;

% Parameters for changing
shouldOutputAmplitudesTxtFile = true;

useNaiveMeanPos = true; %if median pos extract fails, use simple mean of tracker values instead
    min_v = 0;
    max_v =4e-3;
    % bins
    binOffset = 0;


figureWidth=800;
figureHeight=500;

num_bins = 200;

% Global cut
glbl_cut = mm_max_y>0.00 & mm_max_y<0.95*max(mm_max_y)

figure(1)
x_fig=10;
y_fig=10;
set(gcf,'position',[x_fig,y_fig,figureWidth,figureHeight]);
h=histogram(mm_max_y(glbl_cut),num_bins)
set(gca, 'YScale', 'log') %plot on log scale
movegui(gcf,'south');
hold on
xlabel('Signal amplitude, V');
ylabel('Events');
ylim([0.5 max(h.Values)+1000]);
ax = gca;
ax.FontSize = 20;
if trackerExist == 1
    % xlim([-0.005 0.3]);
else
    %xlim([-0.005 0.1]);
end
grid on
if trackerExist == 1
    title_str = sprintf('PICOSEC beam test - Run %s - Max e-peak amplitude',run.id);
else
    title_str = sprintf('PICOSEC LED test - Run %s - Max e-peak amplitude',run.id);
end
title(title_str)


% set the upper and lower bounds for the lower cutoff of the polya fit that
% we would like to try
min_lower = 3e-3;

% If min and max for fit are not within dataset range pick edges of dataset
if (min(h.BinEdges) > min_lower)
    min_lower = min(h.BinEdges);
end
if (max(h.BinEdges(1:num_bins)) < max_v)
    max_v = max(h.BinEdges(1:num_bins)) % cut the last bin edge because it's on the right side of the bin to make array lengths match
end


% min cut value should be the closest value to min_lower that is greater than

[min_cut_val,min_cut_idx] = min(abs(h.BinEdges-min_v));
if (min_cut_val >= min_v)
else
    min_cut_idx=min_cut_idx+1;
end
[diff,cutUp_idx] = min(abs(h.BinEdges-max_v));
cutUp_val = h.BinEdges(cutUp_idx);
if (cutUp_val <= max_v)
else
    cutUp_idx=cutUp_idx-1;
end

% array of starting voltages is the first n_bins starting with
min_arr = h.BinEdges(min_cut_idx:min_cut_idx+binOffset);

% initialize arrays
mean_arr = [];
chi2_arr = [];

% count exceptions
exceptionCounter = 0;

% initialize error array
mean_err_arr = [];

% colors for plotting
colors = distinguishable_colors(binOffset+1);
c = 1;
%P = ["N","theta","nBar"];
disp("ciao")
for i=1:length(min_arr)
    cut_idx = min_cut_idx+i-1;
   
    fitted_curve = fit((h.BinEdges(1:length(h.Values))+h.BinWidth/2)',h.Values','gauss2')

          %%save fit params
 gaus_mean = [];
 gaus_sigma = [];

 gaus_mean = [gaus_mean; fitted_curve.b1];
 gaus_sigma = [gaus_sigma; fitted_curve.c1];

 %%plot the fitted curve
 plot(h.BinEdges(cut_idx:cutUp_idx),fitted_curve(h.BinEdges(cut_idx:cutUp_idx)),'LineWidth',3,'Color',[colors(c,:)])
 xlim([0 4e-3] );
 
 dlmwrite(['noiseGausseFit.txt'],[fitted_curve.b1; fitted_curve.c1] )

% extract fit errors
cfit_curve = cfit(fitted_curve);
param_errors = confint(cfit_curve);
err_b1 = param_errors(2, 2)-fitted_curve.b1;
err_c1 = param_errors(2, 3)-fitted_curve.c1;


%%save fit

str1 = sprintf('Mean = %0.5f +/- %0.5f\n Sigma = %0.5f +/- %0.5f',fitted_curve.b1,err_b1,fitted_curve.c1,err_c1);
annotation('textbox',[.50 .6 .35 .3], 'String',str1,'FitBoxToText','on','FontSize',18);
saveas(figure(1),[store_folder 'RUN' run.id ' - Channel' channel.id '- E-PeakAmplitude-Noise.png'])

end
%ProcessRawFiles_PEAnalysis_noise_plus_sig

%% PROCESS

close all

%addpath 'C:\Users\GDD\Documents\MATLAB\Picosec\CommonFunctions';
%addpath '..\Functions'

%shouldPlotWaveform = false;
%shouldSaveWaveformSamples = 0;

%numberBGSamples = 100;

%% check if offsets for file numbers exist from oscilloscope saving
%offsetTrackerFileNumber = 0;

%% Tracker setup
%trackerExist = 0 %set 0 for SPE runs

trigIntervalTimes = [];
    
%channel.id = 'C3';
       
%run.name = ['PICOSEC - RUN ' run.id ' - Channel ' channel.id];
run.path = sprintf('\\\\eosproject-smb\\eos\\project\\p\\picosec\\lab\\ResistiveMM-SinglePad\\RMM-SPE-SlowPream-Elec-100624\\%dmVdiv\\%dA-%dC\\signal\\', mVdiv, anodeV, cathodeV);

%     run.nfiles = find_fileNo(run.path);
run.lecroy_name = '--Trace--'; %['Run' run.id];

run.nfiles = 10; % override number of file sets to process (see how many trc files are in the folder)

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

    dutColArray = [1 10 11; 2 13 14; 3 19 20; 4 22 23; 5 4 5; 6 25 26; 7 16 17] %July 2023
   

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

%
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

    ch_mm = ReadLeCroyBinaryWaveform(ch_mm_str);
        
    if trackerExist == 1
        ch_tr = ReadLeCroyBinaryWaveform(ch_tr_str);
    end

    str_disp=sprintf('Processing file set No. %d', ff);
    disp(str_disp);

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
   
        % entire event signal
        eventXvec = ch_mm.x(1:end-1,i);
        eventYvec = -ch_mm.y(1:end-1,i);
        
    
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


%         if ( exist('noiseFile') && (ff >= noiseFile ))
%             debug = sprintf('Noise File No. %d', ff);
%             disp(debug);
% 
%             % entire event signal
%             noise_eventXvec = ch_mm.x(1:end-1,i);
%             noise_eventYvec = -ch_mm.y(1:end-1,i);
% 
%             bg_level = mean(noise_eventYvec(1:20));
%             noise_eventYvec = noise_eventYvec - bg_level ;
% 
%             %don't need max value
%             [noise_max_y,noise_iMaxY] = max(noise_eventYvec);
% 
%             max_position = noise_eventXvec(iMaxY);
% 
%             noise_mm_max_y = [noise_mm_max_y; noise_max_y];
% 
%             close all
%         end
    end
end

meanTimeInterval = mean(trigIntervalTimes);
meanRate = 1/meanTimeInterval

%% ANALYSIS 2
   

close all

% histogram and polya fit
format long; % print 15 decimals intead of 4delat

%geometric cut
radius = 2; % 2 mm radius
spacial_cut = 1;

%% Parameters for changing

shouldOutputAmplitudesTxtFile = true;

useNaiveMeanPos = true; %if median pos extract fails, use simple mean of tracker values instead
gaussFitParameterNoise = dlmread(['noiseGausseFit.txt']);
%gauss fit parameters
b1 = gaussFitParameterNoise(1) %mean
c1 = gaussFitParameterNoise(2); %sigma


% start / end voltages for fitting
if trackerExist == 1
     min_v = 30e-3;
     max_v = 0.22;
    % bins
    binOffset = 0;
else
     min_v = 0.0;
     if mVdiv == 10
        max_v = 0.075;
     elseif mVdiv == 20
        max_V = 0.13;
     else 
        max_V = 0.05;
     end
    % bins
    binOffset = 0;
    
end

% polya starting parameters for fitting   
%x0 = [1 1 0.2 3];   %3 polya parameters (N normalization, theta shape, nBar mean)+noise amplitude
x0 = [1 1 0.2 3 b1 c1]; %[1 1 0.3  b1 c1]; %3 parameters 3LEDs test ortec+3 noise  Deafult: [1 1 0.2 3 b1 c1]

figureWidth=800;
figureHeight=500;

num_bins = 400;

if trackerExist == 1

%sampling area for 2D maps
% calculate
area.step = 0.25;   % set grid resolution in mm
area.size = 8;      % set half-size of the observed square area, mm - for large MPC
%pads and small MCP
area.radius = 1;    % set radius of averaging circular window cut, mm
end


% Global cut
if trackerExist == 1
    glbl_cut = mm_max_y>0.00 & mm_max_y<0.95*max(mm_max_y) & trackerX' ~= 0 & trackerY' ~= 0; %remove saturated datapoints
else
    glbl_cut = mm_max_y>0.00 & mm_max_y<0.95*max(mm_max_y)
end
            
if trackerExist == 1
    %% Pad center
    pad.cut_len = 8;   % move +/-8mm from median to find center
    clear mean
    % calculate pad center
        pad.xc_med = median(trackerX(glbl_cut));
        pad.yc_med = median(trackerY(glbl_cut));
        pad.x_idx_cut = trackerX' > (pad.xc_med-pad.cut_len) & trackerX' < (pad.xc_med+pad.cut_len) & glbl_cut;
        pad.xc_n = mean(trackerX(pad.x_idx_cut)); % naive mean for x
        pad.y_idx_cut = trackerY' > (pad.yc_med-pad.cut_len) & trackerY' < (pad.yc_med+pad.cut_len) & glbl_cut;
        pad.yc_n = mean(trackerY(pad.y_idx_cut)); % naive mean for y

        figure;
        title(['%s - Hits on DUT detector (X, Y projections) ']);
        subplot(2,2,1);
        pad.h_x = histogram(trackerX(pad.x_idx_cut),50);
        subplot(2,2,4);
        pad.h_y = histogram(trackerY(pad.y_idx_cut),50);
        set(gca,'view',[90 -90])
        subplot(2,2,3);
        scatter(trackerX(glbl_cut), trackerY(glbl_cut),'.'); % plot hits with global cut
        axis equal
        xlim(pad.h_x.BinLimits);
        ylim(pad.h_y.BinLimits);
         pause(1);

    % find x center from electron peak mean charge
    for i = 1:length(pad.h_x.Values)
        tmp_cut = trackerX>pad.h_x.BinEdges(i) & trackerX<pad.h_x.BinEdges(i)+ pad.h_x.BinWidth;
        %tmp_cut = trackerX>pad.h_x.BinEdges(i) & trackerX<pad.h_x.BinEdges(i)+ pad.h_x.BinWidth & MM_maxy>0.3*max(MM_maxy);
        pad.epeak_x(i) = mean(mm_max_y(tmp_cut));
    end

    fit_data = [];
    fit_data(1,:) = pad.h_x.BinEdges(1:end-1)+pad.h_x.BinWidth/2;
    fit_data(2,:) = pad.epeak_x;
    % fit_data(3,:) = yerr;
    p0=[];
    p0(1) = pad.xc_n;
    p0(2) = 1;
    p0(3) = 0.01;
    p0(4) = 0.01;
    cmd='min; ret';
    [p, err, chi] = fminuit('parabola4_minuit',p0,fit_data(:,1:end),'-b','-c',cmd);
    % store pad center x
    pad.xc = p(1);
    pad.xc_err = err(1);
    % plot to see how parabolic fit looks like
    figure(10)
    bar(fit_data(1,:),fit_data(2,:));
    hold on
    plot(fit_data(1,:),parabola4_minuit(p, fit_data(1,:)),'LineWidth',2);
    xlabel('x-axis, mm');
    ylabel('Charge, pC');
    legend('RAW', 'Fit')
    title_str = sprintf('%s \n E-peak mean over x-axis');
    title(title_str)
    grid on
  

    % find y center from electron peak mean
    for i = 1:length(pad.h_y.Values)
        tmp_cut = trackerY>pad.h_y.BinEdges(i) & trackerY<pad.h_y.BinEdges(i)+pad.h_y.BinWidth;
        pad.epeak_y(i) = mean(mm_max_y(tmp_cut));
    end

    fit_data = [];
    fit_data(1,:) = pad.h_y.BinEdges(1:end-1)+pad.h_y.BinWidth/2;
    fit_data(2,:) = pad.epeak_y;
    % fit_data(3,:) = yerr;
    p0=[];
    p0(1) = pad.yc_n;
    p0(2) = 1;
    p0(3) = 0.01;
    p0(4) = 0.01;
    cmd='min; ret';
    [p, err, chi] = fminuit('parabola4_minuit',p0,fit_data(:,1:end),'-b','-c',cmd);
    % store pad center y
    pad.yc = p(1);
    pad.yc_err = err(1);
    % plot to see how parabolic fit looks like
    figure(11)
    bar(fit_data(1,:),fit_data(2,:));
    hold on
    plot(fit_data(1,:),parabola4_minuit(p, fit_data(1,:)),'LineWidth',2);
    xlabel('y-axis, mm');
    ylabel('Charge, pC');
    legend('RAW', 'Fit')
    title_str = sprintf('%s \n E-peak mean over y-axis');
    title(title_str)

    grid on
         if useNaiveMeanPos 
             pad.xc = pad.xc_n;
             pad.yc = pad.yc_n;
         end


    %% Plotting
    close all
    spacial_cut = glbl_cut == 1 &  sqrt((trackerX'-pad.xc).^2 + (trackerY'-pad.yc).^2) < radius ;

    noise_glbl_cut = noise_mm_max_y>0.00 & noise_mm_max_y<0.99*max(noise_mm_max_y);

        pad.h_x = histogram(trackerX(pad.x_idx_cut),50);
        pad.h_y = histogram(trackerY(pad.y_idx_cut),50);

        figure(5);
        title(['%s - Hits on DUT detector (X, Y projections) ']);
        scatter(trackerX, trackerY,'.b'); % plot hits with global cut
        hold on;
        scatter(trackerX(spacial_cut), trackerY(spacial_cut),'.r'); % plot hits with global cut
        axis equal
        %xlim(pad.h_x.BinLimits);
        %ylim(pad.h_y.BinLimits);
        pause(1);
        movegui(gcf,'north');

else
    spacial_cut = glbl_cut;
end 



figure(1)
x_fig=10;
y_fig=10;
set(gcf,'position',[x_fig,y_fig,figureWidth,figureHeight]);
h=histogram(mm_max_y(spacial_cut),num_bins)
set(gca, 'YScale', 'log') %plot on log scale
movegui(gcf,'south');
hold on
xlabel('Signal amplitude, V');
ylabel('Events');
ylim([0.5 max(h.Values)+1000]);
ax = gca;
ax.FontSize = 20;

grid on
if trackerExist == 1
    title_str = sprintf('PICOSEC beam test - Run %s - Max e-peak amplitude',run.id);
else
    title_str = sprintf('PICOSEC LED test - Run %s - Max e-peak amplitude',run.id);
end
title(title_str)

% set the upper and lower bounds for the lower cutoff of the polya fit that
% we would like to try
min_lower = 0;
%min_upper = 3.5e-3;

% If min and max for fit are not within dataset range pick edges of dataset
if (min(h.BinEdges) > min_lower)
    min_lower = min(h.BinEdges);
end
if (max(h.BinEdges(1:num_bins)) < max_v)
    max_v = max(h.BinEdges(1:num_bins)) % cut the last bin edge because it's on the right side of the bin to make array lengths match
end
%max_v = max(h.BinEdges(1:num_bins));


% min cut value should be the closest value to min_lower that is greater than

[min_cut_val,min_cut_idx] = min(abs(h.BinEdges-min_v));
if (min_cut_val >= min_v)
else
    min_cut_idx=min_cut_idx+1;
end
[diff,cutUp_idx] = min(abs(h.BinEdges-max_v));
cutUp_val = h.BinEdges(cutUp_idx);
if (cutUp_val <= max_v)
else
    cutUp_idx=cutUp_idx-1;
end

% array of starting voltages is the first n_bins starting with
min_arr = h.BinEdges(min_cut_idx:min_cut_idx+binOffset);

% initialize arrays
mean_arr = [];
chi2_arr = [];

% take upper cut as very end of array
%cutUp_idx = length(h.Values);

% count exceptions
exceptionCounter = 0;

% initialize error array
mean_err_arr = [];

% colors for plotting
colors = distinguishable_colors(binOffset+1);
c = 1;
%P = ["N","theta","nBar", "a1"];
P = ["N","theta","nBar", "a1", "b1", "c1"];

% Indices of the bins we want to exlude
exludeBins = [15:20];

% For polya+gaus fit, one is with fixed mean and sigma of the gaus and one with all params free
for i=1:length(min_arr)
    cut_idx = min_cut_idx+i-1; 
    % define polya+gauss
    %fitfun = fittype( @(N,theta,nBar,a1,x) (N./nBar).*((theta+1).^(theta+1)).*((x./nBar).^theta).*exp(-(theta+1).*x./nBar)./gamma(theta+1) + a1*exp(-((x-b1)/c1).^2));
    fitfun = fittype( @(N,theta,nBar,a1,b1,c1,x) (N./nBar).*((theta+1).^(theta+1)).*((x./nBar).^theta).*exp(-(theta+1).*x./nBar)./gamma(theta+1) + a1*exp(-((x-b1)/c1).^2));

    options = fitoptions(fitfun);
    %options.Upper = [300 30 1 10000];
    options.Upper = [1000 150 0.05 100000 b1+b1*0.5 c1+c1*2.5];
    %options.Lower = [0 0 0 0];
    options.Lower = [0 0 0 0 0 0];

   % x0 = [2.7 0.7 0.03 4700]; 
    %options.StartPoint = x0;
  
    try
        [fitted_curve,gof] = fit((h.BinEdges(cut_idx:cutUp_idx)+h.BinWidth/2)',h.Values(cut_idx:cutUp_idx)',fitfun,options)
    catch exception
        exceptionCounter = exceptionCounter+1;
        mean_arr(i) = 0;
        chi2_arr(i) = 0;
        mean_err_arr(i,:) = [0;0];
        continue
    end

    %P = [P;fitted_curve.N,fitted_curve.theta,fitted_curve.nBar,fitted_curve.a1];
    P = [P;fitted_curve.N,fitted_curve.theta,fitted_curve.nBar,fitted_curve.a1,fitted_curve.b1,fitted_curve.c1];


    %plot the fitted curve
   if (fitted_curve.theta > 0)
       hold on
       plot(h.BinEdges(cut_idx:cutUp_idx),fitted_curve(h.BinEdges(cut_idx:cutUp_idx)),'LineWidth',3,'Color',[colors(c,:)])
       % used color, increment color counter
       %saveas(figure(1),[store_folder 'RUN' run.id ' - Channel' channel.id '- E-PeakAmplitude.png'])
        
       %calculate chi squared and degrees of freedom
       ch2 = sum(((h.Values(cut_idx:cutUp_idx)-fitted_curve(h.BinEdges(cut_idx:cutUp_idx))').^2)./fitted_curve(h.BinEdges(cut_idx:cutUp_idx))');
       dof = size(h.Values(cut_idx:cutUp_idx),2)-3;
       nch2 = ch2/dof
       np = 1-chi2cdf(ch2,dof) % P(\chi^2>ch2)
        
       mean_arr(i) = fitted_curve.nBar;
       chi2_arr(i) = nch2;
       mean_err = confint(cfit(fitted_curve),0.68);
       mean_err_arr(i,:) = mean_err(:,3);
        
       c = c+1;
   else
       exceptionCounter = exceptionCounter+1;
       mean_arr(i) = 0;
       chi2_arr(i) = 0;
       mean_err_arr(i,:) = [0;0];
   end
 
end


%% calculate the final mean from all means
format shortE

% remove 0s which correspond to fits that threw exceptions, do not want 0s
% biasing the mean
filter = mean_arr ~= 0;
mean_arr = mean_arr(filter);
chi2_arr = chi2_arr(filter);
mean_err_arr_temp = mean_err_arr(filter,:);
mean_err_arr = mean_err_arr_temp;
min_arr = min_arr(filter);
meanCalculated = sum(mean_arr)/length(mean_arr);

%error propagation
delta_arr(1) = (1/2)*(sum((mean_err_arr(:,1)'-mean_arr).^2))^(1/2); % lower
delta_arr(2) = (1/2)*(sum((mean_err_arr(:,2)'-mean_arr).^2))^(1/2); % upper

% chi2_arr

str1 = sprintf('Mean = %0.5f +/- %0.5f',meanCalculated,delta_arr(2));
annotation('textbox',[.50 .6 .35 .3], 'String',str1,'FitBoxToText','on','FontSize',18);
saveas(figure(1),[store_folder 'RUN' run.id ' - Channel' channel.id '- E-PeakAmplitude-Signal.png'])
% saveas(figure(1),[store_folder ' Channel' channel.id '- E-PeakAmplitude.png'])


%% plot the chi2, mean and error in the mean

hold off;
figure(2)
x_fig=10;
y_fig=10;
set(gcf,'position',[x_fig,y_fig,figureWidth,figureHeight]);
movegui(gcf,'northwest');

%plot(min_arr,chi2_arr, 'Marker','o','Color',[colors(:)]);
scatter(min_arr,chi2_arr,35,[colors(1:c-1,:)],'filled');
title('\chi^2/DOF vs Polya Mean');
xlabel('nBar [V]');
ylabel('\chi^2/DOF');
ax = gca;
ax.FontSize = 20;
grid on;
saveas(figure(2),[store_folder 'RUN' run.id ' - Channel' channel.id '- MeanVsChi2.png'])
% saveas(figure(2),[store_folder ' Channel' channel.id '- MeanVsChi2.png'])


figure(3)
x_fig=10;
y_fig=10;
set(gcf,'position',[x_fig,y_fig,figureWidth,figureHeight]);
%movegui(gcf,'e');

hold on;

for i = 1:length(min_arr)
    e = errorbar(min_arr(i),mean_arr(i),abs(mean_arr(i)-mean_err_arr(i,1)'),abs(mean_arr(i)-mean_err_arr(i,2)'),'Marker','o','Color',[colors(i,:)]);
    e.LineWidth = 1
end

title('Polya Mean vs Lower Cutoff Voltage');
ylabel('nBar [V]');
xlabel('Lower Polya Cutoff [V]');
ax = gca;
ax.FontSize = 20;
grid on;
saveas(figure(3),[store_folder 'RUN' run.id ' - Channel' channel.id '- MeanVsCutoffError.png'])
% saveas(figure(3),[store_folder ' Channel' channel.id '- MeanVsCutoffError.png'])
movegui(gcf,'northwest');

hold off;
