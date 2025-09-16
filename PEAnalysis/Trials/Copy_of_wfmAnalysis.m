clear all
close all

addpath 'C:\Users\GDD\Documents\MATLAB\Picosec\CommonFunctions';
addpath '..\Functions'

shouldPlotWaveform = false;
shouldSaveWaveformSamples = 0;

numberBGSamples = 10;

% Check if offsets for file numbers exist from oscilloscope saving
offsetTrackerFileNumber = 0;

%% DATA PROCESSING
trackerExist = 0; % 0 for SPE runs and lab

% Variable to include or not the gaussians on the analysis of the signal
only_polya = false;
% Define run values
cathodeV = 500;
anodeV = 275;
mVdiv = 10;
folder = 'RMM-SPE-SlowPream-CsI-STD-240624';
% 500 a 10mVdiv, ...
% Specify run info
run.id = sprintf('%dV', cathodeV);
run.year = '2024Lab';   
channel.id = 'C3';
run.name = ['PICOSEC - RUN ' run.id ' - Channel ' channel.id];
% run.path = sprintf('\\\\eosproject-smb\\eos\\project\\p\\picosec\\lab\\ResistiveMM-SinglePad\\%s\\%dmVdiv\\%dA-%dC\\signal\\', folder, mVdiv, anodeV, cathodeV);
run.path = sprintf('\\\\eosproject-smb\\eos\\project\\p\\picosec\\lab\\ResistiveMM-SinglePad\\%s\\%dmVdiv\\%dA-%dC\\', folder, mVdiv, anodeV, cathodeV);
run.lecroy_name = '--Trace--';
run.nfiles = 15; % Override number of file sets to process (see how many trc files are in the folder)
store_folder = sprintf('\\\\eosproject-smb\\eos\\project\\p\\picosec\\lab\\ResistiveMM-SinglePad\\%s\\%dmVdiv\\%dA-%dC\\Result\\', folder, mVdiv, anodeV, cathodeV);
mkdir(store_folder);

str_disp=sprintf('- - - - - ANALYSING RUN %s - - - - - ', run.id);
disp(str_disp);

trigIntervalTimes = [];

shiftX = 0;
shiftY = 0;

if trackerExist == 1
    tracker.path = ['\\eosproject-smb\eos\project\p\picosec\testbeam\2023_July_h4\tracker\reconstructed\asciiRun' run.id '.dat'];
    tracker.en = 1; % Match eventIDs to tracking data and add XY to output
    
    tracker.dutIndex = 2; % DUT index: 1: MM1 (Multipad), 2: MM2(VacChamber), 3: MM3 (ElectronsetupMM), 4: MM4 (Picolarge main tracker), 5:MM5 (SiPM)

    dutColArray = [1 10 11; 2 16 17; 3 4 5; 4 13 14; 5 7 8]; %[dID colXID colYID; ] -> July 2022 MM run

    if tracker.en
        trackerFile = fopen(tracker.path,'rt');
        D = textscan(trackerFile, '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f', 'Delimiter','\t', 'HeaderLines',2, 'CollectOutput',1);
        tracker.data = cell2mat(D);
    end

    % Options for extracting tracker ID from the bitstream
    opts_TR.ch_name = ['C3' run.lecroy_name];
    opts_TR.baud_rate = 40e6; % 40Mbps baudrate
    opts_TR.n_bits = 16;      % Bits after start bit

end

opts_MM.ch_name = [channel.id run.lecroy_name]; % LeCroy file name format
opts_MM.en_plot = 1;     % Enable debugging plots
event_id_ov = 0;

% Vectors to store data of pretrigger region and full wfm (full_...)
mm_max_y = [];          full_mm_max_y = [];
mm_bgLevel = [];        full_mm_bgLevel = [];
mm_bgRMS = [];          full_mm_bgRMS = [];
noise_mm_max_y = [];    full_noise_mm_max_y = [];
int_vector = [];        full_int_vector = [];

k=1; 
event_id_prev = -1;

for (ff=0:run.nfiles)
    % Generate file name strings
    ch_mm_str=sprintf('%s%s%05d.trc', run.path, opts_MM.ch_name,ff);
    if trackerExist == 1
        ch_tr_str=sprintf('%s%s%05d.trc', run.path, opts_TR.ch_name,ff+offsetTrackerFileNumber); 
    end

    % Check if the file number exists, if not, go to next number
    if not(isfile(ch_mm_str))
        continue
    end

    % Display current file number
    str_disp=sprintf('Loading file set No. %d', ff);
    disp(str_disp);

    % Read files using third party function
    ch_mm = ReadLeCroyBinaryWaveform(ch_mm_str);
        
    if trackerExist == 1
        ch_tr = ReadLeCroyBinaryWaveform(ch_tr_str);
    end

    str_disp=sprintf('Processing file set No. %d', ff);
    disp(str_disp);

    % Retrieve segments (events) in one file
    nTRCseg = ch_mm.info.nbSegments;
    lTRCseg = size(ch_mm.y,1); % Event lenght
    Ts = ch_mm.x(2,1) - ch_mm.x(1,1); % Sampling time

 
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
 
        if tracker.en
            % Process tracker ID channel
            event_id = process_tr_bitstream(t_vec_mm, ch_tr.y(:,i), opts_TR);
            
            % Count overflows
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

            % Find corresponding tracker entry and save XY info
            trackIdx = find(tracker.data(:,1) == MM_temp.event_id);
            if trackIdx>0
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
    
        % Store valid data into array of structures
        MM_data(k)= MM_temp;
        trackerX(k) = MM_temp.x;
        trackerY(k) = MM_temp.y;
        k=k+1;
    end
        
        
    % Pretrigger region for noise analysis
    eventXvec = ch_mm.x(1:end-end*0.8,i);
    eventYvec = -ch_mm.y(1:end-end*0.8,i);    
   
    % Full range events
    full_eventXvec = ch_mm.x(1:end-1,i);
    full_eventYvec = -ch_mm.y(1:end-1,i);
   
    % Pretrigger region
    bg_level = mean(eventYvec(1:numberBGSamples));
    eventYvec = eventYvec - bg_level ;
    bg_levelCorrected = mean(eventYvec(1:numberBGSamples));
    bg_rms = std(eventYvec(1:numberBGSamples));
    
    % Full range
    full_bg_level = mean(full_eventYvec(1:numberBGSamples));
    full_eventYvec = full_eventYvec - full_bg_level ;
    full_bg_levelCorrected = mean(full_eventYvec(1:numberBGSamples));
    full_bg_rms = std(full_eventYvec(1:numberBGSamples));
  
    % Pretrigger region
    [max_y,iMaxY] = max(eventYvec);
    max_position = eventXvec(iMaxY);
  
    % Full range
    [full_max_y,full_iMaxY] = max(full_eventYvec);
    full_max_position = full_eventXvec(full_iMaxY);
    
    if shouldSaveWaveformSamples >0 && length(mm_bgLevel)<shouldSaveWaveformSamples
        plot(eventXvec-min(eventXvec), eventYvec); 
        ylim([-0.005 0.1] );
        pause(1);
        close all
    end
    
    % Pretrigger region
    iStart = iMaxY - 200;
    iEnd = iMaxY + 200;
   
    % Full range
    full_iStart = full_iMaxY - 200;
    full_iEnd = full_iMaxY + 200;

    if (iStart < 1)
        iStart = 1;
    end

    if (iEnd > length(eventXvec))
        iEnd = length(eventXvec);
    end

    if (full_iStart < 1)
        full_iStart = 1;
    end

    if (full_iEnd > length(full_eventXvec))
        full_iEnd = length(full_eventXvec);
    end

    % Define peak region
    iEventXvec = eventXvec(iStart:iEnd);
    iEventYvec = eventYvec(iStart:iEnd);
    full_iEventXvec = full_eventXvec(full_iStart:full_iEnd);
    full_iEventYvec = full_eventYvec(full_iStart:full_iEnd);

    % Save to data vectors
    int_iEventYvec = sum(iEventYvec); % sum of the amps in the peak. Histogram of the integral of the peak?
    int_vector= [int_vector; int_iEventYvec];
    mm_max_y = [mm_max_y; max_y]; % concatenate onto mm_max_y for each file
    mm_bgLevel = [mm_bgLevel; bg_level]; 
    mm_bgRMS = [mm_bgRMS; bg_rms]; 
    
    full_int_iEventYvec = sum(full_iEventYvec); % sum of the amps in the peak. Histogram of the integral of the peak?
    full_int_vector= [full_int_vector; full_int_iEventYvec];
    full_mm_max_y = [full_mm_max_y; full_max_y]; % concatenate onto mm_max_y for each file
    full_mm_bgLevel = [full_mm_bgLevel; full_bg_level]; 
    full_mm_bgRMS = [full_mm_bgRMS; full_bg_rms]; 

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
    
    % entire event signal
    full_noise_eventXvec = ch_mm.x(1:end-1,i);
    full_noise_eventYvec = -ch_mm.y(1:end-1,i);

    full_bg_level = mean(full_noise_eventYvec(1:20));
    full_noise_eventYvec = full_noise_eventYvec - full_bg_level ;

    %don't need max value
    [full_noise_max_y,full_noise_iMaxY] = max(full_noise_eventYvec);

    full_max_position = full_noise_eventXvec(iMaxY);
    full_noise_mm_max_y = [noise_mm_max_y; full_noise_max_y];

    close all

    end

end

meanTimeInterval = mean(trigIntervalTimes);
meanRate = 1/meanTimeInterval;

close all


%% NOISE ANALYSIS

% Hist and polya fit
format long;
spacial_cut = 1;

% Parameters for changing
shouldOutputAmplitudesTxtFile = true;
useNaiveMeanPos = true; % If median pos extract fails, use simple mean of tracker values instead
min_v = 0;
max_v = 4e-3;
binOffset = 0;
figureWidth=800;
figureHeight=500;
num_bins = 200;

% Global cut
glbl_cut = mm_max_y>0.00 & mm_max_y<0.95*max(mm_max_y);

% Define figure params
figure(1);
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
grid on
if trackerExist == 1
    title_str = sprintf('PICOSEC beam test - Run %s - Max e-peak amplitude',run.id);
else
    title_str = sprintf('PICOSEC LED test - Run %s - Max e-peak amplitude',run.id);
end
title(title_str)


% Set the upper and lower bounds for the lower cutoff of the polya fit to be tried
min_lower = 3e-3;

% If min and max for fit are not within dataset range pick edges of dataset
if (min(h.BinEdges) > min_lower)
    min_lower = min(h.BinEdges);
end
if (max(h.BinEdges(1:num_bins)) < max_v)
    max_v = max(h.BinEdges(1:num_bins)) % cut the last bin edge because it's on the right side of the bin to make array lengths match
end


% Min cut value should be the closest value to min_lower that is greater than
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

% Starting voltages is the first n_bins starting with
min_arr = h.BinEdges(min_cut_idx:min_cut_idx+binOffset);

% Arrays to save the data
mean_arr = [];
chi2_arr = [];

% Variable to count exceptions
exceptionCounter = 0;

% Initialize error array
mean_err_arr = [];

% Colors for plotting
colors = distinguishable_colors(binOffset+1);
c = 1;

for i=1:length(min_arr)
    cut_idx = min_cut_idx+i-1;
   
    fitted_curve = fit((h.BinEdges(1:length(h.Values))+h.BinWidth/2)',h.Values','gauss2')

    % Store fit params
    gaus_mean = [];
    gaus_sigma = [];    
    gaus_mean = [gaus_mean; fitted_curve.b1];
    gaus_sigma = [gaus_sigma; fitted_curve.c1];
    
    % Plot the fitted curve
    plot(h.BinEdges(cut_idx:cutUp_idx),fitted_curve(h.BinEdges(cut_idx:cutUp_idx)),'LineWidth',3,'Color',[colors(c,:)])
    xlim([0 4e-3]);
    dlmwrite(['noiseGausseFit.txt'],[fitted_curve.b1; fitted_curve.c1] )

    % Fit errors
    cfit_curve = cfit(fitted_curve);
    param_errors = confint(cfit_curve);
    err_b1 = param_errors(2, 2)-fitted_curve.b1;
    err_c1 = param_errors(2, 3)-fitted_curve.c1;
    
    % Save fit
    str1 = sprintf('Mean = %0.5f +/- %0.5f\n Sigma = %0.5f +/- %0.5f',fitted_curve.b1,err_b1,fitted_curve.c1,err_c1);
    annotation('textbox',[.50 .6 .35 .3], 'String',str1,'FitBoxToText','on','FontSize',18);
    saveas(figure(1),[store_folder 'RUN' run.id ' - Channel' channel.id '- E-PeakAmplitude-Noise.png'])

end

%% SIGNAL ANALYSIS 
   
close all

% Hist and polya fit
format long; % print 15 decimals intead of 4 delat

%geometric cut
radius = 2; % 2 mm radius
spacial_cut = 1;

% Parameters for changing

shouldOutputAmplitudesTxtFile = true;
useNaiveMeanPos = true; %if median pos extract fails, use simple mean of tracker values instead
gaussFitParameterNoise = dlmread(['noiseGausseFit.txt']);

% Define gauss fit params
b1 = gaussFitParameterNoise(1); % Mean
c1 = gaussFitParameterNoise(2); % Sigma


% start / end voltages for fitting
if trackerExist == 1
     min_v = 30e-3;
     max_v = 0.22;
     binOffset = 0;
else
     if only_polya == false 
        min_v = 0.0;
     else 
        min_v = 0.001;
     end
     if mVdiv == 10
        max_v = 0.075;
     elseif mVdiv == 20
        max_v = 0.15;
     elseif mVdiv == 50
        max_v = 0.5;
     else 
        max_v = 0.05;
     end
     binOffset = 0;
    
end

% Polya initial params  
x0 = [1 1 0.2 3 b1 c1]; %3 parameters 3LEDs test ortec+3 noise

figureWidth=800;
figureHeight=500;
num_bins = 400;

if trackerExist == 1    
    % Sample area for 2D maps
    area.step = 0.25;   % Set grid resolution in mm
    area.size = 8;      % Set half-size of the observed square area, mm - for large MPC
    % For pads and small MCP
    area.radius = 1;    % Set radius of averaging circular window cut, mm
end


% Global cut
if trackerExist == 1
    glbl_cut = full_mm_max_y>0.00 & full_mm_max_y<0.95*max(full_mm_max_y) & trackerX' ~= 0 & trackerY' ~= 0; %remove saturated datapoints
else
    glbl_cut = full_mm_max_y>0.00 & full_mm_max_y<0.95*max(full_mm_max_y)
end
            
if trackerExist == 1
    % Pad center
    pad.cut_len = 8;   % Move +/-8mm from median to find center
    clear mean
    % Calculate pad center
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

    % Find x center from electron peak mean charge
    for i = 1:length(pad.h_x.Values)
        tmp_cut = trackerX>pad.h_x.BinEdges(i) & trackerX<pad.h_x.BinEdges(i)+ pad.h_x.BinWidth;
        pad.epeak_x(i) = mean(full_mm_max_y(tmp_cut));
    end

    fit_data = [];
    fit_data(1,:) = pad.h_x.BinEdges(1:end-1)+pad.h_x.BinWidth/2;
    fit_data(2,:) = pad.epeak_x;
    p0=[];
    p0(1) = pad.xc_n;
    p0(2) = 1;
    p0(3) = 0.03;
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
  

    % Find y center from electron peak mean
    for i = 1:length(pad.h_y.Values)
        tmp_cut = trackerY>pad.h_y.BinEdges(i) & trackerY<pad.h_y.BinEdges(i)+pad.h_y.BinWidth;
        pad.epeak_y(i) = mean(full_mm_max_y(tmp_cut));
    end

    fit_data = [];
    fit_data(1,:) = pad.h_y.BinEdges(1:end-1)+pad.h_y.BinWidth/2;
    fit_data(2,:) = pad.epeak_y;

    p0=[];
    p0(1) = pad.yc_n;
    p0(2) = 1;
    p0(3) = 0.01;
    p0(4) = 0.01;
    cmd='min; ret';
    [p, err, chi] = fminuit('parabola4_minuit',p0,fit_data(:,1:end),'-b','-c',cmd);
    
    % Store pad center y
    pad.yc = p(1);
    pad.yc_err = err(1);

    % Plot parabolic fit
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

    % Plot
    close all
    spacial_cut = glbl_cut == 1 &  sqrt((trackerX'-pad.xc).^2 + (trackerY'-pad.yc).^2) < radius ;

    noise_glbl_cut = full_noise_mm_max_y>0.00 & full_noise_mm_max_y<0.99*max(full_noise_mm_max_y);

    pad.h_x = histogram(trackerX(pad.x_idx_cut),50);
    pad.h_y = histogram(trackerY(pad.y_idx_cut),50);

    figure(5);
    title(['%s - Hits on DUT detector (X, Y projections) ']);
    scatter(trackerX, trackerY,'.b'); % Plot hits with global cut
    hold on;
    scatter(trackerX(spacial_cut), trackerY(spacial_cut),'.r'); % plot hits with global cut
    axis equal
    pause(1);
    movegui(gcf,'north');

else
    spacial_cut = glbl_cut;
end 

figure(1)
x_fig=10;
y_fig=10;
set(gcf,'position',[x_fig,y_fig,figureWidth,figureHeight]);
h=histogram(full_mm_max_y(spacial_cut),num_bins)
% h.FaceColor = 'red'; h.EdgeColor = 'black';
set(gca, 'YScale', 'log') %plot on log scale
movegui(gcf,'south');
hold on
xlabel('Signal amplitude [V]');
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

% colors for plotting
colors = distinguishable_colors(binOffset+1);
c = 1;
if only_polya == false
    P = ["N","theta","nBar", "a1", "b1", "c1"];
else
    P = ["N","theta","nBar"];
end


% Indices of the bins we want to exclude
exludeBins = [15:20];

% For polya+gaus fit, one is with fixed mean and sigma of the gaus and one with all params free
for i=1:length(min_arr)
    cut_idx = min_cut_idx+i-1; 
    
    if only_polya == false
    % Define polya+gauss
        fitfun = fittype( @(N,theta,nBar,a1,b1,c1,x) (N./nBar).*((theta+1).^(theta+1)).*((x./nBar).^theta).*exp(-(theta+1).*x./nBar)./gamma(theta+1) + a1*exp(-((x-b1)/c1).^2));
        options = fitoptions(fitfun);
        options.Upper = [1000 150 0.05 100000 b1+b1*0.5 c1+c1*2.5];
        options.Lower = [0 0 0 0 0 0];

    else
        fitfun = fittype( @(N,theta,nBar,x) (N./nBar).*((theta+1).^(theta+1)).*((x./nBar).^theta).*exp(-(theta+1).*x./nBar)./gamma(theta+1));
        options = fitoptions(fitfun);
        options.Upper = [1000 150 0.05];
        options.Lower = [0 0 0 ];
    end

    try
        [fitted_curve,gof] = fit((h.BinEdges(cut_idx:cutUp_idx)+h.BinWidth/2)',h.Values(cut_idx:cutUp_idx)',fitfun,options)
    catch exception
        exceptionCounter = exceptionCounter+1;
        mean_arr(i) = 0;
        chi2_arr(i) = 0;
        mean_err_arr(i,:) = [0;0];
        continue
    end
    
    if only_polya == false
        P = [P;fitted_curve.N,fitted_curve.theta,fitted_curve.nBar,fitted_curve.a1,fitted_curve.b1,fitted_curve.c1];
    else 
        P = [P;fitted_curve.N,fitted_curve.theta,fitted_curve.nBar];
    end
    % Plot fit
    if (fitted_curve.theta > 0)
       hold on;
       %%%%%%%%%%%%%%%%%%%%%%
       plot(h.BinEdges(cut_idx:cutUp_idx),fitted_curve(h.BinEdges(cut_idx:cutUp_idx)),'LineWidth',3,'Color',[colors(c,:)]);
        
       % Calculate chi-squared and degrees of freedom
       ch2 = sum(((h.Values(cut_idx:cutUp_idx)-fitted_curve(h.BinEdges(cut_idx:cutUp_idx))').^2)./fitted_curve(h.BinEdges(cut_idx:cutUp_idx))');
       dof = size(h.Values(cut_idx:cutUp_idx),2)-3;
       nch2 = ch2/dof;
       np = 1-chi2cdf(ch2,dof); % P(\chi^2>ch2)
        
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


% Calculate the final mean from all means
format shortE

% Remove 0s associated to exceptions (to not bias the mean)
filter = mean_arr ~= 0;
mean_arr = mean_arr(filter);
chi2_arr = chi2_arr(filter);
mean_err_arr_temp = mean_err_arr(filter,:);
mean_err_arr = mean_err_arr_temp;
min_arr = min_arr(filter);
meanCalculated = sum(mean_arr)/length(mean_arr);

% Error propagation
delta_arr(1) = (1/2)*(sum((mean_err_arr(:,1)'-mean_arr).^2))^(1/2); % lower
delta_arr(2) = (1/2)*(sum((mean_err_arr(:,2)'-mean_arr).^2))^(1/2); % upper

chi2_arr

str1 = sprintf('Mean = %0.5f +/- %0.5f',meanCalculated,delta_arr(2));
annotation('textbox',[.50 .6 .35 .3], 'String',str1,'FitBoxToText','on','FontSize',18);
saveas(figure(1),[store_folder 'RUN' run.id ' - Channel' channel.id '- E-PeakAmplitude-Signal.png']);

hold off;
