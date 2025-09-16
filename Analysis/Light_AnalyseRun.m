%clear all
close all

%load input file
%in_file = 'C:\Users\GDD\Documents\Picosec\May22\Analysed\Run0123-GDD.mat';

shouldSaveMat = false;

shouldUseCoarseMedianEstimate = true; %prefilter data with max counts for median time diff estiamation

%%configuration
twalk.en = 1;                                   % enable timewalk correction
enable_cut_twalk = true;           %plot the timewalk correction in the sampling area

% center resolution circular - radius
small.r = 2; %used for only circle
%small.r = 2.5; %used for only circle - 15mm MM + 5mm MgF2
medium.r = 2.5; %used for ring only

%shift sampling point for time res relative to pad center determined by
%alignmnet
%account for misalignment
shiftX = 0;
shiftY = 0;

alignToDUT = true; %if false, alignment to REF MCP

%use fit of amplitude to center sampling area
shouldCenterPadWithFit = false;

%use with large trigger to calculate correct arrival time median for global
%cut - for larger trigger or scanning runs
shouldDetermineSATMedianWithGeoCut = false;

%pad size (visualisation + efficiency calc)
pad.size = 15;     % MCP size in mm - Hamamatsu
%pad.size = 15;     % MM sensor size
pad.cut_len = 16;   % move +/-8mm from median to find center
pad.isCircle = 0;   %plotting round if circle, otherwise square

%July22: DUT index: 1: MM1 (Multipad), 2: MM2(VacChamber), 3: MM3
%(ElectronsetupMM), 4: MM4 (Picolarge main tracker), 5:MM5 (SiPM)
if tracker.dutIndex == 1
    pad.isCircle = 0;
    %pad.isCircle = 1;
elseif tracker.dutIndex == 2
    pad.isCircle = 1;
elseif tracker.dutIndex == 3
    pad.isCircle = 1;
elseif tracker.dutIndex == 4
    pad.isCircle = 1;
elseif tracker.dutIndex == 5
    pad.isCircle = 1;
elseif tracker.dutIndex == 6
    pad.isCircle = 1;
elseif tracker.dutIndex == 7
    pad.isCircle = 0;
end


%sampling area for 2D maps
%% calculate
area.step = 0.25;   % set grid resolution in mm
area.size = 8;      % set half-size of the observed square area, mm - for
%pads and small MCP
%area.size = 20;      % set half-size of the observed square area, mm - for large MCP
area.radius = 1;    % set radius of averaging circular window cut, mm


%load(in_file);

runID = str2num(run.id);
runTitleString = [run.name ' ' run.oscilloscope ' DUT:C' opts_MM.chID ' REF:C' opts_MCP.chID ' ' runInfoString ' '];

dataFilePath = append('/eos/project/p/picosec/testbeam/2025_July_h4/info/MatlabRunInfo');
runDataFileName=append('RunData',string(runID),'.mat');

%finish config
pad.curvature = [1,1]; %circle
if pad.isCircle == 0
    pad.curvature = [0,0]; %rect
end


%% define folder for storing data
store_folder = ['\\eosproject-smb\eos\project\p\picosec\testbeam\2025_July_h4\Results\Run' run.id '-' run.oscilloscope '-dut' opts_MM.chID '-ref' opts_MCP.chID '-' runInfoString];
mkdir(store_folder);

resultsTablePath = ['\\eosproject-smb\eos\project\p\picosec\testbeam\2025_July_h4\Results\resultsTable.txt'];

%% finished setting up, select events to include

%%used only for electron run - match with scinitllators on other scope
%run through all files and implement new array
%scintillatorAccepted =  zeros(length(time_MCP),1);

%% do some initial resolution measurement
% make simple cut with respect to median time and minimum amplitude magic
% numbers :-)
%time_avg_raw = median(time_diff_sigmoid(MM_maxy>0.01));  % assume mu from the median vaule
%time_avg_raw = median(time_diff_sigmoid(MM_maxy>0.0001));  % assume mu from the median vaule

[counts,centers] = hist(time_diff_sigmoid,1000);
[maxCounts,maxIdx] = max(counts);
coarseMedianEstimate = centers(maxIdx);
time_minCoarse = coarseMedianEstimate - 0.5;     % predicted resolution 100ps cut 3 sigma
time_maxCoarse = coarseMedianEstimate + 0.5;     % left and 3 sigma right from median

%override time_minCoarse and time_maxCoarse
% time_minCoarse = -8.5;
% time_maxCoarse = -7;

coarseCut = time_diff_sigmoid>time_minCoarse & time_diff_sigmoid<time_maxCoarse & MM_maxy>0.0001;
time_diff_sigmoidCoarseFiltered = time_diff_sigmoid(coarseCut);



if shouldDetermineSATMedianWithGeoCut
    AnalyseRun_DetermineSATMedian
else
    
    if shouldUseCoarseMedianEstimate
        time_avg_raw = median(time_diff_sigmoidCoarseFiltered);  % assume mu from the median vaule
    else
        time_avg_raw = median(time_diff_sigmoid(MM_maxy>0.0001));  % assume mu from the median vaule
    end
end
%time_avg_raw = 12.5;
%time_avg_raw
time_min = time_avg_raw - 0.3;     % predicted resolution 100ps cut 3 sigma
time_max = time_avg_raw + 0.3;     % left and 3 sigma right from median

%time_avg_raw = median(time_diff_sigmoid);  % assume mu from the median vaule
  %time_min = time_avg_raw - 3;     % predicted resolution 100ps cut 3 sigma
  %time_max = time_avg_raw + 3;     % left and 3 sigma right from median

% find events within the cut and cut with respect to amplirudes and
% existance of the tracker data (GLOBAL CUT)
%with tracker

 glbl_cut = time_diff_sigmoid>time_min & time_diff_sigmoid<time_max & ...
     MCP_maxy>0.01*max(MCP_maxy) & MM_maxy>0.01*max(MM_maxy)&...
     MM_maxy<0.99*max(MM_maxy)& MCP_maxy<0.95*max(MCP_maxy) & ...
     trackerX~=0 & trackerY~=0 & time_diff_sigmoid>-200;



  
%glbl_cut = MCP_maxy>0.01*max(MCP_maxy) & MM_maxy>0.01*max(MM_maxy)&...
 %   MM_maxy<0.99*max(MM_maxy)& MCP_maxy<0.95*max(MCP_maxy) & ...
  % trackerX~=0 & trackerY~=0 ;

%additional cuts
glbl_cut = glbl_cut & MCP_maxy>0.15;

%notracker, with S1+S2
%glbl_cut = time_diff_sigmoid>time_min & time_diff_sigmoid<time_max & MCP_maxy>0.01*max(MCP_maxy) & MM_maxy>0.01*max(MM_maxy)& MM_maxy<0.95*max(MM_maxy)& MCP_maxy<0.95*max(MCP_maxy) & s1Signal>0.1 & s2Signal>0.1;

%only observeMM
%glbl_cut =  trackerX~=0 & trackerY~=0;

%without tracker
%glbl_cut = time_diff_sigmoid>time_min & time_diff_sigmoid<time_max & MCP_maxy>0.01*max(MCP_maxy) & MM_maxy>0.01*max(MM_maxy)& MM_maxy<0.95*max(MM_maxy)& MCP_maxy<0.95*max(MCP_maxy);

% make cut on time difference
time_diff_cut = time_diff_sigmoid(glbl_cut);
mean(time_diff_cut)
std(time_diff_cut)

% Time walk analysis (this needs to be improoved)
twalk.n_epk = 50;                               % number of e-peak bins
e_peak_srt = sort(e_peak_MM(glbl_cut));         % make temporary sort of e-peaks
% try to distrubute e_peak vector evenly over the amplitude range
twalk.epeak_vec = e_peak_srt(1:round(length(e_peak_srt)/twalk.n_epk):end); % 50 bins
for i=1:length(twalk.epeak_vec)-1
    temp_cut = e_peak_MM > twalk.epeak_vec(i) & e_peak_MM < twalk.epeak_vec(i+1) & glbl_cut;
    twalk.mean_sat(i) = mean(time_diff_sigmoid(temp_cut));
    twalk.e_peak(i) =  mean(e_peak_MM(temp_cut));
    twalk.npts(i) = sum(temp_cut);
    twalk.rms(i) = std(time_diff_sigmoid(temp_cut));
    twalk.err(i) = std(time_diff_sigmoid(temp_cut))./sqrt(twalk.npts(i)); % mean error
    twalk.e_peak_err_p(i) = twalk.e_peak(i)-(twalk.epeak_vec(i+1)); % e charge limits
    twalk.e_peak_err_n(i) = -twalk.e_peak(i)+(twalk.epeak_vec(i));
end

% fit correction function using minuit
fit_data = [];
fit_data(1,:) = twalk.e_peak;
fit_data(2,:) = twalk.mean_sat;
fit_data(3,:) = twalk.err;
p0=[];
p0(1) = min(twalk.mean_sat);
p0(2) = 1;
p0(3) = 0.5;
cmd='min; ret';
[p, err, chi] = fminuit('twalk_fn_minuit',p0,fit_data,'-b','-c',cmd);
twalk.p = p;
twalk.chi = chi;

% plot time walk correction and fit
figure
errorbar(twalk.e_peak,twalk.mean_sat,twalk.err,twalk.err,twalk.e_peak_err_n,twalk.e_peak_err_p,'o');
hold on
plot(twalk.e_peak,twalk_fn_minuit(p,twalk.e_peak),'LineWidth',1.5);
xlabel('Electron peak charge, pC')
ylabel('SAT, ns')
title_str = sprintf('%s \n SAT vs. Electron peak charge',runTitleString);
title(title_str)

grid
movegui(gcf,'northwest');
saveas(gcf,[store_folder '\Run' run.id '_timewalk.png'])

%plot resolution vs. e-charge
figure
errorbar(twalk.e_peak,twalk.rms*1000,[],[],twalk.e_peak_err_n,twalk.e_peak_err_p,'o');
hold on
xlabel('Electron peak charge, pC')
ylabel('Resolution, ps')
title_str = sprintf('%s \n Resolution vs. Electron peak charge',runTitleString);
title(title_str)

grid
saveas(gcf,[store_folder '\Run' run.id '_res_vs_charge.png'])

% %to avoid negative e-peak charge problem
% negative_charge_exclusion_cut = e_peak_MM < 0;
% e_peak_MM(negative_charge_exclusion_cut) = 0.0000001;

% make time walk correction
if(twalk.en == 1)
    time_diff_sigmoidTWCorr = time_diff_sigmoid - twalk_fn_minuit(twalk.p, e_peak_MM);
else
    time_diff_sigmoidTWCorr = time_diff_sigmoid;
end

%% calculate DUT  centre if alignToDUT == false

if alignToDUT == true

    % calculate pad center
    pad.xc_med = median(trackerX(glbl_cut));
    pad.yc_med = median(trackerY(glbl_cut));

    %alternative naive mean with cut on MM pad for scanning runs
    pad.x_idx_cut = trackerX > (pad.xc_med-pad.cut_len) & trackerX < (pad.xc_med+pad.cut_len) & glbl_cut & MM_maxy>0.5*max(MM_maxy);
    pad.xc_n = mean(trackerX(pad.x_idx_cut)); % naive mean for x
    pad.y_idx_cut = trackerY > (pad.yc_med-pad.cut_len) & trackerY < (pad.yc_med+pad.cut_len) & glbl_cut & MM_maxy>0.5*max(MM_maxy);
    pad.yc_n = mean(trackerY(pad.y_idx_cut)); % naive mean for y

    %std naive mean
    %     pad.x_idx_cut = trackerX > (pad.xc_med-pad.cut_len) & trackerX < (pad.xc_med+pad.cut_len) & glbl_cut;
    %     pad.xc_n = mean(trackerX(pad.x_idx_cut)); % naive mean for x
    %     pad.y_idx_cut = trackerY > (pad.yc_med-pad.cut_len) & trackerY < (pad.yc_med+pad.cut_len) & glbl_cut;
    %     pad.yc_n = mean(trackerY(pad.y_idx_cut)); % naive mean for y




    % Plot x and y hit histograms
    figure;
    title(['%s - Hits on DUT detector (X, Y projections) ' runTitleString]);
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
    
    saveas(gcf,[store_folder '\Run' run.id '_hits_DUT.png'])

    % find x center from electron peak mean charge
    for i = 1:length(pad.h_x.Values)
        tmp_cut = trackerX>pad.h_x.BinEdges(i) & trackerX<pad.h_x.BinEdges(i)+ pad.h_x.BinWidth;
        %tmp_cut = trackerX>pad.h_x.BinEdges(i) & trackerX<pad.h_x.BinEdges(i)+ pad.h_x.BinWidth & MM_maxy>0.3*max(MM_maxy);
        pad.epeak_x(i) = mean(e_peak_MM(tmp_cut));
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
    figure
    bar(fit_data(1,:),fit_data(2,:));
    hold on
    plot(fit_data(1,:),parabola4_minuit(p, fit_data(1,:)),'LineWidth',2);
    xlabel('x-axis, mm');
    ylabel('Charge, pC');
    legend('RAW', 'Fit')
    title_str = sprintf('%s \n E-peak mean over x-axis',runTitleString);
    title(title_str)

    grid on
    saveas(gcf,[store_folder '/Run' run.id '_DUT_epeak_X_proj.png'])
    p
    % pause(10);

    % find y center from electron peak mean
    for i = 1:length(pad.h_y.Values)
        tmp_cut = trackerY>pad.h_y.BinEdges(i) & trackerY<pad.h_y.BinEdges(i)+pad.h_y.BinWidth;
        pad.epeak_y(i) = mean(e_peak_MM(tmp_cut));
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
    figure
    bar(fit_data(1,:),fit_data(2,:));
    hold on
    plot(fit_data(1,:),parabola4_minuit(p, fit_data(1,:)),'LineWidth',2);
    xlabel('y-axis, mm');
    ylabel('Charge, pC');
    legend('RAW', 'Fit')
    title_str = sprintf('%s \n E-peak mean over y-axis',runTitleString);
    title(title_str)

    grid on
    saveas(gcf,[store_folder '\Run' run.id '_DUT_epeak_Y_proj.png'])


    if shouldCenterPadWithFit==false
        pad.xc = pad.xc_n;
        pad.yc = pad.yc_n;
    end
end

%% calculate REF MCP centre

if alignToDUT == false

    % calculate pad center
    pad.xc_med = median(trackerX(glbl_cut));
    pad.yc_med = median(trackerY(glbl_cut));
    pad.x_idx_cut = trackerX > (pad.xc_med-pad.cut_len) & trackerX < (pad.xc_med+pad.cut_len) & glbl_cut;
    pad.xc_n = mean(trackerX(pad.x_idx_cut)); % naive mean for x
    pad.y_idx_cut = trackerY > (pad.yc_med-pad.cut_len) & trackerY < (pad.yc_med+pad.cut_len) & glbl_cut;
    pad.yc_n = mean(trackerY(pad.y_idx_cut)); % naive mean for y

    % Plot x and y hit histograms
    figure;
    title(['Hits on REF detector (X, Y projections)' runTitleString]);
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

    saveas(gcf,[store_folder '\Run' run.id '_hits_REF.png'])

    % find x center from electron peak mean charge
    for i = 1:length(pad.h_x.Values)
        tmp_cut = trackerX>pad.h_x.BinEdges(i) & trackerX<pad.h_x.BinEdges(i)+ pad.h_x.BinWidth;
        pad.epeak_x_mcp(i) = mean(e_peak_MCP(tmp_cut));
    end

    fit_data = [];
    fit_data(1,:) = pad.h_x.BinEdges(1:end-1)+pad.h_x.BinWidth/2;
    fit_data(2,:) = pad.epeak_x_mcp;
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
    fitResultX = p(1)
    p
    pad.xc_err = err(1);
    % plot to see how parabolic fit looks like
    figure
    bar(fit_data(1,:),fit_data(2,:));
    hold on
    plot(fit_data(1,:),parabola4_minuit(p, fit_data(1,:)),'LineWidth',2);
    xlabel('x-axis, mm');
    ylabel('Charge, pC');
    legend('RAW', 'Fit')
    title_str = sprintf('%s \n E-peak REF MCP mean over x-axis',runTitleString);
    title(title_str)
    grid on
    saveas(gcf,[store_folder '/Run' run.id '_REF_epeak_X_proj.png'])


    % find y center from electron peak mean
    for i = 1:length(pad.h_y.Values)
        tmp_cut = trackerY>pad.h_y.BinEdges(i) & trackerY<pad.h_y.BinEdges(i)+pad.h_y.BinWidth;
        pad.epeak_y_mcp(i) = mean(e_peak_MCP(tmp_cut));
    end

    fit_data = [];
    fit_data(1,:) = pad.h_y.BinEdges(1:end-1)+pad.h_y.BinWidth/2;
    fit_data(2,:) = pad.epeak_y_mcp;
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

    fitResultX = p(1);
    p;
    pad.yc_err = err(1);
    % plot to see how parabolic fit looks like
    figure
    bar(fit_data(1,:),fit_data(2,:));
    hold on
    plot(fit_data(1,:),parabola4_minuit(p, fit_data(1,:)),'LineWidth',2);
    xlabel('y-axis, mm');
    ylabel('Charge, pC');
    legend('RAW', 'Fit')
    title_str = sprintf('%s \n E-peak REF MCP mean over y-axis',runTitleString);
    title(title_str)

    grid on
    saveas(gcf,[store_folder '\Run' run.id '_REF_epeak_Y_proj.png'])

    %pad.xc = pad.xc_n;
    %pad.yc = pad.yc_n;

end

if isnan(pad.xc)
    pad.xc=25;
end

if isnan(pad.yc)
    pad.yc=25;
end

% %override pad centering
%   pad.xc=76;
%   pad.yc=78.8;


lineAnglesArray = [0,45,90];

%% 2D PLOTS FOR SAT VS E_PEAK
figure
histogram2(e_peak_MM(glbl_cut),time_diff_sigmoid(glbl_cut), 'NumBins', [100 100], 'DisplayStyle', 'tile', 'ShowEmptyBins', 'on');
xlabel('Electron peak charge, pC')
ylabel('SAT, ns')
title_str = sprintf('%s \n SAT vs. Electron peak charge (glbl cut - 2D plot)',runTitleString);
title(title_str)

grid
movegui(gcf,'northwest');
saveas(gcf,[store_folder '\Run' run.id '_2Dplot_e_peak_vs_SAT_glbl_cut.png']);

figure
low_e_peak_cut = e_peak_MM < 4 & glbl_cut;
medium_e_peak_cut = e_peak_MM > 4 & e_peak_MM < 10 & glbl_cut;
high_e_peak_cut = e_peak_MM > 10 & glbl_cut;
histogram(time_diff_sigmoid(low_e_peak_cut), 50);
hold on;
histogram(time_diff_sigmoid(medium_e_peak_cut), 50);
histogram(time_diff_sigmoid(high_e_peak_cut), 50);
xlabel('SAT, ns')
ylabel('Events')
title_str = sprintf('%s \n Projections of SAT vs. Electron peak charge (glbl cut - 2D plot)',runTitleString);
title(title_str)
legend('below 4 pC','between 4 and 10 pC','over 10 pC');

grid
movegui(gcf,'northwest');
saveas(gcf,[store_folder '\Run' run.id '_2Dplot_e_peak_vs_SAT_glbl_cut_projections_below4_middle_over10pC.png']);


if enable_cut_twalk
    %circular cut_small - sampling in center
    cut_sampling = (((trackerX - pad.xc-shiftX).^2 + (trackerY - pad.yc-shiftY).^2) < small.r^2) & glbl_cut;
    
    % Time walk analysis (this needs to be improoved)
    twalk.n_epk = 50;                               % number of e-peak bins
    e_peak_srt = sort(e_peak_MM(cut_sampling));         % make temporary sort of e-peaks
    % try to distrubute e_peak vector evenly over the amplitude range
    twalk.epeak_vec = e_peak_srt(1:round(length(e_peak_srt)/twalk.n_epk):end); % 50 bins
    for i=1:length(twalk.epeak_vec)-1
        temp_cut = e_peak_MM > twalk.epeak_vec(i) & e_peak_MM < twalk.epeak_vec(i+1) & cut_sampling;
        twalk.mean_sat(i) = mean(time_diff_sigmoid(temp_cut));
        twalk.e_peak(i) =  mean(e_peak_MM(temp_cut));
        twalk.npts(i) = sum(temp_cut);
        twalk.rms(i) = std(time_diff_sigmoid(temp_cut));
        twalk.err(i) = std(time_diff_sigmoid(temp_cut))./sqrt(twalk.npts(i)); % mean error
        twalk.e_peak_err_p(i) = twalk.e_peak(i)-(twalk.epeak_vec(i+1)); % e charge limits
        twalk.e_peak_err_n(i) = -twalk.e_peak(i)+(twalk.epeak_vec(i));
    end
    
    % fit correction function using minuit
    fit_data = [];
    fit_data(1,:) = twalk.e_peak;
    fit_data(2,:) = twalk.mean_sat;
    fit_data(3,:) = twalk.err;
    p0=[];
    p0(1) = min(twalk.mean_sat);
    p0(2) = 1;
    p0(3) = 0.5;
    cmd='min; ret';
    [p, err, chi] = fminuit('twalk_fn_minuit',p0,fit_data,'-b','-c',cmd);
    twalk.p = p;
    twalk.chi = chi;
    
    % plot time walk correction and fit
    figure
    errorbar(twalk.e_peak,twalk.mean_sat,twalk.err,twalk.err,twalk.e_peak_err_n,twalk.e_peak_err_p,'o');
    hold on
    plot(twalk.e_peak,twalk_fn_minuit(p,twalk.e_peak),'LineWidth',1.5);
    xlabel('Electron peak charge, pC')
    ylabel('SAT, ns')
    title_str = sprintf('%s \n SAT vs. Electron peak charge',runTitleString);
    title(title_str)
    
    grid
    movegui(gcf,'northwest');
    saveas(gcf,[store_folder '\Run' run.id '_timewalk_cut_sampling.png'])


        % make time walk correction
    if(twalk.en == 1)
        time_diff_sigmoidTWCorr = time_diff_sigmoid - twalk_fn_minuit(twalk.p, e_peak_MM);
    else
        time_diff_sigmoidTWCorr = time_diff_sigmoid;
    end

    figure
    histogram2(e_peak_MM(cut_sampling),time_diff_sigmoid(cut_sampling), 'NumBins', [100 100], 'DisplayStyle', 'tile', 'ShowEmptyBins', 'on');
    xlabel('Electron peak charge, pC')
    ylabel('SAT, ns')
    title_str = sprintf('%s \n SAT vs. Electron peak charge (cut sampling - 2D plot)',runTitleString);
    title(title_str)
    
    grid
    movegui(gcf,'northwest');
    saveas(gcf,[store_folder '\Run' run.id '_2Dplot_e_peak_vs_SAT_cut_sampling.png']);
    
    figure
    low_e_peak_cut = e_peak_MM < 4 & cut_sampling;
    medium_e_peak_cut = e_peak_MM > 4 & e_peak_MM < 10 & cut_sampling;
    high_e_peak_cut = e_peak_MM > 10 & cut_sampling;
    histogram(time_diff_sigmoid(low_e_peak_cut), 50);
    hold on;
    histogram(time_diff_sigmoid(medium_e_peak_cut), 50);
    histogram(time_diff_sigmoid(high_e_peak_cut), 50);
    xlabel('SAT, ns')
    ylabel('Events')
    title_str = sprintf('%s \n Projections of SAT vs. Electron peak charge (cut sampling - 2D plot)',runTitleString);
    title(title_str)
    legend('below 4 pC','between 4 and 10 pC','over 10 pC');
    
    grid
    movegui(gcf,'northwest');
    saveas(gcf,[store_folder '\Run' run.id '_2Dplot_e_peak_vs_SAT_cut_sampling_projections_below4_middle_over10pC.png']);
end


%% efficiency measurement
%data for efficiency is stored in arrays containing all events (all that
%were saved by oscilloscope)
%eventsValidREF, eventsValidDUT, eventsValidTracker, eventsTrackerX, eventsTrackerY

%select events that have valid tracker data
geoCutActiveArea = (((eventsTrackerX - pad.xc-shiftX).^2 + (eventsTrackerY - pad.yc-shiftY).^2) < (pad.size/2)^2);
geoCutSamplingArea = (((eventsTrackerX - pad.xc-shiftX).^2 + (eventsTrackerY - pad.yc-shiftY).^2) < small.r^2);

efficiencyTrackerCut = eventsValidTracker & eventsTrackerX~=0 & eventsTrackerY~=0;

numberAllEvents = length(eventsValidTracker)
numberEventsWithTrackerActiveArea = length(eventsValidTracker(efficiencyTrackerCut & geoCutActiveArea))
numberEventsWithTrackerSamplingArea = length(eventsValidTracker(efficiencyTrackerCut & geoCutSamplingArea))

refEventsTrackerValidCutActiveArea = efficiencyTrackerCut & eventsValidREF & geoCutActiveArea;
refEventsTrackerValidActiveArea = eventsValidREF(refEventsTrackerValidCutActiveArea);
numberAcceptedRefEventsActiveArea = length(refEventsTrackerValidActiveArea);

dutEventsTrackerValidCutActiveArea = efficiencyTrackerCut & eventsValidDUT & geoCutActiveArea;
dutEventsTrackerValidActiveArea = eventsValidDUT(dutEventsTrackerValidCutActiveArea);
numberAcceptedDUTEventsActiveArea = length(dutEventsTrackerValidActiveArea);
dutEventsTrackerValidXActiveArea = eventsTrackerX(dutEventsTrackerValidCutActiveArea);
dutEventsTrackerValidYActiveArea = eventsTrackerY(dutEventsTrackerValidCutActiveArea);
;
efficiencyREFActiveArea = numberAcceptedRefEventsActiveArea/numberEventsWithTrackerActiveArea
efficiencyDUTActiveArea = numberAcceptedDUTEventsActiveArea/numberEventsWithTrackerActiveArea

refEventsTrackerValidCutSamplingArea = efficiencyTrackerCut & eventsValidREF & geoCutSamplingArea;
refEventsTrackerValidSamplingArea = eventsValidREF(refEventsTrackerValidCutSamplingArea);
numberAcceptedRefEventsSamplingArea = length(refEventsTrackerValidSamplingArea);

dutEventsTrackerValidCutSamplingArea = efficiencyTrackerCut & eventsValidDUT & geoCutSamplingArea;
dutEventsTrackerValidSamplingArea = eventsValidDUT(dutEventsTrackerValidCutSamplingArea);
numberAcceptedDUTEventsSamplingArea = length(dutEventsTrackerValidSamplingArea);
dutEventsTrackerValidXSamplingArea = eventsTrackerX(dutEventsTrackerValidCutSamplingArea);
dutEventsTrackerValidYSamplingArea = eventsTrackerY(dutEventsTrackerValidCutSamplingArea);

efficiencyREFSamplingArea = numberAcceptedRefEventsSamplingArea/numberEventsWithTrackerSamplingArea
efficiencyDUTSamplingArea = numberAcceptedDUTEventsSamplingArea/numberEventsWithTrackerSamplingArea

%% repeat efficiencyCalculation with cuts on timing and amp as glb cut

%acceptance cut incl. cut on REF - used previously
% acceptanceCut = eventsTimeDiffSigmoid>time_min & eventsTimeDiffSigmoid<time_max & ...
%     eventsREFmaxY>0.01*max(eventsREFmaxY) & eventsDUTmaxY>0.01*max(eventsDUTmaxY)&...
%     eventsDUTmaxY<0.95*max(eventsDUTmaxY)& eventsREFmaxY<0.95*max(eventsREFmaxY);

%acceptance cut without cut on REF
acceptanceCut = eventsTimeDiffSigmoid>time_min & eventsTimeDiffSigmoid<time_max & ...
    eventsDUTmaxY>0.01*max(eventsDUTmaxY);

%acceptance cut - only timing cut
acceptanceCut = eventsTimeDiffSigmoid>time_min & eventsTimeDiffSigmoid<time_max;


geoCutActiveArea = (((eventsTrackerX - pad.xc-shiftX).^2 + (eventsTrackerY - pad.yc-shiftY).^2) < (pad.size/2)^2);
geoCutSamplingArea = (((eventsTrackerX - pad.xc-shiftX).^2 + (eventsTrackerY - pad.yc-shiftY).^2) < small.r^2);
efficiencyTrackerCut = eventsValidTracker & eventsTrackerX~=0 & eventsTrackerY~=0 ;


numberAllEvents = length(eventsValidTracker)
numberEventsWithTrackerActiveArea = length(eventsValidTracker(efficiencyTrackerCut & geoCutActiveArea))
numberEventsWithTrackerSamplingArea = length(eventsValidTracker(efficiencyTrackerCut & geoCutSamplingArea))

refEventsTrackerValidCutActiveArea = efficiencyTrackerCut & eventsValidREF & geoCutActiveArea & acceptanceCut;
refEventsTrackerValidActiveArea = eventsValidREF(refEventsTrackerValidCutActiveArea);
numberAcceptedRefEventsActiveArea = length(refEventsTrackerValidActiveArea);

dutEventsTrackerValidCutActiveArea = efficiencyTrackerCut & eventsValidDUT & geoCutActiveArea & acceptanceCut;
dutEventsTrackerValidActiveArea = eventsValidDUT(dutEventsTrackerValidCutActiveArea);
numberAcceptedDUTEventsActiveArea = length(dutEventsTrackerValidActiveArea);
dutEventsTrackerValidXActiveArea = eventsTrackerX(dutEventsTrackerValidCutActiveArea);
dutEventsTrackerValidYActiveArea = eventsTrackerY(dutEventsTrackerValidCutActiveArea);
;
efficiencyREFActiveAreaAccepted = numberAcceptedRefEventsActiveArea/numberEventsWithTrackerActiveArea
efficiencyDUTActiveAreaAccepted = numberAcceptedDUTEventsActiveArea/numberEventsWithTrackerActiveArea

refEventsTrackerValidCutSamplingArea = efficiencyTrackerCut & eventsValidREF & geoCutSamplingArea & acceptanceCut;
refEventsTrackerValidSamplingArea = eventsValidREF(refEventsTrackerValidCutSamplingArea);
numberAcceptedRefEventsSamplingArea = length(refEventsTrackerValidSamplingArea);

dutEventsTrackerValidCutSamplingArea = efficiencyTrackerCut & eventsValidDUT & geoCutSamplingArea & acceptanceCut;
dutEventsTrackerValidSamplingArea = eventsValidDUT(dutEventsTrackerValidCutSamplingArea);
numberAcceptedDUTEventsSamplingArea = length(dutEventsTrackerValidSamplingArea);
dutEventsTrackerValidXSamplingArea = eventsTrackerX(dutEventsTrackerValidCutSamplingArea);
dutEventsTrackerValidYSamplingArea = eventsTrackerY(dutEventsTrackerValidCutSamplingArea);

efficiencyREFSamplingAreaAccepted = numberAcceptedRefEventsSamplingArea/numberEventsWithTrackerSamplingArea
efficiencyDUTSamplingAreaAccepted = numberAcceptedDUTEventsSamplingArea/numberEventsWithTrackerSamplingArea


%% plot events which were excluded
exludedEventsCut = not(dutEventsTrackerValidCutActiveArea) & eventsValidTracker;



% plot hits
figure;
scatter(eventsTrackerX(exludedEventsCut), eventsTrackerY(exludedEventsCut),'.')   % plot all hits
hold on
rectangle('Position',[pad.xc-pad.size/2 pad.yc-pad.size/2 pad.size pad.size], 'Curvature',pad.curvature, 'LineWidth',1);
rectangle('Position',[pad.xc-shiftX-(2*small.r)/2 pad.yc-shiftY-(2*small.r)/2 (2*small.r) (2*small.r)], 'Curvature',[1,1], 'LineWidth',2,'EdgeColor','red');    %plot selected area
axis equal
xlabel('x-axis, mm');
ylabel('y-axis, mm');
grid on
xlim([pad.xc-pad.size pad.xc+pad.size]);
ylim([pad.yc-pad.size pad.yc+pad.size]);
str_title = sprintf('%s \n inefficient hits', runTitleString);
title(str_title);
legend('Trigger','Hits','Observation path','DUT center','Circle for resolution measurement');
movegui(gcf,'north');
saveas(gcf,[store_folder '\Run' run.id '_inefficientHits.png'])




pause(1);








% plot hits
figure;
scatter(trackerX, trackerY,'.')   % plot all hits
hold on;
scatter(trackerX(glbl_cut), trackerY(glbl_cut),'.'); % plot hits with global cut
scatter(pad.xc,pad.yc);                            % plot center point
rectangle('Position',[pad.xc-pad.size/2 pad.yc-pad.size/2 pad.size pad.size], 'Curvature',pad.curvature, 'LineWidth',1);    %plot active area
rectangle('Position',[pad.xc-shiftX-(2*small.r)/2 pad.yc-shiftY-(2*small.r)/2 (2*small.r) (2*small.r)], 'Curvature',[1,1], 'LineWidth',2,'EdgeColor','red');    %plot selected area
axis equal
xlabel('x-axis, mm');
ylabel('y-axis, mm');
grid on
str_title = sprintf('%s \n tracker data', runTitleString);
title(str_title);
legend('Trigger','Hits','Observation path','DUT center','Circle for resolution measurement');
if length(trackerX(glbl_cut))>0
    xlim([min(trackerX(glbl_cut)) max(trackerX(glbl_cut))]);
    ylim([min(trackerY(glbl_cut)) max(trackerY(glbl_cut))]);
end
movegui(gcf,'north');
saveas(gcf,[store_folder '\Run' run.id '_hits_overlay.png'])



%%loop over different angles and produce plots at 0/45/90deg

lineAnglesArray = [0,45,90,135,180,225,270,315];
satLinesX = [];
satLinesY = [];
satLinesErr = [];

for linePos = 1:length(lineAnglesArray)

    %% generate observation line and shift to MCP center position
    line.angle = lineAnglesArray(linePos);     % line angle
    line.step = 0.5;     % spacing in mm
    line.half_len = 5.5;   % half length of the line
    line.pts = -line.half_len:line.step:line.half_len;
    line.x = cos(line.angle*pi/180)*line.pts + pad.xc;
    line.y = sin(line.angle*pi/180)*line.pts + pad.yc;
    line.n = length(line.pts);
    plot(line.x, line.y, 'LineWidth', 2);              % plot observation line

end



for linePos = 1:length(lineAnglesArray)

    %% generate observation line and shift to MCP center position
    line.angle = lineAnglesArray(linePos);     % line angle
    line.step = 0.5;     % spacing in mm
    line.half_len = 5.5;   % half length of the line
    line.pts = -line.half_len:line.step:line.half_len;
    line.x = cos(line.angle*pi/180)*line.pts + pad.xc;
    line.y = sin(line.angle*pi/180)*line.pts + pad.yc;
    line.n = length(line.pts);

    % calculate resolution over the line
    line.radius = 1;
    for i = 1:line.n
        cut_circ = ((trackerX - line.x(i)).^2 + (trackerY - line.y(i)).^2) < line.radius^2 & glbl_cut;
        line.sat(i) = mean(time_diff_sigmoidTWCorr(cut_circ));
        line.sat_pts(i) = sum(cut_circ);
        line.sat_err(i) = std(time_diff_sigmoidTWCorr(cut_circ))./sqrt(line.sat_pts(i));
    end

    % plot SAT over observation line
    figure
    errorbar(line.pts,1000*(line.sat),line.sat_err*1000,'o', 'LineWidth',1);
    str_title = sprintf('%s:\n SAT line at %2.1f°',runTitleString, line.angle);
    title(str_title);
    xlabel('Distance from pad centre, mm');
    ylabel('SAT, ps');
    satLines.legend{linePos} = ['Cutting angle ' num2str(lineAnglesArray(linePos)) '°'];

    grid on;
    saveas(gcf,[store_folder '\Run' run.id '_SAT_over_' num2str(line.angle) 'deg_line.png'])

    satLinesX = [satLinesX;line.pts];
    satLinesY = [satLinesY;1000*(line.sat)];
    satLinesErr = [satLinesErr;1000*(line.sat_err)];
end


%combined plot
for linePos = 1:length(lineAnglesArray)
    % plot SAT over observation line
    if linePos==1

        figure
    end
    errorbar(satLinesX(linePos,:),satLinesY(linePos,:),satLinesErr(linePos,:),'o', 'LineWidth',1);
    str_title = sprintf('%s:\n SAT Comparison',runTitleString);
    title(str_title);
    xlabel('Distance from pad centre, mm');
    ylabel('SAT, ps');
    grid on;

    if linePos==1
        hold on
    end

end
legend(satLines.legend);
saveas(gcf,[store_folder '\Run' run.id '_SAT_linesComparison.png'])



area.x_vec =-area.size:area.step:area.size+area.step;    % define area of interest for 2D plots in X axis
area.y_vec =-area.size:area.step:area.size+area.step;    % define area of interest for 2D plots in X axis

% move rectangle to pad center
area.x_vec = area.x_vec + pad.xc - area.step/2;
area.y_vec = area.y_vec + pad.yc - area.step/2;

% ndgrid plot
[area.xx, area.yy] = ndgrid(area.x_vec, area.y_vec);
for i=1:length(area.x_vec)
    for j=1:length(area.y_vec)
        % make moving circular cut that respects the global cut (idx_cut)
        cut_circ = ((trackerX - area.x_vec(i)).^2 + (trackerY - area.y_vec(j)).^2) < area.radius^2 & glbl_cut;
        area.sat(i,j) = mean(time_diff_sigmoidTWCorr(cut_circ));
        area.rms(i,j) = std(time_diff_sigmoidTWCorr(cut_circ));
        area.amp(i,j) = mean(MM_maxy(cut_circ));
        area.amp_MCP(i,j) = mean(MCP_maxy(cut_circ));
        area.e_peak(i,j) = mean(e_peak_MM(cut_circ));
        area.sat_pts(i,j) = sum(cut_circ);
        area.sat_err(i,j) = std(time_diff_sigmoidTWCorr(cut_circ))./sqrt(area.sat_pts(i,j));

        %plot for efficiencyMap
        cut_circEfficiency = ((eventsTrackerX - area.x_vec(i)).^2 + (eventsTrackerY - area.y_vec(j)).^2) < area.radius^2 & efficiencyTrackerCut;
        cut_circEfficiencyAcceptance = ((eventsTrackerX - area.x_vec(i)).^2 + (eventsTrackerY - area.y_vec(j)).^2) < area.radius^2 & efficiencyTrackerCut & acceptanceCut & eventsValidDUT==1;
       % area.dutEfficiency(i,j) = mean(eventsValidDUT(cut_circEfficiency));

        area.dutEfficiency(i,j) = length(eventsValidTracker(cut_circEfficiencyAcceptance)) / length(eventsValidTracker(cut_circEfficiency));

        
        area.refEfficiency(i,j) = mean(eventsValidREF(cut_circEfficiency));
    end
end


%% plot efficiency
figure
h=pcolor(area.xx,area.yy,area.dutEfficiency);
hold on
rectangle('Position',[pad.xc-pad.size/2 pad.yc-pad.size/2 pad.size pad.size], 'Curvature',pad.curvature, 'LineWidth',1);
rectangle('Position',[pad.xc-shiftX-(2*small.r)/2 pad.yc-shiftY-(2*small.r)/2 (2*small.r) (2*small.r)], 'Curvature',[1,1], 'LineWidth',2,'EdgeColor','red');    %plot selected area
set(h, 'EdgeColor', 'none');
%set(gca, 'XDir','reverse')
axis equal
colorbar
xlabel('x-axis, mm');
ylabel('y-axis, mm');
h = colorbar;
h.Label.String = 'Efficiency DUT';
caxis([0 1]);
h.Label.Position(1) = 3;
str_title = sprintf('%s: DUT Efficiency \n All Events: ActiveArea: %.1f%% | CenterSampling: %.1f%% \n Included TimeRes: ActiveArea: %.1f%% | CenterSampling: %.1f%%', runTitleString, efficiencyDUTActiveArea*100, efficiencyDUTSamplingArea*100, efficiencyDUTActiveAreaAccepted*100, efficiencyDUTSamplingAreaAccepted*100);
title(str_title);
saveas(gcf,[store_folder '\Run' run.id '_EfficiencyDUT.png'])
movegui(gcf,'northeast');
pause(1);

%%plot efficiency ref
figure
h=pcolor(area.xx,area.yy,area.refEfficiency);
hold on
rectangle('Position',[pad.xc-pad.size/2 pad.yc-pad.size/2 pad.size pad.size], 'Curvature',pad.curvature, 'LineWidth',1);
rectangle('Position',[pad.xc-shiftX-(2*small.r)/2 pad.yc-shiftY-(2*small.r)/2 (2*small.r) (2*small.r)], 'Curvature',[1,1], 'LineWidth',2,'EdgeColor','red');    %plot selected area
set(h, 'EdgeColor', 'none');
%set(gca, 'XDir','reverse')
axis equal
colorbar
xlabel('x-axis, mm');
ylabel('y-axis, mm');
h = colorbar;
h.Label.String = 'Efficiency REF';
caxis([0 1]);
h.Label.Position(1) = 3;
str_title = sprintf('%s: REF Efficiency \n All Events: ActiveArea: %.1f%% | CenterSampling: %.1f%% \n Included TimeRes: ActiveArea: %.1f%% | CenterSampling: %.1f%%', runTitleString, efficiencyREFActiveArea*100, efficiencyREFSamplingArea*100, efficiencyREFActiveAreaAccepted*100, efficiencyREFSamplingAreaAccepted*100);
title(str_title);
saveas(gcf,[store_folder '\Run' run.id '_EfficiencyREF.png'])
pause(1);

%% plot SAT over area
figure
h=pcolor(area.xx,area.yy,area.sat);
hold on
rectangle('Position',[pad.xc-pad.size/2 pad.yc-pad.size/2 pad.size pad.size], 'Curvature',pad.curvature, 'LineWidth',1);
rectangle('Position',[pad.xc-shiftX-(2*small.r)/2 pad.yc-shiftY-(2*small.r)/2 (2*small.r) (2*small.r)], 'Curvature',[1,1], 'LineWidth',2,'EdgeColor','red');    %plot selected area
set(h, 'EdgeColor', 'none');
%set(gca, 'XDir','reverse')
axis equal
colorbar
xlabel('x-axis, mm');
ylabel('y-axis, mm');
h = colorbar;
h.Label.String = 'SAT, ns';
h.Label.Position(1) = 3;
str_title = sprintf('%s:\n SAT map (\\phi_{avg} = %2.1f mm)', runTitleString, 2*area.radius);
title(str_title);
saveas(gcf,[store_folder '\Run' run.id '_SAT_map.png'])

%% plot time resolution over area
figure
h=pcolor(area.xx,area.yy,area.rms*1000);
hold on
rectangle('Position',[pad.xc-pad.size/2 pad.yc-pad.size/2 pad.size pad.size], 'Curvature',pad.curvature, 'LineWidth',1);
rectangle('Position',[pad.xc-shiftX-(2*small.r)/2 pad.yc-shiftY-(2*small.r)/2 (2*small.r) (2*small.r)], 'Curvature',[1,1], 'LineWidth',2,'EdgeColor','red');    %plot selected area
set(h, 'EdgeColor', 'none');
axis equal
colorbar
xlabel('x-axis, mm');
ylabel('y-axis, mm');
h = colorbar;
caxis([20 60]);
h.Label.String = 'Time resolution, ps';
h.Label.Position(1) = 3;
%str_title = sprintf('%s:\n Time resolution map (\\phi_{avg} = %2.1f mm)', runTitleString, 2*area.radius);
str_title = sprintf('Time resolution map') % (\\phi_{avg} = %2.1f mm)', runTitleString, 2*area.radius);
title(str_title);

movegui(gcf,'south');
saveas(gcf,[store_folder '\Run' run.id '_RES_map.png'])

%% plot SAT mean errors
figure
h=pcolor(area.xx,area.yy,area.sat_err*1000);
hold on
rectangle('Position',[pad.xc-pad.size/2 pad.yc-pad.size/2 pad.size pad.size],  'Curvature',pad.curvature,'LineWidth',1);
rectangle('Position',[pad.xc-shiftX-(2*small.r)/2 pad.yc-shiftY-(2*small.r)/2 (2*small.r) (2*small.r)], 'Curvature',[1,1], 'LineWidth',2,'EdgeColor','red');    %plot selected area
set(h, 'EdgeColor', 'none');
axis equal
colorbar
xlabel('x-axis, mm');
ylabel('y-axis, mm');
h = colorbar;
h.Label.String = 'Error, ps';
h.Label.Position(1) = 3;
str_title = sprintf('%s: \n Mean error map (\\phi_{avg} = %2.1f mm)', runTitleString, 2*area.radius);
title(str_title);
saveas(gcf,[store_folder '\Run' run.id '_SAT_err.png'])


%% plot DUT amplitude over area
figure
h=pcolor(area.xx,area.yy,area.amp);
hold on
rectangle('Position',[pad.xc-pad.size/2 pad.yc-pad.size/2 pad.size pad.size],  'Curvature',pad.curvature,'LineWidth',2);
rectangle('Position',[pad.xc-shiftX-(2*small.r)/2 pad.yc-shiftY-(2*small.r)/2 (2*small.r) (2*small.r)], 'Curvature',[1,1], 'LineWidth',2,'EdgeColor','red');    %plot selected area
set(h, 'EdgeColor', 'none');
axis equal
colorbar
xlabel('x-axis, mm');
ylabel('y-axis, mm');
h = colorbar;
h.Label.String = 'Amplitude, V';
h.Label.Position(1) = 3;
str_title = sprintf('%s:\n DUT MM map (\\phi_{avg} = %2.1f mm)', runTitleString, 2*area.radius);
title(str_title);
movegui(gcf,'west');
saveas(gcf,[store_folder '\Run' run.id '_DUT_amp.png'])

%% plot DUT amplitude over area as the 3D plot
figure
h=surf(area.xx,area.yy,area.amp);
hold on
rectangle('Position',[pad.xc-pad.size/2 pad.yc-pad.size/2 pad.size pad.size],  'Curvature',pad.curvature,'LineWidth',2);
rectangle('Position',[pad.xc-shiftX-(2*small.r)/2 pad.yc-shiftY-(2*small.r)/2 (2*small.r) (2*small.r)], 'Curvature',[1,1], 'LineWidth',2,'EdgeColor','red');    %plot selected area
set(h, 'EdgeColor', 'none');
%axis equal
colorbar
xlabel('x-axis, mm');
ylabel('y-axis, mm');
h = colorbar;
h.Label.String = 'Amplitude, V';
h.Label.Position(1) = 3;
str_title = sprintf('%s:\n DUT MM map (\\phi_{avg} = %2.1f mm)', runTitleString, 2*area.radius);
title(str_title);
saveas(gcf,[store_folder '\Run' run.id '_DUT_amp_3D.png'])

%% plot E-peak charge over area
figure
h=pcolor(area.xx,area.yy,area.e_peak);
hold on
rectangle('Position',[pad.xc-pad.size/2 pad.yc-pad.size/2 pad.size pad.size],  'Curvature',pad.curvature, 'LineWidth',2);
rectangle('Position',[pad.xc-shiftX-(2*small.r)/2 pad.yc-shiftY-(2*small.r)/2 (2*small.r) (2*small.r)], 'Curvature',[1,1], 'LineWidth',2,'EdgeColor','red');    %plot selected area
set(h, 'EdgeColor', 'none');
axis equal
colorbar
xlabel('x-axis, mm');
ylabel('y-axis, mm');
h = colorbar;
h.Label.String = 'Charge, pC';
h.Label.Position(1) = 3;
str_title = sprintf('%s:\n DUT E-peak charge (\\phi_{avg} = %2.1f mm)', runTitleString, 2*area.radius);
title(str_title);
saveas(gcf,[store_folder '\Run' run.id '_charge_map.png'])

%% plot REF MCP Amplitude over area
figure
h = pcolor(area.xx,area.yy,area.amp_MCP);
hold on
rectangle('Position',[pad.xc-pad.size/2 pad.yc-pad.size/2 pad.size pad.size], 'Curvature',pad.curvature, 'LineWidth',2);
rectangle('Position',[pad.xc-shiftX-(2*small.r)/2 pad.yc-shiftY-(2*small.r)/2 (2*small.r) (2*small.r)], 'Curvature',[1,1], 'LineWidth',2,'EdgeColor','red');    %plot selected area
set(h, 'EdgeColor', 'none');
axis equal
colorbar
xlabel('x-axis, mm');
ylabel('y-axis, mm');
h = colorbar;
h.Label.String = 'Amplitude, V';
h.Label.Position(1) = 3;
str_title = sprintf('%s:\n REF MCP amplitude map (\\phi_{avg} = %2.1f mm)', runTitleString, 2*area.radius);
title(str_title);
movegui(gcf,'center');
saveas(gcf,[store_folder '\Run' run.id '_REF_MCP_map.png'])

%% plot REF MCP Amplitude over area 3D
figure
h = surf(area.xx,area.yy,area.amp_MCP);
hold on
rectangle('Position',[pad.xc-pad.size/2 pad.yc-pad.size/2 pad.size pad.size], 'Curvature',pad.curvature,'LineWidth',2);
rectangle('Position',[pad.xc-shiftX-(2*small.r)/2 pad.yc-shiftY-(2*small.r)/2 (2*small.r) (2*small.r)], 'Curvature',[1,1], 'LineWidth',2,'EdgeColor','red');    %plot selected area
set(h, 'EdgeColor', 'none');
%axis equal
colorbar
xlabel('x-axis, mm');
ylabel('y-axis, mm');
h = colorbar;
h.Label.String = 'Amplitude, V';
h.Label.Position(1) = 3;
str_title = sprintf('%s:\n REF MCP amplitude map (\\phi_{avg} = %2.1f mm)', runTitleString, 2*area.radius);
title(str_title);
saveas(gcf,[store_folder '\Run' run.id '_REF_MCP_map_3D.png'])



%% plot peak amplitude histogram cut included
figure
histogram(MCP_maxy(glbl_cut),100)
hold on
h=histogram(MM_maxy(glbl_cut),100);
hold on
xbins = h.BinEdges(1:end-1)+h.BinWidth/2;
fit_data = [];
fit_data(1,:)=1:length(xbins);
fit_data(2,:) = h.Values;
% fit_data(3,:) = yerr;
p0=[];
p0(1) = sum(h.Values)*h.BinWidth;   % normalization factor
p0(2) = 1;
p0(3) = 0.1;
cmd='min; ret';
[p, err, chi] = fminuit('polya_minuit',p0,fit_data(:,1:end),'-b','-c',cmd);
polya_cnt_fit = polya_minuit(p,fit_data(1,:));
plot(xbins,polya_cnt_fit,'Linewidth',2);
e_peak_amp.mean = sum(polya_cnt_fit.*xbins)/sum(polya_cnt_fit);
[dummy, e_peak_amp.max_idx] = max(polya_cnt_fit);

xlabel('Signal amplitude, V')
ylabel('Events');
grid on
legend('MCP','MM','MM fit');
title_str = sprintf('%s \n e-peak amplitude \\mu = %4.4f V U_{max} = %4.4f V',runTitleString, e_peak_amp.mean,xbins(e_peak_amp.max_idx));
title(title_str)
movegui(gcf,'east');

saveas(gcf,[store_folder '\Run' run.id '_signalAmplitude_Hist_comparison.png'])


%% plot peak amplitude histogram cut included - Only MM
figure
h=histogram(MM_maxy(glbl_cut),100);
hold on
xbins = h.BinEdges(1:end-1)+h.BinWidth/2;
fit_data = [];
fit_data(1,:)=1:length(xbins);
fit_data(2,:) = h.Values;
% fit_data(3,:) = yerr;
p0=[];
p0(1) = sum(h.Values)*h.BinWidth;   % normalization factor
p0(2) = 1;
p0(3) = 0.1;
cmd='min; ret';
[p, err, chi] = fminuit('polya_minuit',p0,fit_data(:,1:end),'-b','-c',cmd);
polya_cnt_fit = polya_minuit(p,fit_data(1,:));
plot(xbins,polya_cnt_fit,'Linewidth',2);
e_peak_amp.mean = sum(polya_cnt_fit.*xbins)/sum(polya_cnt_fit);
[dummy, e_peak_amp.max_idx] = max(polya_cnt_fit);

xlabel('Signal amplitude, V')
ylabel('Events');
grid on
legend('MM','MM fit');
title_str = sprintf('%s \n e-peak amplitude \\mu = %4.4f V U_{max} = %4.4f V',runTitleString, e_peak_amp.mean,xbins(e_peak_amp.max_idx));
title(title_str)
saveas(gcf,[store_folder '\Run' run.id '_signalAmplitude_Hist_MM.png'])


riseTimeCut = abs(riseTime)<10;

%% plot rise time distribution
figure
h=histogram(riseTime(glbl_cut&riseTimeCut),100);
hold on
xbins = h.BinEdges(1:end-1)+h.BinWidth/2;

xlabel('Rise time (ns)')
ylabel('Events');
grid on
xlim([0 2]);
%legend('MM','MM fit');
title_str = sprintf('%s \n Signal rise time mean = %4.4f ns',runTitleString, mean(riseTime(riseTimeCut&glbl_cut)));
title(title_str)
saveas(gcf,[store_folder '\Run' run.id '_riseTimeHist.png'])





%% plot electron charge
figure;
plot(MM_maxy(glbl_cut), e_peak_MM(glbl_cut),'.');
grid on
xlabel('Signal amplitude, V')
ylabel('Electron charge, pC');
title_str = sprintf('%s \n Electron charge plot',runTitleString);
title(title_str)
saveas(gcf,[store_folder '\Run' run.id '_echarge_vs_ampl.png'])

%% plot electron lead charge histogram

figure;
h=histogram(e_lead_MM(glbl_cut),100);
hold on
xbins = h.BinEdges(1:end-1)+h.BinWidth/2;
fit_data = [];
fit_data(1,:)=1:length(xbins);
fit_data(2,:) = h.Values;
% fit_data(3,:) = yerr;
p0=[];
p0(1) = sum(h.Values)*h.BinWidth;   % normalization factor
p0(2) = 1;
p0(3) = 0.1;
cmd='min; ret';
[p, err, chi] = fminuit('polya_minuit',p0,fit_data(:,1:end),'-b','-c',cmd);
polya_cnt_fit = polya_minuit(p,fit_data(1,:));
plot(xbins,polya_cnt_fit,'Linewidth',2);
e_peak_amp.mean = sum(polya_cnt_fit.*xbins)/sum(polya_cnt_fit);
[~, e_peak_amp.max_idx] = max(polya_cnt_fit);

xlabel('Electron lead charge, pC');
ylabel('Events');
grid on
legend('Hist','Polya fit');
title_str = sprintf('%s \n e-lead charge \\mu = %4.4f pC Q_{max} = %4.4f pC',runTitleString, e_peak_amp.mean,xbins(e_peak_amp.max_idx));
title(title_str)
saveas(gcf,[store_folder '\Run' run.id '_eleadcharge_hist.png'])


%% plot electron peak charge histogram
figure;
h=histogram(e_peak_MM(glbl_cut),100);
hold on
xbins = h.BinEdges(1:end-1)+h.BinWidth/2;
fit_data = [];
fit_data(1,:)=1:length(xbins);
fit_data(2,:) = h.Values;
% fit_data(3,:) = yerr;
p0=[];
p0(1) = sum(h.Values)*h.BinWidth;   % normalization factor
p0(2) = 1;
p0(3) = 0.1;
cmd='min; ret';
[p, err, chi] = fminuit('polya_minuit',p0,fit_data(:,1:end),'-b','-c',cmd);
polya_cnt_fit = polya_minuit(p,fit_data(1,:));
plot(xbins,polya_cnt_fit,'Linewidth',2);
e_peak_amp.mean = sum(polya_cnt_fit.*xbins)/sum(polya_cnt_fit);
[dummy, e_peak_amp.max_idx] = max(polya_cnt_fit);
saveas(gcf,[store_folder '\Run' run.id '_echarge_hist.png'])




%%sampling in small area

%circular cut_small - sampling in center
cut_sampling = (((trackerX - pad.xc-shiftX).^2 + (trackerY - pad.yc-shiftY).^2) < small.r^2) & glbl_cut;

%ring-shaped cut small
%cut_small = (((trackerX - pad.xc).^2 + (trackerY - pad.yc).^2) > small.r^2) & (((trackerX - pad.xc).^2 + (trackerY - pad.yc).^2) < medium.r^2) & glbl_cut;

mean(time_diff_sigmoidTWCorr(cut_sampling))
std(time_diff_sigmoidTWCorr(cut_sampling))

%fractionEventsSelectedInPad = length(time_diff_sigmoidTWCorr(cut_small))/length(time_diff_sigmoidTWCorr);
fractionEventsSelectedForTiming = length(time_diff_sigmoidTWCorr(cut_sampling))/length(time_diff_sigmoidTWCorr);



%% plot peak amplitude histogram cut included - ONLY SAMPLING AREA
figure
histogram(MCP_maxy(cut_sampling),100)
hold on
h=histogram(MM_maxy(cut_sampling),100);
hold on
xbins = h.BinEdges(1:end-1)+h.BinWidth/2;
fit_data = [];
fit_data(1,:)=1:length(xbins);
fit_data(2,:) = h.Values;
% fit_data(3,:) = yerr;
p0=[];
p0(1) = sum(h.Values)*h.BinWidth;   % normalization factor
p0(2) = 1;
p0(3) = 0.1;
cmd='min; ret';
[p, err, chi] = fminuit('polya_minuit',p0,fit_data(:,1:end),'-b','-c',cmd);
polya_cnt_fit = polya_minuit(p,fit_data(1,:));
plot(xbins,polya_cnt_fit,'Linewidth',2);
e_peak_amp.mean = sum(polya_cnt_fit.*xbins)/sum(polya_cnt_fit);
[dummy, e_peak_amp.max_idx] = max(polya_cnt_fit);

xlabel('Signal amplitude, V')
ylabel('Events');
grid on
legend('MCP','MM','MM fit');
title_str = sprintf('%s \n e-peak amplitude \\mu = %4.4f V U_{max} = %4.4f V',runTitleString, e_peak_amp.mean,xbins(e_peak_amp.max_idx));
title(title_str)
movegui(gcf,'east');

saveas(gcf,[store_folder '\Run' run.id '_signalAmplitude_Hist_comparison-samplingArea.png'])

%% plot peak amplitude histogram cut included - Only MM - ONLY SAMPLING AREA
figure
h=histogram(MM_maxy(cut_sampling),100);
hold on
xbins = h.BinEdges(1:end-1)+h.BinWidth/2;
fit_data = [];
fit_data(1,:)=1:length(xbins);
fit_data(2,:) = h.Values;
% fit_data(3,:) = yerr;
p0=[];
p0(1) = sum(h.Values)*h.BinWidth;   % normalization factor
p0(2) = 1;
p0(3) = 0.1;
cmd='min; ret';
[p, err, chi] = fminuit('polya_minuit',p0,fit_data(:,1:end),'-b','-c',cmd);
polya_cnt_fit = polya_minuit(p,fit_data(1,:));
plot(xbins,polya_cnt_fit,'Linewidth',2);
e_peak_amp.mean = sum(polya_cnt_fit.*xbins)/sum(polya_cnt_fit);
[dummy, e_peak_amp.max_idx] = max(polya_cnt_fit);

xlabel('Signal amplitude, V')
ylabel('Events');
grid on
legend('MM','MM fit');
title_str = sprintf('%s \n e-peak amplitude \\mu = %4.4f V U_{max} = %4.4f V',runTitleString, e_peak_amp.mean,xbins(e_peak_amp.max_idx));
title(title_str)
saveas(gcf,[store_folder '\Run' run.id '_signalAmplitude_Hist_MM-samplingArea.png'])




%% plot electron charge
figure;
plot(MM_maxy(cut_sampling), e_peak_MM(cut_sampling),'.');
grid on
xlabel('Signal amplitude, V')
ylabel('Electron charge, pC');
title_str = sprintf('%s \n Electron charge plot',runTitleString);
title(title_str)
saveas(gcf,[store_folder '\Run' run.id '_echarge_vs_ampl-samplingArea.png'])

%% plot electron lead charge histogram - ONLY SAMPLING AREA

figure;
h=histogram(e_lead_MM(glbl_cut),100);
hold on
xbins = h.BinEdges(1:end-1)+h.BinWidth/2;
fit_data = [];
fit_data(1,:)=1:length(xbins);
fit_data(2,:) = h.Values;
% fit_data(3,:) = yerr;
p0=[];
p0(1) = sum(h.Values)*h.BinWidth;   % normalization factor
p0(2) = 1;
p0(3) = 0.1;
cmd='min; ret';
[p, err, chi] = fminuit('polya_minuit',p0,fit_data(:,1:end),'-b','-c',cmd);
polya_cnt_fit = polya_minuit(p,fit_data(1,:));
plot(xbins,polya_cnt_fit,'Linewidth',2);
e_peak_amp.mean = sum(polya_cnt_fit.*xbins)/sum(polya_cnt_fit);
[dummy, e_peak_amp.max_idx] = max(polya_cnt_fit);

xlabel('Electron charge, pC');
ylabel('Events');
grid on
legend('Hist','Polya fit');
title_str = sprintf('%s - e-lead charge \\mu = %4.4f pC Q_{max} = %4.4f pC',runTitleString, e_peak_amp.mean,xbins(e_peak_amp.max_idx));
title(title_str)
saveas(gcf,[store_folder '\Run' run.id '_eleadcharge_hist-samplingArea.png'])


%% plot electron peak charge histogram - ONLY SAMPLING AREA

figure;
h=histogram(e_peak_MM(cut_sampling),100);
hold on
xbins = h.BinEdges(1:end-1)+h.BinWidth/2;
fit_data = [];
fit_data(1,:)=1:length(xbins);
fit_data(2,:) = h.Values;
% fit_data(3,:) = yerr;
p0=[];
p0(1) = sum(h.Values)*h.BinWidth;   % normalization factor
p0(2) = 1;
p0(3) = 0.1;
cmd='min; ret';
[p, err, chi] = fminuit('polya_minuit',p0,fit_data(:,1:end),'-b','-c',cmd);
polya_cnt_fit = polya_minuit(p,fit_data(1,:));
plot(xbins,polya_cnt_fit,'Linewidth',2);
e_peak_amp.mean = sum(polya_cnt_fit.*xbins)/sum(polya_cnt_fit);
[dummy, e_peak_amp.max_idx] = max(polya_cnt_fit);
saveas(gcf,[store_folder '\Run' run.id '_echarge_hist-samplingArea.png'])


%% timing in small area

% plot single gauss histogram
figure
%rng(abs(pd.mu));
h = histogram(time_diff_sigmoidTWCorr(cut_sampling),100);
hold on
fit_data = [];
fit_data(1,:) = h.BinEdges(1:end-1) + h.BinWidth/2;
fit_data(2,:) = h.Values;
% fit_data(3,:) = yerr;

% try different combinations factors in order to find the min chi-square

p0=[];
p0(1) = sum(fit_data(2,:))*h.BinWidth/std(h.Data)/sqrt(2*pi);   % total normalization factor
p0(2) = std(time_diff_sigmoidTWCorr(cut_sampling));                       % sigma1
p0(3) = mean(time_diff_sigmoidTWCorr(cut_sampling));                      % mean1
%   p0(4) = comb(i);                                  % combination factotr
% p0(5) = std(time_diff_sigmoidTWCorr(cut_sampling));                       % sigma2
%p0(6) = mean(fit_data(1,:));                     % mean2 (not used) same as p3
step = [1 2 3];
l_b = [p0(1)*0.8 0.1*p0(2) p0(3)-0.1];
u_b = [p0(1)*1.2 10*p0(2) p0(3)+0.1];
StepBounds= [step; l_b; u_b]';
MinuitCommands = 'min; ret;';
chi_idx = fit_data(2,:) > 0;  % take into account only non-empty bins

[p, err, chi2_min] = fminuit('gauss1_minuit', p0,fit_data(:,chi_idx), '-c', MinuitCommands,'-s', StepBounds)
%    [p, err, chi2_min] = fminuit('gauss1_minuit', p0,fit_data(:,chi_idx), '-c', MinuitCommands,'-s');
%    g2.chi2_vec(i) = chi2_min;

% find minimum chi-square
%[chi2_min,idx_min_chi] = min(g2.chi2_vec);
%p0(4) = comb(idx_min_chi);
%[p, err, chi2_min] = fminuit('gauss1_minuit', p0,fit_data(:,chi_idx), '-c', MinuitCommands,'-s', StepBounds);

%g2.chi2 = chi2_min;
%g2.ndf = sum(chi_idx) - length(p);
plot(fit_data(1,:),gauss1_minuit(p,fit_data(1,:)),'Linewidth',2);
%plot(fit_data(1,:),gauss1_minuit([p(1) p(2) p(3)],fit_data(1,:)),'Linewidth',2);
%plot(fit_data(1,:),gauss1_minuit([p(1)*(1-p(4)) p(5) p(3)],fit_data(1,:)),'Linewidth',2);
%g2.p=p(4);
%g2.mu = [p(3) p(3)];
%g2.mu_err = [err(3) err(3)];
%g2.sigma = [p(2) p(5)];
%g2.sigma_err = [err(2) err(5)];
%g2.sigma_all = sqrt(g2.p*g2.sigma(1)^2+(1-g2.p)*g2.sigma(2)^2+g2.p*(1-g2.p)*(g2.mu(1)-g2.mu(2))^2);
%g2.mu_all = g2.p*g2.mu(1)+(1-g2.p)*g2.mu(2);

message =  sprintf('  \\mu = %2.3f ns \\pm %2.3f ps\n',p(3),1000*err(3));
message = [message sprintf('  \\sigma_{fit} = %2.1f ps \\pm %2.3f ps\n',1000*p(2),1000*err(2))];
%message = [message sprintf('  \\sigma_2 = %2.1f ps \\pm %2.3f ps\n',1000*g2.sigma(2),1000*g2.sigma_err(2))];
%message = [message sprintf('  \\sigma_{tot} = %2.1f ps \n',1000*g2.sigma_all)];
message = [message sprintf('  RMS_{hist} = %2.1f ps ',1000*std(time_diff_sigmoidTWCorr(cut_sampling)))];

xlabel('Time difference: PICOSEC vs Reference, ns');
ylabel('Events');
xlim ([-0.2 0.2])
legend('RAW hist','Gauss fit');
grid
title_str = sprintf('Time resolution - Single Gauss fit');
title(title_str)
y_pos=get(gca,'ylim');
x_pos=get(gca,'xlim');
text(x_pos(1),0.88*y_pos(2),message)
movegui(gcf,'southeast');
saveas(gcf,[store_folder '\Run' run.id '_time_res-samplingArea-singleGauss.png'])



time_diff_sigmoidTWCorrMaskOutliers = time_diff_sigmoidTWCorr>0.12 | time_diff_sigmoidTWCorr<-0.12 ;
outliersMask = time_diff_sigmoidTWCorrMaskOutliers & cut_sampling;
evendIDOutliers = eventIDArray(outliersMask);
riseTimeOutliers = riseTime(outliersMask);


%% timing in small area

% plot dual gauss histogram
figure
%rng(abs(pd.mu));
h = histogram(time_diff_sigmoidTWCorr(cut_sampling),100);
hold on
fit_data = [];
fit_data(1,:) = h.BinEdges(1:end-1) + h.BinWidth/2;
fit_data(2,:) = h.Values;
% fit_data(3,:) = yerr;

% try different combinations factors in order to find the min chi-square
comb = 0.05:0.05:0.95;
for i = 1:length(comb)
    p0=[];
    p0(1) = sum(fit_data(2,:))*h.BinWidth/std(h.Data)/sqrt(2*pi);   % total normalization factor
    p0(2) = std(time_diff_sigmoidTWCorr(cut_sampling));                       % sigma1
    p0(3) = mean(time_diff_sigmoidTWCorr(cut_sampling));                      % mean1
    p0(4) = comb(i);                                  % combination factotr
    p0(5) = std(time_diff_sigmoidTWCorr(cut_sampling));                       % sigma2
    %p0(6) = mean(fit_data(1,:));                     % mean2 (not used) same as p3
    step = [1 2 3 4 5];
    l_b = [p0(1)*0.8 0.1*p0(2) p0(3)-0.1 0.5 0.1*p0(5)];
    u_b = [p0(1)*1.2 10*p0(2) p0(3)+0.1 0.999 10*p0(5)];
    StepBounds= [step; l_b; u_b]';
    MinuitCommands = 'min; ret;';
    chi_idx = fit_data(2,:) > 0;  % take into account only non-empty bins
    [p, err, chi2_min] = fminuit('gauss2_minuit', p0,fit_data(:,chi_idx), '-c', MinuitCommands,'-s', StepBounds);
    g2.chi2_vec(i) = chi2_min;
end
% find minimum chi-square
[chi2_min,idx_min_chi] = min(g2.chi2_vec);
p0(4) = comb(idx_min_chi);
[p, err, chi2_min] = fminuit('gauss2_minuit', p0,fit_data(:,chi_idx), '-c', MinuitCommands,'-s', StepBounds);

g2.chi2 = chi2_min;
g2.ndf = sum(chi_idx) - length(p);
plot(fit_data(1,:),gauss2_minuit(p,fit_data(1,:)),'Linewidth',2);
plot(fit_data(1,:),gauss1_minuit([p(1)*p(4) p(2) p(3)],fit_data(1,:)),'Linewidth',2);
plot(fit_data(1,:),gauss1_minuit([p(1)*(1-p(4)) p(5) p(3)],fit_data(1,:)),'Linewidth',2);
g2.p=p(4);
g2.mu = [p(3) p(3)];
g2.mu_err = [err(3) err(3)];
g2.sigma = [p(2) p(5)];
g2.sigma_err = [err(2) err(5)];
g2.sigma_all = sqrt(g2.p*g2.sigma(1)^2+(1-g2.p)*g2.sigma(2)^2+g2.p*(1-g2.p)*(g2.mu(1)-g2.mu(2))^2);
g2.mu_all = g2.p*g2.mu(1)+(1-g2.p)*g2.mu(2);

message = sprintf('  \\chi^2 / NDF = %2.1f / %d\n\n',g2.chi2,g2.ndf);
message = [message sprintf('  \\mu = %2.3f ns \\pm %2.3f ps\n',g2.mu_all,1000*g2.mu_err(1))];
message = [message sprintf('  \\sigma_1 = %2.1f ps \\pm %2.3f ps\n',1000*g2.sigma(1),1000*g2.sigma_err(1))];
message = [message sprintf('  \\sigma_2 = %2.1f ps \\pm %2.3f ps\n',1000*g2.sigma(2),1000*g2.sigma_err(2))];
message = [message sprintf('  \\sigma_{tot} = %2.1f ps \n',1000*g2.sigma_all)];
message = [message sprintf('  RMS_{tot} = %2.1f ps ',1000*std(time_diff_sigmoidTWCorr(cut_sampling)))];

xlabel('Time difference, ns');
ylabel('Events');
legend('RAW hist','Gauss combined','Gauss 1', 'Gauss 2');
grid
title_str = sprintf('2Gauss %s (\\phi = %2.1f mm)', runTitleString, 2*small.r );
title(title_str)
y_pos=get(gca,'ylim');
x_pos=get(gca,'xlim');
text(x_pos(1),0.75*y_pos(2),message)
movegui(gcf,'southeast');
saveas(gcf,[store_folder '\Run' run.id '_time_res-samplingArea.png'])



%% analysis of symmetry of REF MCP
sym_angles = [0 45 90 135];
figure;
for j=1:length(sym_angles)
    line_amp.angle = sym_angles(j);     % line angle
    line_amp.step = 0.5;     % spacing in mm
    line_amp.half_len = 10;   % half length of the line
    line_amp.pts = -line_amp.half_len:line_amp.step:line_amp.half_len;
    line_amp.x = cos(line_amp.angle*pi/180)*line_amp.pts + pad.xc;
    line_amp.y = sin(line_amp.angle*pi/180)*line_amp.pts + pad.yc;
    line_amp.n = length(line_amp.pts);

    line_amp.radius = 1;
    for i = 1:line_amp.n
        cut_circ = ((trackerX - line_amp.x(i)).^2 + (trackerY - line_amp.y(i)).^2) < line_amp.radius^2 & glbl_cut;
        line_amp.avg(i) = mean(MM_maxy(cut_circ));
        line_amp.n_hits(i) = sum(cut_circ);
        line_amp.err(i) = std(MM_maxy(cut_circ))./sqrt(line_amp.n_hits(i));
    end
    % plot amplitude over observation line
    line_amp.legend{j} = ['Cutting angle ' num2str(sym_angles(j)) '°'];
    errorbar(line_amp.pts,line_amp.avg,line_amp.err,'.', 'LineWidth',1);
    hold on;

end
str_title = sprintf('%s:\n DUT amplitude cross-sections ',runTitleString);
title(str_title);
legend(line_amp.legend, 'Location', 'Best');
xlabel('Distance from pad centre, mm');
ylabel('Amplitude, V');
grid on;
saveas(gcf,[store_folder '\Run' run.id '_DUT_amp_cross.png'])



figure;
subplot(3,3,1); %timewalk plot
fit_data = [];
fit_data(1,:) = twalk.e_peak;
fit_data(2,:) = twalk.mean_sat;
fit_data(3,:) = twalk.err;
p0=[];
p0(1) = min(twalk.mean_sat);
p0(2) = 1;
p0(3) = 0.5;
cmd='min; ret';
[p, err, chi] = fminuit('twalk_fn_minuit',p0,fit_data,'-b','-c',cmd);
twalk.p = p;
twalk.chi = chi;
% plot time walk correction and fit
errorbar(twalk.e_peak,twalk.mean_sat,twalk.err,twalk.err,twalk.e_peak_err_n,twalk.e_peak_err_p,'o');
hold on
plot(twalk.e_peak,twalk_fn_minuit(p,twalk.e_peak),'LineWidth',1.5);
xlabel('Electron peak charge, pC')
ylabel('SAT, ns')
title_str = sprintf('SAT vs. Electron peak charge');%,runTitleString);
title(title_str)
grid
hold off

subplot(3,3,2); % hits
scatter(trackerX, trackerY,'.')   % plot all hits
hold on;
scatter(trackerX(glbl_cut), trackerY(glbl_cut),'.'); % plot hits with global cut
scatter(pad.xc,pad.yc);                            % plot center point
rectangle('Position',[pad.xc-pad.size/2 pad.yc-pad.size/2 pad.size pad.size], 'Curvature',pad.curvature, 'LineWidth',1);    %plot active area
rectangle('Position',[pad.xc-shiftX-(2*small.r)/2 pad.yc-shiftY-(2*small.r)/2 (2*small.r) (2*small.r)], 'Curvature',[1,1], 'LineWidth',2,'EdgeColor','red');    %plot selected area
axis equal
xlabel('x-axis, mm');
ylabel('y-axis, mm');
grid on
str_title = sprintf('Tracker data');
title(str_title);
legend('Trigger','Hits','Observation path','DUT center','Circle for resolution measurement');
if length(trackerX(glbl_cut))>0
    xlim([min(trackerX(glbl_cut)) max(trackerX(glbl_cut))]);
    ylim([min(trackerY(glbl_cut)) max(trackerY(glbl_cut))]);
end


subplot(3,3,3); % SAT lines
for linePos = 1:length(lineAnglesArray)
    % plot SAT over observation line
    errorbar(satLinesX(linePos,:),satLinesY(linePos,:),satLinesErr(linePos,:),'o', 'LineWidth',1);
    str_title = sprintf('SAT Comparison');
    title(str_title);
    xlabel('Distance from pad centre, mm');
    ylabel('SAT, ps');
    grid on;

    if linePos==1
        hold on
    end
end
legend(satLines.legend);


subplot(3,3,4); % amplitude
histogram(MCP_maxy(cut_sampling),100)
hold on
h=histogram(MM_maxy(cut_sampling),100);
xbins = h.BinEdges(1:end-1)+h.BinWidth/2;
fit_data = [];
fit_data(1,:)=1:length(xbins);
fit_data(2,:) = h.Values;
% fit_data(3,:) = yerr;
p0=[];
p0(1) = sum(h.Values)*h.BinWidth;   % normalization factor
p0(2) = 1;
p0(3) = 0.1;
cmd='min; ret';
[p, err, chi] = fminuit('polya_minuit',p0,fit_data(:,1:end),'-b','-c',cmd);
polya_cnt_fit = polya_minuit(p,fit_data(1,:));
plot(xbins,polya_cnt_fit,'Linewidth',2);
e_peak_amp.mean = sum(polya_cnt_fit.*xbins)/sum(polya_cnt_fit);
[dummy, e_peak_amp.max_idx] = max(polya_cnt_fit);
xlabel('Signal amplitude, V')
ylabel('Events');
grid on
legend('MCP','MM','MM fit');
title_str = sprintf('e-peak amplitude \\mu = %4.4f, V U_{max} = %4.4f V', e_peak_amp.mean,xbins(e_peak_amp.max_idx));
title(title_str)


subplot(3,3,5); % EfficiencyDUT
h=pcolor(area.xx,area.yy,area.dutEfficiency);
hold on
rectangle('Position',[pad.xc-pad.size/2 pad.yc-pad.size/2 pad.size pad.size], 'Curvature',pad.curvature, 'LineWidth',1);
rectangle('Position',[pad.xc-shiftX-(2*small.r)/2 pad.yc-shiftY-(2*small.r)/2 (2*small.r) (2*small.r)], 'Curvature',[1,1], 'LineWidth',2,'EdgeColor','red');    %plot selected area
set(h, 'EdgeColor', 'none');
%set(gca, 'XDir','reverse')
axis equal
colorbar
xlabel('x-axis, mm');
ylabel('y-axis, mm');
h = colorbar;
h.Label.String = 'Efficiency DUT';
caxis([0 1]);
h.Label.Position(1) = 3;
str_title = sprintf('DUT Efficiency \n All Events: ActiveArea: %.1f%% | CenterSampling: %.1f%% \n Included TimeRes: ActiveArea: %.1f%% | CenterSampling: %.1f%%', efficiencyDUTActiveArea*100, efficiencyDUTSamplingArea*100, efficiencyDUTActiveAreaAccepted*100, efficiencyDUTSamplingAreaAccepted*100);
title(str_title);


subplot(3,3,6); % MapAmplitude
h=histogram(riseTime(glbl_cut&riseTimeCut),100);
hold on
xbins = h.BinEdges(1:end-1)+h.BinWidth/2;
xlabel('Rise time (ns)')
ylabel('Events');
grid on
xlim([0 2]);
title_str = sprintf('Signal rise time \\mu = %4.4f ns', mean(riseTime(riseTimeCut&glbl_cut)));
title(title_str)


subplot(3,3,7); % MapAmplitude
h = pcolor(area.xx,area.yy,area.amp_MCP);
hold on
rectangle('Position',[pad.xc-pad.size/2 pad.yc-pad.size/2 pad.size pad.size], 'Curvature',pad.curvature, 'LineWidth',2);
rectangle('Position',[pad.xc-shiftX-(2*small.r)/2 pad.yc-shiftY-(2*small.r)/2 (2*small.r) (2*small.r)], 'Curvature',[1,1], 'LineWidth',2,'EdgeColor','red');    %plot selected area
set(h, 'EdgeColor', 'none');
axis equal
colorbar
xlabel('x-axis, mm');
ylabel('y-axis, mm');
h = colorbar;
h.Label.String = 'Amplitude, V';
h.Label.Position(1) = 3;
str_title = sprintf('REF MCP amplitude map (\\phi_{avg} = %2.1f mm)', 2*area.radius);
title(str_title);


subplot(3,3,8); % MapAmplitude
h=pcolor(area.xx,area.yy,area.amp);
hold on
rectangle('Position',[pad.xc-pad.size/2 pad.yc-pad.size/2 pad.size pad.size],  'Curvature',pad.curvature,'LineWidth',2);
rectangle('Position',[pad.xc-shiftX-(2*small.r)/2 pad.yc-shiftY-(2*small.r)/2 (2*small.r) (2*small.r)], 'Curvature',[1,1], 'LineWidth',2,'EdgeColor','red');    %plot selected area
set(h, 'EdgeColor', 'none');
axis equal
colorbar
xlabel('x-axis, mm');
ylabel('y-axis, mm');
h = colorbar;
h.Label.String = 'Amplitude, V';
h.Label.Position(1) = 3;
str_title = sprintf('DUT MM map (\\phi_{avg} = %2.1f mm)', 2*area.radius);
title(str_title);


% subplot(3,3,8); % MapAmplitude
% h=surf(area.xx,area.yy,area.amp);
% hold on
% rectangle('Position',[pad.xc-pad.size/2 pad.yc-pad.size/2 pad.size pad.size],  'Curvature',pad.curvature,'LineWidth',2);
% rectangle('Position',[pad.xc-shiftX-(2*small.r)/2 pad.yc-shiftY-(2*small.r)/2 (2*small.r) (2*small.r)], 'Curvature',[1,1], 'LineWidth',2,'EdgeColor','red');    %plot selected area
% set(h, 'EdgeColor', 'none');
% %axis equal
% colorbar
% xlabel('x-axis, mm');
% ylabel('y-axis, mm');
% h = colorbar;
% h.Label.String = 'Amplitude, V';
% h.Label.Position(1) = 3;
% str_title = sprintf('DUT MM map (\\phi_{avg} = %2.1f mm)', 2*area.radius);
% title(str_title);
% 

%% timing plot subplot - either double or single gauss


% subplot(3,3,9); % Timing res
% 
% h = histogram(time_diff_sigmoidTWCorr(cut_sampling),100);
% hold on
% fit_data = [];
% fit_data(1,:) = h.BinEdges(1:end-1) + h.BinWidth/2;
% fit_data(2,:) = h.Values;
% % fit_data(3,:) = yerr;
% 
% % try different combinations factors in order to find the min chi-square
% comb = 0.05:0.05:0.95;
% for i = 1:length(comb)
%     p0=[];
%     p0(1) = sum(fit_data(2,:))*h.BinWidth/std(h.Data)/sqrt(2*pi);   % total normalization factor
%     p0(2) = std(time_diff_sigmoidTWCorr(cut_sampling));                       % sigma1
%     p0(3) = mean(time_diff_sigmoidTWCorr(cut_sampling));                      % mean1
%     p0(4) = comb(i);                                  % combination factotr
%     p0(5) = std(time_diff_sigmoidTWCorr(cut_sampling));                       % sigma2
%     %p0(6) = mean(fit_data(1,:));                     % mean2 (not used) same as p3
%     step = [1 2 3 4 5];
%     l_b = [p0(1)*0.8 0.1*p0(2) p0(3)-0.1 0.5 0.1*p0(5)];
%     u_b = [p0(1)*1.2 10*p0(2) p0(3)+0.1 0.999 10*p0(5)];
%     StepBounds= [step; l_b; u_b]';
%     MinuitCommands = 'min; ret;';
%     chi_idx = fit_data(2,:) > 0;  % take into account only non-empty bins
%     [p, err, chi2_min] = fminuit('gauss2_minuit', p0,fit_data(:,chi_idx), '-c', MinuitCommands,'-s', StepBounds);
%     g2.chi2_vec(i) = chi2_min;
% end
% % find minimum chi-square
% [chi2_min,idx_min_chi] = min(g2.chi2_vec);
% p0(4) = comb(idx_min_chi);
% [p, err, chi2_min] = fminuit('gauss2_minuit', p0,fit_data(:,chi_idx), '-c', MinuitCommands,'-s', StepBounds);
% 
% g2.chi2 = chi2_min;
% g2.ndf = sum(chi_idx) - length(p);
% plot(fit_data(1,:),gauss2_minuit(p,fit_data(1,:)),'Linewidth',2);
% plot(fit_data(1,:),gauss1_minuit([p(1)*p(4) p(2) p(3)],fit_data(1,:)),'Linewidth',2);
% plot(fit_data(1,:),gauss1_minuit([p(1)*(1-p(4)) p(5) p(3)],fit_data(1,:)),'Linewidth',2);
% g2.p=p(4);
% g2.mu = [p(3) p(3)];
% g2.mu_err = [err(3) err(3)];
% g2.sigma = [p(2) p(5)];
% g2.sigma_err = [err(2) err(5)];
% g2.sigma_all = sqrt(g2.p*g2.sigma(1)^2+(1-g2.p)*g2.sigma(2)^2+g2.p*(1-g2.p)*(g2.mu(1)-g2.mu(2))^2);
% g2.mu_all = g2.p*g2.mu(1)+(1-g2.p)*g2.mu(2);
% 
% message = sprintf('  \\chi^2 / NDF = %2.1f / %d\n\n',g2.chi2,g2.ndf);
% message = [message sprintf('  \\mu = %2.3f ns \\pm %2.3f ps\n',g2.mu_all,1000*g2.mu_err(1))];
% message = [message sprintf('  \\sigma_1 = %2.1f ps \\pm %2.3f ps\n',1000*g2.sigma(1),1000*g2.sigma_err(1))];
% message = [message sprintf('  \\sigma_2 = %2.1f ps \\pm %2.3f ps\n',1000*g2.sigma(2),1000*g2.sigma_err(2))];
% message = [message sprintf('  \\sigma_{tot} = %2.1f ps \n',1000*g2.sigma_all)];
% message = [message sprintf('  RMS_{tot} = %2.1f ps ',1000*std(time_diff_sigmoidTWCorr(cut_sampling)))];
% 
% xlabel('Time difference, ns');
% ylabel('events');
% legend('RAW hist','Gauss combined','Gauss 1', 'Gauss 2');
% grid
% title_str = sprintf('2Gauss (\\phi = %2.1f mm)', 2*small.r );
% title(title_str)
% y_pos=get(gca,'ylim');
% x_pos=get(gca,'xlim');
% text(x_pos(1),0.75*y_pos(2),message)
% 
% title_str = sprintf('%s',runTitleString );
% 
% 
% defect = time_diff_sigmoidTWCorr < -0.2 & cut_sampling;
% time_diff_sigmoidTWCorr(defect) = -0.21;

% plot single gauss histogram
subplot(3,3,9); % Timing res
h = histogram(time_diff_sigmoidTWCorr(cut_sampling),100);
hold on
fit_data = [];
fit_data(1,:) = h.BinEdges(1:end-1) + h.BinWidth/2;
fit_data(2,:) = h.Values;
% fit_data(3,:) = yerr;



% try different combinations factors in order to find the min chi-square

p0=[];
p0(1) = sum(fit_data(2,:))*h.BinWidth/std(h.Data)/sqrt(2*pi);   % total normalization factor
p0(2) = std(time_diff_sigmoidTWCorr(cut_sampling));                       % sigma1
p0(3) = mean(time_diff_sigmoidTWCorr(cut_sampling));                      % mean1
%   p0(4) = comb(i);                                  % combination factotr
% p0(5) = std(time_diff_sigmoidTWCorr(cut_sampling));                       % sigma2
%p0(6) = mean(fit_data(1,:));                     % mean2 (not used) same as p3
step = [1 2 3];
l_b = [p0(1)*0.8 0.1*p0(2) p0(3)-0.1];
u_b = [p0(1)*1.2 10*p0(2) p0(3)+0.1];
StepBounds= [step; l_b; u_b]';
MinuitCommands = 'min; ret;';
chi_idx = fit_data(2,:) > 0;  % take into account only non-empty bins

[p, err, chi2_min] = fminuit('gauss1_minuit', p0,fit_data(:,chi_idx), '-c', MinuitCommands,'-s', StepBounds)
%    [p, err, chi2_min] = fminuit('gauss1_minuit', p0,fit_data(:,chi_idx), '-c', MinuitCommands,'-s');
%    g2.chi2_vec(i) = chi2_min;

% find minimum chi-square
%[chi2_min,idx_min_chi] = min(g2.chi2_vec);
%p0(4) = comb(idx_min_chi);
%[p, err, chi2_min] = fminuit('gauss1_minuit', p0,fit_data(:,chi_idx), '-c', MinuitCommands,'-s', StepBounds);

%g2.chi2 = chi2_min;
%g2.ndf = sum(chi_idx) - length(p);
plot(fit_data(1,:),gauss1_minuit(p,fit_data(1,:)),'Linewidth',2);
%plot(fit_data(1,:),gauss1_minuit([p(1) p(2) p(3)],fit_data(1,:)),'Linewidth',2);
%plot(fit_data(1,:),gauss1_minuit([p(1)*(1-p(4)) p(5) p(3)],fit_data(1,:)),'Linewidth',2);
%g2.p=p(4);
%g2.mu = [p(3) p(3)];
%g2.mu_err = [err(3) err(3)];
%g2.sigma = [p(2) p(5)];
%g2.sigma_err = [err(2) err(5)];
%g2.sigma_all = sqrt(g2.p*g2.sigma(1)^2+(1-g2.p)*g2.sigma(2)^2+g2.p*(1-g2.p)*(g2.mu(1)-g2.mu(2))^2);
%g2.mu_all = g2.p*g2.mu(1)+(1-g2.p)*g2.mu(2);

message =  sprintf('  \\mu = %2.3f ns \\pm %2.3f ps\n',p(3),1000*err(3));
message = [message sprintf('  \\sigma_{fit} = %2.1f ps \\pm %2.3f ps\n',1000*p(2),1000*err(2))];
%message = [message sprintf('  \\sigma_2 = %2.1f ps \\pm %2.3f ps\n',1000*g2.sigma(2),1000*g2.sigma_err(2))];
%message = [message sprintf('  \\sigma_{tot} = %2.1f ps \n',1000*g2.sigma_all)];
message = [message sprintf('  RMS_{hist} = %2.1f ps ',1000*std(time_diff_sigmoidTWCorr(cut_sampling)))];

xlabel('Time difference: PICOSEC vs Reference, ns');
ylabel('Events');
xlim ([-0.2 0.2])
legend('RAW hist','Gauss fit');
grid
title_str = sprintf('Time resolution - Single Gauss fit');
title(title_str)
y_pos=get(gca,'ylim');
x_pos=get(gca,'xlim');
text(x_pos(1),0.88*y_pos(2),message)


title_str = sprintf('%s',runTitleString );
sgtitle(title_str)
%text(5,50, title_str)

set(gcf, 'Position', get(0, 'Screensize'));
saveas(gcf,[store_folder '\Run' run.id '.png'])


%% save main result numbers in txt file
amplitudeMean =  e_peak_amp.mean
meanRiseTime = mean(riseTime(riseTimeCut&glbl_cut));
rmsTimeRes = 1000*std(time_diff_sigmoidTWCorr(cut_sampling));
efficiencySamplingArea = efficiencyDUTSamplingAreaAccepted*100;

 txt_tab.names = {'RunID','RunDesc','amplitudeMean','meanRiseTime','rmsTimeRes', 'efficiencySamplingArea'};
 fid = fopen([store_folder '/Run' run.id '_results.txt'],'w');
 % print header with names
for i=1:length(txt_tab.names)
     fprintf(fid,'%s\t',txt_tab.names{i});
 end
 fprintf(fid,'\n');

 fprintf(fid,'%d\t',runID);
 fprintf(fid,'%s\t',runTitleString);
 fprintf(fid,'%.4f\t',efficiencySamplingArea);
 fprintf(fid,'%.4f\t',meanRiseTime);
 fprintf(fid,'%.4f\t',amplitudeMean);
 fprintf(fid,'%.4f\t',rmsTimeRes);
 fprintf(fid,'%.4f\t',1000*p(2)); %timing from single gauss
 

 fclose(fid);

 
 fid = fopen(resultsTablePath,'a+');
 fprintf(fid,'\n');
 fprintf(fid,'%d\t',runID);
 fprintf(fid,'%.4f\t',efficiencySamplingArea);
 fprintf(fid,'%.4f\t',meanRiseTime);
 fprintf(fid,'%.4f\t',amplitudeMean);
 fprintf(fid,'%.4f\t',rmsTimeRes);
 fprintf(fid,'%.4f\t',1000*p(2)); %timing from single gauss
 fprintf(fid,'%s\t',runTitleString);
 fclose(fid);


 %% For Marta 

SATlist = time_diff_sigmoidTWCorr(cut_sampling); %for further analysis in ROOT
writematrix(SATlist(:), [store_folder '/Run' run.id '-' runInfoString '-SATlist.txt']);

%% save data to mat file and store all m files used for processing

%save([store_folder '\Run' run.id '_' run.oscilloscope '_plots.mat'], 'run', 'pad', 'line', 'area', 'small', 'g2', 'twalk'); % save interesting data for plots

% %% store tabular data
% txt_tab.names = {'No','EventID','REF_MCP_timestamp','DUT_MCP_timestamp','REF_MCP_maxY', 'DUT_MCP_maxY','REF_MCP_e_peak', 'DUT_MCP_e_peak', 'TrackerX', 'TrackerY'};
% fid = fopen([store_folder '/Run' run.id '_tabular.txt'],'w');
% % print header with names
% for i=1:length(txt_tab.names)
%     fprintf(fid,'%s\t',txt_tab.names{i});
% end
% fprintf(fid,'\n');
% for k=1:length(MCP_data)
%     fprintf(fid,'%d\t%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n', k, MM_data(k).event_id, MCP_data(k).sigmoid.timepoint, MM_data(k).sigmoid.timepoint, MCP_data(k).sig.max.y,...
%         MM_data(k).sig.max.y, MCP_data(k).sig.charge.e_peak, MM_data(k).sig.charge.e_peak, trackerX(k),trackerY(k));
% end
% fclose(fid);


%%upload results to EOS
%
% if exist('batchProcess','var') == 1
%     uploadCommand = append('pscp -r -pw ','Win_Admin',' C:\Users\GDD\Documents\MATLAB\Picosec\Analysis\analyzed\Run',run.id,'-',run.oscilloscope,'-ref',MCPnameString,'-dut',DUTnameString,' ','gdd','@lxplus.cern.ch:/eos/project/p/picosec/testbeam/2022_May_h4/Results');
%     uploadResult = system(convertStringsToChars(uploadCommand))
%
% else
%     uploadCommand = append('pscp -r -pw ','Win_Admin',' C:\Users\GDD\Documents\MATLAB\Picosec\Analysis\analyzed\Run',run.id,'-',run.oscilloscope,'-ref',opts_MCP.ch_name,'-dut',opts_MM.ch_name,' ','gdd','@lxplus.cern.ch:/eos/project/p/picosec/testbeam/2022_May_h4/Results');
%     uploadResult = system(convertStringsToChars(uploadCommand))
%
% end


%zip([store_folder '/Run' run.id '_' run.oscilloscope '_mfiles.zip'],{dir('*.m').name}); % save all m files used for procesing as a zip file
%zip([store_folder '/Run' run.id '_' run.oscilloscope '_input.zip'],in_file); % ssave input file as zip
if shouldSaveMat
    save([store_folder '\Run' run.id '-' run.oscilloscope '.mat'], 'run', 'time_diff_sigmoid', 'time_diff_sigmoidTWCorr', 'MCP_maxy', 'MM_maxy', 'trackerX', 'trackerY', 'e_peak_MM', 'eventIDArray');
    time_diff_sigmoid_samplingCut = time_diff_sigmoid(cut_sampling);
    time_diff_sigmoidTWCorr_samplingCut = time_diff_sigmoidTWCorr(cut_sampling);
    MCP_maxy_samplingCut = MCP_maxy(cut_sampling);
    MM_maxy_samplingCut = MM_maxy(cut_sampling);
    trackerX_samplingCut = trackerX(cut_sampling);
    trackerY_samplingCut = trackerY(cut_sampling);
    e_peak_MM_samplingCut = e_peak_MM(cut_sampling);
    eventIDArray_samplingCut = eventIDArray(cut_sampling);
    save([store_folder '\Run' run.id '-' run.oscilloscope 'SamplingCut.mat'], 'run', 'time_diff_sigmoid_samplingCut', 'time_diff_sigmoidTWCorr_samplingCut', 'MCP_maxy_samplingCut', 'MM_maxy_samplingCut', 'trackerX_samplingCut', 'trackerY_samplingCut', 'e_peak_MM_samplingCut', 'eventIDArray_samplingCut');
end