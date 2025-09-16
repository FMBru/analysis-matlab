%clear all
close all

addpath 'C:\Users\GDD\Documents\MATLAB\Picosec\SPEAnalysisPolyaMeanFinalVersion';


%load input file
%in_file = 'C:\Users\GDD\Documents\Picosec\May22\Analysed\Run0123-GDD.mat';

shouldSaveMat = false;

%%configuration
twalk.en = 1;                                   % enable timewalk correction


% center resolution circular - radius
small.r = 2.0; %used for only circle
medium.r = 2.5; %used for ring only

%shift sampling point for time res relative to pad center determined by
%alignmnet
%account for misalignment
shiftX = 0;
shiftY = 0;

alignToDUT = true; %if false, alignment to REF MCP

%use fit of amplitude to center sampling area
shouldCenterPadWithFit = true;

%use with large trigger to calculate correct arrival time median for global
%cut - for larger trigger or scanning runs
shouldDetermineSATMedianWithGeoCut = false;


%pad size (visualisation)
pad.size = 10;     % MCP size in mm - Hamamatsu
%pad.size = 12;     % MCP size in mm - Planacon, 2x2 central pads
pad.cut_len = 8;   % move +/-8mm from median to find center
pad.isCircle = 0;   %plotting round if circle, otherwise square

%May22: DUT index: 2: MM3 (VacChamber), 3: MM1(Multipad), 4: MM2 (on support plate), 5: MM4 (electron setup)
% if tracker.dutIndex == 2
%     pad.isCircle = 1;
% elseif tracker.dutIndex == 3
%     pad.isCircle = 0;
% elseif tracker.dutIndex == 4
%     pad.isCircle = 1;
% elseif tracker.dutIndex == 5
%     pad.isCircle = 1;
% end 


%July22: DUT index: 1: MM1 (Multipad), 2: MM2(VacChamber), 3: MM3
%(ElectronsetupMM), 4: MM4 (Picolarge main tracker), 5:MM5 (SiPM)
if tracker.dutIndex == 1
    pad.isCircle = 0;
elseif tracker.dutIndex == 2
    pad.isCircle = 1;
elseif tracker.dutIndex == 3
    pad.isCircle = 1;
elseif tracker.dutIndex == 4
    pad.isCircle = 1;
elseif tracker.dutIndex == 5
    pad.isCircle = 1;
end 

%for planacon
%pad.size = 25;     % MCP size in mm - Hamamatsu
%pad.cut_len = 25;   % move +/-8mm from median to find center
%pad.isCircle = 0;   %plotting round if circle, otherwise square


%sampling area for 2D maps
%% calculate
area.step = 0.25;   % set grid resolution in mm
area.size = 8;      % set half-size of the observed square area, mm - for
%pads and small MCP
%area.size = 20;      % set half-size of the observed square area, mm - for large MCP
area.radius = 1;    % set radius of averaging circular window cut, mm


%load(in_file);

runID = str2num(run.id);

runTitleString = [run.name ' ' run.oscilloscope ' DUT:C' opts_MM.chID ' REF:C' opts_MCP.chID ' '];

addpath 'C:\Users\GDD\Documents\MATLAB\Picosec\SCP'
connObj  =  ssh2_config('lxplus.cern.ch','gdd','Win_Admin',22);
%download matlab run data file
dataFileLocal='C:\Users\GDD\Documents\Picosec\May22\MatlabRunData';
dataFilePath = append('/eos/project/p/picosec/testbeam/2022_July_h4/info/MatlabRunInfo');
runDataFileName=append('RunData',string(runID),'.mat');
dataFileLocalPath = append(dataFileLocal,'\',runDataFileName);
% ssh2_struct = scp_get(connObj, convertStringsToChars(runDataFileName), convertStringsToChars(dataFileLocal), convertStringsToChars(dataFilePath))


%load(dataFileLocalPath);


%finish config
pad.curvature = [1,1]; %circle
if pad.isCircle == 0
    pad.curvature = [0,0]; %rect
end






%% define folder for storing data
%store_folder = ['\\eosproject-smb\eos\project\p\picosec\testbeam\2022_May_h4\Results\Run' run.id '-' run.oscilloscope '\CH' num2str(analysis.dutChannel) '-Pad' num2str(getPadForChannelNumber(analysis.dutChannel))];
%with CH folder
%store_folder = ['\\eosproject-smb\eos\project\p\picosec\testbeam\2022_July_h4\Results\Run' run.id '-' run.oscilloscope '\CH' num2str(analysis.dutChannel)];
%store_folder = ['\\eosproject-smb\eos\project\p\picosec\testbeam\2022_July_h4\Results\Run' run.id '-' run.oscilloscope '\CH' num2str(analysis.dutChannel)];


%no CH folder
store_folder = ['\\eosproject-smb\eos\project\p\picosec\testbeam\2022_July_h4\Results\Run' run.id '-' run.oscilloscope '-dut' opts_MM.chID '-ref' opts_MCP.chID];

% store_folder = ['analyzed/Run' run.id '-' run.oscilloscope];
% if exist('batchProcess','var') == 1
%     store_folder = ['analyzed/Run',run.id,'-',run.oscilloscope,'-ref',MCPnameString,'-dut',DUTnameString];
% else
%     store_folder = ['analyzed/Run',run.id,'-',run.oscilloscope,'-ref',opts_MCP.ch_name,'-dut',opts_MM.ch_name];
%     
% end


mkdir(store_folder);

%% extract from structurearrays
k=1;
for i=1:length(MCP_data)
    time_MCP(i) = MCP_data(i).sigmoid.timepoint;
    time_MM(i) = MM_data(i).sigmoid.timepoint;
    blavg_MM(i) = MM_data(i).sig.blavg;   % baseline average
    blavg_MCP(i) = MCP_data(i).sig.blavg; % baseline average MCP
    blrms_MM(i) = MM_data(i).sig.blrms;   % baseline RMS MM (noise)
    blrms_MCP(i) = MCP_data(i).sig.blrms; % baseline RMS MCP (noise)
    ymax_MCP(i)=MCP_data(i).sig.max.y;    % take maximum of MCP
    ymax_MM(i)=MM_data(i).sig.max.y;      % take maximum of MM
    ymin_MM(i)=MM_data(i).sig.min.y;      % take minimum of MM
    e_peak_MM(i)=MM_data(i).sig.charge.e_peak;
    e_peak_MCP(i)=MCP_data(i).sig.charge.e_peak;% extract electrom peak charge
    riseTime(i)=MM_data(i).sigmoid.timepoint90-MM_data(i).sigmoid.timepoint10;% extract electrom peak charge
    %if(exitflag_MCP(i)+exitflag_MM(i)==0)
end

%%used only for electron run - match with scinitllators on other scope
%run through all files and implement new array
%scintillatorAccepted =  zeros(length(time_MCP),1);

%% do some initial resolution measurement
% make simple cut with respect to median time and minimum amplitude magic
% numbers :-)
time_avg_raw = median(time_diff_sigmoid);  % assume mu from the median vaule
if shouldDetermineSATMedianWithGeoCut
   AnalyseRun_DetermineSATMedian
else 
    time_avg_raw = median(time_diff_sigmoid);  % assume mu from the median vaule
end
%time_avg_raw = 12.5;

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
    MM_maxy<0.95*max(MM_maxy)& MCP_maxy<0.95*max(MCP_maxy) & ...
    trackerX~=0 & trackerY~=0;

%additional cuts
%glbl_cut = glbl_cut &MCP_maxy>0.35;

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
% for i=1:length(twalk.epeak_vec)-1
%     temp_cut = e_peak_MM > twalk.epeak_vec(i) & e_peak_MM < twalk.epeak_vec(i+1) & glbl_cut;
%     twalk.mean_sat(i) = mean(time_diff_sigmoid(temp_cut));
%     twalk.e_peak(i) =  mean(e_peak_MM(temp_cut));
%     twalk.npts(i) = sum(temp_cut);
%     twalk.rms(i) = std(time_diff_sigmoid(temp_cut));
%     twalk.err(i) = std(time_diff_sigmoid(temp_cut))./sqrt(twalk.npts(i)); % mean error
%     twalk.e_peak_err_p(i) = twalk.e_peak(i)-(twalk.epeak_vec(i+1)); % e charge limits
%     twalk.e_peak_err_n(i) = -twalk.e_peak(i)+(twalk.epeak_vec(i));
% end
% 
% % fit correction function using minuit
% fit_data = [];
% fit_data(1,:) = twalk.e_peak;
% fit_data(2,:) = twalk.mean_sat;
% fit_data(3,:) = twalk.err;
% p0=[];
% p0(1) = min(twalk.mean_sat);
% p0(2) = 1;
% p0(3) = 0.5;
% cmd='min; ret';
% [p, err, chi] = fminuit('twalk_fn_minuit',p0,fit_data,'-b','-c',cmd);
% twalk.p = p;
% twalk.chi = chi;
% 
% % plot time walk correction and fit
% figure
% errorbar(twalk.e_peak,twalk.mean_sat,twalk.err,twalk.err,twalk.e_peak_err_n,twalk.e_peak_err_p,'o');
% hold on
% plot(twalk.e_peak,twalk_fn_minuit(p,twalk.e_peak),'LineWidth',1.5);
% xlabel('Electron peak charge, pC')
% ylabel('SAT, ns')
%             title_str = sprintf('%s \n SAT vs. Electron peak charge',runTitleString);
%     title(title_str)
% 
% grid
% movegui(gcf,'northwest');
% saveas(gcf,[store_folder '\Run' run.id '_timewalk.png'])

% plot resolution vs. e-charge
% figure
% errorbar(twalk.e_peak,twalk.rms*1000,[],[],twalk.e_peak_err_n,twalk.e_peak_err_p,'o');
% hold on
% xlabel('Electron peak charge, pC')
% ylabel('Resolution, ps')
%             title_str = sprintf('%s \n Resolution vs. Electron peak charge',runTitleString);
%     title(title_str)
% 
% grid
% saveas(gcf,[store_folder '\Run' run.id '_res_vs_charge.png'])
% 
% % make time walk correction
% if(twalk.en == 1)
%     time_diff_sigmoid = time_diff_sigmoid - twalk_fn_minuit(twalk.p, e_peak_MM);
% end


%% calculate DUT  centre if alignToDUT == false

if alignToDUT == true
    
    % calculate pad center
    pad.xc_med = median(trackerX(glbl_cut));
    pad.yc_med = median(trackerY(glbl_cut));
    
    %alternative naive mean with cut on MM pad for scanning runs
    pad.x_idx_cut = trackerX > (pad.xc_med-pad.cut_len) & trackerX < (pad.xc_med+pad.cut_len) & glbl_cut& MM_maxy>0.5*max(MM_maxy);
    pad.xc_n = mean(trackerX(pad.x_idx_cut)); % naive mean for x
    pad.y_idx_cut = trackerY > (pad.yc_med-pad.cut_len) & trackerY < (pad.yc_med+pad.cut_len) & glbl_cut& MM_maxy>0.5*max(MM_maxy);
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

%override pad centering
   % pad.xc=26;
   % pad.yc=26;


lineAnglesArray = [0,45,90];




% plot hits
figure;
scatter(trackerX, trackerY,'.')   % plot all hits
hold on;
scatter(trackerX(glbl_cut), trackerY(glbl_cut),'.'); % plot hits with global cut

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


scatter(pad.xc,pad.yc);                            % plot center point
rectangle('Position',[pad.xc-pad.size/2 pad.yc-pad.size/2 pad.size pad.size], 'Curvature',pad.curvature, 'LineWidth',1);    %plot active area
small.r = 2;
rectangle('Position',[pad.xc-shiftX-(2*small.r)/2 pad.yc-shiftY-(2*small.r)/2 (2*small.r) (2*small.r)], 'Curvature',[1,1], 'LineWidth',2,'EdgeColor','red');    %plot selected area
axis equal
xlabel('x-axis, mm');
ylabel('y-axis, mm');
grid on
str_title = sprintf('%s \n tracker data', runTitleString);
title(str_title);
legend('Trigger','Hits','Observation path','DUT center','Circle for resolution measurement');
xlim([min(trackerX(glbl_cut)) max(trackerX(glbl_cut))]);
ylim([min(trackerY(glbl_cut)) max(trackerY(glbl_cut))]);

movegui(gcf,'north');
saveas(gcf,[store_folder '\Run' run.id '_hits_overlay.png'])



%%loop over different angles and produce plots at 0/45/90deg

lineAnglesArray = [0,45,90];
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

% calculate resolution over the line
line.radius = 1;
for i = 1:line.n
    cut_circ = ((trackerX - line.x(i)).^2 + (trackerY - line.y(i)).^2) < line.radius^2 & glbl_cut;
    line.sat(i) = mean(time_diff_sigmoid(cut_circ));
    line.sat_pts(i) = sum(cut_circ);
    line.sat_err(i) = std(time_diff_sigmoid(cut_circ))./sqrt(line.sat_pts(i));
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
        area.sat(i,j) = mean(time_diff_sigmoid(cut_circ));
        area.rms(i,j) = std(time_diff_sigmoid(cut_circ));
        area.amp(i,j) = mean(ymax_MM(cut_circ));
        area.amp_MCP(i,j) = mean(ymax_MCP(cut_circ));
        area.e_peak(i,j) = mean(e_peak_MM(cut_circ));
        area.sat_pts(i,j) = sum(cut_circ);
        area.sat_err(i,j) = std(time_diff_sigmoid(cut_circ))./sqrt(area.sat_pts(i,j));
    end
end

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
caxis([5 30]);
h.Label.String = 'Time resolution, ps';
h.Label.Position(1) = 3;
str_title = sprintf('%s:\n Time resolution map (\\phi_{avg} = %2.1f mm)', runTitleString, 2*area.radius);
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
ylabel('events');
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
ylabel('events');
grid on
legend('MM','MM fit');
title_str = sprintf('%s \n e-peak amplitude \\mu = %4.4f V U_{max} = %4.4f V',runTitleString, e_peak_amp.mean,xbins(e_peak_amp.max_idx));
title(title_str)
saveas(gcf,[store_folder '\Run' run.id '_signalAmplitude_Hist_MM.png'])




%% plot rise time distribution
figure
h=histogram(riseTime(glbl_cut),100);
hold on
xbins = h.BinEdges(1:end-1)+h.BinWidth/2;

xlabel('Rise time (ns)')
ylabel('Events');
grid on
xlim([0 3]);
%legend('MM','MM fit');
title_str = sprintf('%s \n Signal rise time \\mu = %4.4f ns',runTitleString, mean(riseTime(glbl_cut)));
title(title_str)
saveas(gcf,[store_folder '\Run' run.id '_riseTimeHist.png'])





%% plot electron charge
for i=1:length(MM_data)
MM_e_charge(i)=MM_data(i).sig.charge.e_peak;
end
figure;
plot(MM_maxy(glbl_cut), MM_e_charge(glbl_cut),'.');
grid on
xlabel('Signal amplitude, V')
ylabel('Electron charge, pC');
title_str = sprintf('%s \n Electron charge plot',runTitleString);
title(title_str)
saveas(gcf,[store_folder '\Run' run.id '_echarge_vs_ampl.png'])

%% plot electron lead charge histogram
MM_e_lead_charge = zeros(length(MM_data),1);
for i=1:length(MM_data)
MM_e_lead_charge(i)=MM_data(i).sig.charge.lead_edge;
end
figure;
h=histogram(MM_e_lead_charge(glbl_cut),100);
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
ylabel('events');
grid on
legend('Hist','Polya fit');
title_str = sprintf('%s \n e-lead charge \\mu = %4.4f pC Q_{max} = %4.4f pC',runTitleString, e_peak_amp.mean,xbins(e_peak_amp.max_idx));
title(title_str)
saveas(gcf,[store_folder '\Run' run.id '_eleadcharge_hist.png'])


%% plot electron peak charge histogram
MM_e_peak_charge = zeros(length(MM_data),1);
for i=1:length(MM_data)
MM_e_peak_charge(i)=MM_data(i).sig.charge.e_peak;
end
figure;
h=histogram(MM_e_peak_charge(glbl_cut),100);
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

mean(time_diff_sigmoid(cut_sampling))
std(time_diff_sigmoid(cut_sampling))

%fractionEventsSelectedInPad = length(time_diff_sigmoid(cut_small))/length(time_diff_sigmoid);
fractionEventsSelectedForTiming = length(time_diff_sigmoid(cut_sampling))/length(time_diff_sigmoid);



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
ylabel('events');
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
ylabel('events');
grid on
legend('MM','MM fit');
title_str = sprintf('%s \n e-peak amplitude \\mu = %4.4f V U_{max} = %4.4f V',runTitleString, e_peak_amp.mean,xbins(e_peak_amp.max_idx));
title(title_str)
saveas(gcf,[store_folder '\Run' run.id '_signalAmplitude_Hist_MM-samplingArea.png'])




%% plot electron charge
for i=1:length(MM_data)
MM_e_charge(i)=MM_data(i).sig.charge.e_peak;
end
figure;
plot(MM_maxy(cut_sampling), MM_e_charge(cut_sampling),'.');
grid on
xlabel('Signal amplitude, V')
ylabel('Electron charge, pC');
title_str = sprintf('%s \n Electron charge plot',runTitleString);
title(title_str)
saveas(gcf,[store_folder '\Run' run.id '_echarge_vs_ampl-samplingArea.png'])

%% plot electron lead charge histogram - ONLY SAMPLING AREA
MM_e_lead_charge = zeros(length(MM_data),1);
for i=1:length(MM_data)
MM_e_lead_charge(i)=MM_data(i).sig.charge.lead_edge;
end
figure;
h=histogram(MM_e_lead_charge(glbl_cut),100);
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
ylabel('events');
grid on
legend('Hist','Polya fit');
title_str = sprintf('%s - e-lead charge \\mu = %4.4f pC Q_{max} = %4.4f pC',runTitleString, e_peak_amp.mean,xbins(e_peak_amp.max_idx));
title(title_str)
saveas(gcf,[store_folder '\Run' run.id '_eleadcharge_hist-samplingArea.png'])


%% plot electron peak charge histogram - ONLY SAMPLING AREA
MM_e_peak_charge = zeros(length(MM_data),1);
for i=1:length(MM_data)
MM_e_peak_charge(i)=MM_data(i).sig.charge.e_peak;
end
figure;
h=histogram(MM_e_peak_charge(cut_sampling),100);
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
h = histogram(time_diff_sigmoid(cut_sampling),100);
hold on
fit_data = [];
fit_data(1,:) = h.BinEdges(1:end-1) + h.BinWidth/2;
fit_data(2,:) = h.Values;
% fit_data(3,:) = yerr;

% try different combinations factors in order to find the min chi-square

    p0=[];
    p0(1) = sum(fit_data(2,:))*h.BinWidth/std(h.Data)/sqrt(2*pi);   % total normalization factor
    p0(2) = std(time_diff_sigmoid(cut_sampling));                       % sigma1
    p0(3) = mean(time_diff_sigmoid(cut_sampling));                      % mean1
%   p0(4) = comb(i);                                  % combination factotr
   % p0(5) = std(time_diff_sigmoid(cut_sampling));                       % sigma2
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
message = [message sprintf('  \\sigma_{hist} = %2.1f ps ',1000*std(time_diff_sigmoid(cut_sampling)))];

xlabel('Time difference, ns');
ylabel('events');
%xlim ([-17.68 -17.54])
legend('RAW hist','Gauss fit');
grid
%title_str = sprintf('2Gauss %s (\\phi = %2.1f mm)', runTitleString, 2*small.r );
%title(title_str)
y_pos=get(gca,'ylim');
x_pos=get(gca,'xlim');
text(x_pos(1),0.75*y_pos(2),message)
movegui(gcf,'southeast');
saveas(gcf,[store_folder '\Run' run.id '_time_res-samplingArea-singleGauss.png'])




%% timing in small area

% plot dual gauss histogram
figure
%rng(abs(pd.mu));
h = histogram(time_diff_sigmoid(cut_sampling),100);
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
    p0(2) = std(time_diff_sigmoid(cut_sampling));                       % sigma1
    p0(3) = mean(time_diff_sigmoid(cut_sampling));                      % mean1
    p0(4) = comb(i);                                  % combination factotr
    p0(5) = std(time_diff_sigmoid(cut_sampling));                       % sigma2
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
message = [message sprintf('  RMS_{tot} = %2.1f ps ',1000*std(time_diff_sigmoid(cut_sampling)))];

xlabel('Time difference, ns');
ylabel('events');
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
        line_amp.avg(i) = mean(ymax_MM(cut_circ));
        line_amp.n_hits(i) = sum(cut_circ);
        line_amp.err(i) = std(ymax_MM(cut_circ))./sqrt(line_amp.n_hits(i));
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


%%%
mm_max_y = [];
noise_mm_max_y = [];
int_vector = [];

for (ff=0:run.nfiles)
    
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
         
       int_iEventYvec = sum(iEventYvec); % sum of the amps in the peak. Histogram of the integral of the peak?
       int_vector= [int_vector; int_iEventYvec];
       mm_max_y = [mm_max_y; max_y]; %concatenate onto mm_max_y for each file
       
       close all

              
    end
      
    end 
       
    








%% histogram and polya fit
format long; % print 15 decimals intead of 4delat

%% Parameters for chaning

% cutoffs
%min_v = 2.29e-3;
min_v = 2.0e-3;
max_v = 0.150;

% bins
n_bins = 20;

% polya starting parameters for fitting
%x0 = [10 2 0.01]; 
x0 = [10 2 0.01];

%% Plotting
%glbl_cut = mm_max_y>0.00 & mm_max_y<0.95*max(mm_max_y); %remove saturated datapoints
glbl_cut = time_diff_sigmoid>time_min & time_diff_sigmoid<time_max & ...
    MCP_maxy>0.01*max(MCP_maxy) & MM_maxy>0.01*max(MM_maxy)&...
    MM_maxy<0.95*max(MM_maxy)& MCP_maxy<0.95*max(MCP_maxy) & ...
    trackerX~=0 & trackerY~=0;

noise_glbl_cut = noise_mm_max_y>0.00 & noise_mm_max_y<0.99*max(noise_mm_max_y);


close all;

num_bins = 500;
figure(1)
x_fig=10;
y_fig=10;
width=1500;
height=700;
set(gcf,'position',[x_fig,y_fig,width,height]);
h=histogram(mm_max_y(glbl_cut),num_bins)
set(gca, 'YScale', 'log') %plot on log scale
hold on
xlabel('Signal amplitude, V');
ylabel('Events');
ylim([0.5 max(h.Values)+1000]);
ax = gca;
ax.FontSize = 20;
%xlim([-0.5 1.5]);
grid on
%title_str = sprintf('Picosec LED measurement - Run %s - max e-peak amplitude',run.id);
title_str = sprintf('Picosec Beam measurement - Run %s - max e-peak amplitude',run.id);
title(title_str)

% set the upper and lower bounds for the lower cutoff of the polya fit that
% we would like to try
min_lower = 3e-3;
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
% min_v
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
min_arr = h.BinEdges(min_cut_idx:min_cut_idx+n_bins);

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
colors = distinguishable_colors(n_bins+1);
c = 1;
P = ["N","theta","nBar"];

for i=1:length(min_arr)
    cut_idx = min_cut_idx+i-1; 
    % define polya
    fitfun = fittype( @(N,theta,nBar,x) (N./nBar).*((theta+1).^(theta+1)).*((x./nBar).^theta).*exp(-(theta+1).*x./nBar)./gamma(theta+1));

    try
        [fitted_curve,gof] = fit(h.BinEdges(cut_idx:cutUp_idx)',h.Values(cut_idx:cutUp_idx)',fitfun,'StartPoint',x0)
    catch exception
        exceptionCounter = exceptionCounter+1;
        mean_arr(i) = 0;
        chi2_arr(i) = 0;
        mean_err_arr(i,:) = [0;0];
        continue
    end
    
    P = [P;fitted_curve.N,fitted_curve.theta,fitted_curve.nBar];


    % plot the fitted curve
    if (fitted_curve.theta > 0)
        hold on
        plot(h.BinEdges(cut_idx:cutUp_idx),fitted_curve(h.BinEdges(cut_idx:cutUp_idx)),'LineWidth',3,'Color',[colors(c,:)])
        % used color, increment color counter

        % calculate chi squared and degrees of freedom
        ch2 = sum(((h.Values(cut_idx:cutUp_idx)-fitted_curve(h.BinEdges(cut_idx:cutUp_idx))').^2)./fitted_curve(h.BinEdges(cut_idx:cutUp_idx))');
        dof = size(h.Values(cut_idx:cutUp_idx),2);
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

% print out file with params
filename = [run.id, 'SPEoptParams.csv'];
filedir = ['C:\Users\GDD\Documents\MATLAB\Picosec\Michaela\SPEAnalysisPolyaMeanFinalVersion\', filename];
%filedir = ['C:\Users\Michaela\Documents\4th Year\CERN\MATLAB\SinglePEAnalysisPolyaMean\',filename];
%store_folder = ['\\eosproject-smb\eos\project\p\picosec\testbeam\2022_July_h4\Results\',filename];
%mkdir(store_folder);
writematrix(P,filename);
fileattrib(filedir,'+w','a');
% annotate and label plot
%text(min(h.BinEdges)+0.02*max(h.BinEdges),max(h.Values)+2,'\chi^2/dof = '+string(round(ch2)) + '/' + string(size(h.Values(cut_idx:cutUp_idx),2)) + ' = ' + string(nch2),'FontSize',14)
%text();


%% calculate the final mean from all means
format shortE
% mean
% remove 0s which correspond to fits that threw exceptions, do not want 0s
% biasing the mean
filter = mean_arr ~= 0;
mean_arr = mean_arr(filter);
chi2_arr = chi2_arr(filter);
mean_err_arr_temp = mean_err_arr(filter,:);
mean_err_arr = [];
mean_err_arr = mean_err_arr_temp;
min_arr = min_arr(filter);
mean = sum(mean_arr)/length(mean_arr);

%error propagation
delta_arr(1) = (1/2)*(sum((mean_err_arr(:,1)'-mean_arr).^2))^(1/2); % lower
delta_arr(2) = (1/2)*(sum((mean_err_arr(:,2)'-mean_arr).^2))^(1/2); % upper

% print values
str = sprintf('Mean: %f \n Upper bound: %f \n Lower bound: %f',mean,mean+delta_arr(2),mean-delta_arr(1));
display(str);

% str1 = strcat('Mean =  ',{' '},num2str(mean),' +/- ',num2str(delta_arr(2)));
% str2 = strcat('Upper Bound = ',num2str(mean+delta_arr(2)));
% str3 = strcat('Lower Bound = ',num2str(mean-delta_arr(1)));

str1 = ['Mean = ' num2str(mean) ' +/- ' num2str(delta_arr(2))];
str1 = sprintf('Mean = %0.5f V +/- %0.5f V',mean,delta_arr(2));
%str2 = ['Upper Bound = ' num2str(mean+delta_arr(2))];
%str3 = ['Lower Bound = ' num2str(mean-delta_arr(1))];
%str = 'Test'
%str = {str1,str2,str3};
%str = {strA};
annotation('textbox',[.50 .6 .35 .3], 'String',str1,'FitBoxToText','on','FontSize',18);
saveas(figure(1),[run.id, 'E-PeakAmplitude.png']);

%% plot the chi2, mean and error in the mean

hold off;
figure(2)
x_fig=10;
y_fig=10;
width=1200;
height=600;
set(gcf,'position',[x_fig,y_fig,width,height]);
%plot(min_arr,chi2_arr, 'Marker','o','Color',[colors(:)]);
scatter(min_arr,chi2_arr,35,[colors(1:c-1,:)],'filled');
title('\chi^2/DOF vs Polya Mean');
xlabel('nBar [V]');
ylabel('\chi^2/DOF');
ax = gca;
ax.FontSize = 20;
grid on;
saveas(figure(2),[run.id, 'MeanVsChi2.png']);

figure(4)
x_fig=10;
y_fig=10;
width=1200;
height=600;
set(gcf,'position',[x_fig,y_fig,width,height]);
hold on;

for i = 1:length(min_arr)
    e = errorbar(min_arr(i),mean_arr(i),abs(mean_arr(i)-mean_err_arr(i,1)'),abs(mean_arr(i)-mean_err_arr(i,2)'),'Marker','o','Color',[colors(i,:)]);
    e.LineWidth = 1
end

%errorbar(min_arr,mean_arr,abs(mean_arr-mean_err_arr(:,1)'),abs(mean_arr-mean_err_arr(:,2)'),'ob');
title('Polya Mean vs Lower Cutoff Voltage');
ylabel('nBar [V]');
xlabel('Lower Polya Cutoff [V]');
saveas(figure(4),[run.id, 'MeanVsCutoffError.png']);
ax = gca;
ax.FontSize = 20;
grid on;
hold off;

%% save some parameters to a csv file
A = ["Mean","Mean Error Lower","Mean Error Upper","Polya Min Lower Cutoff","Polya Upper Cutoff","nbins","N init","theta init","nBar init","exceptions"];
vars = [mean,delta_arr(1),delta_arr(2),min_v,max_v,n_bins,x0(1),x0(2),x0(3),exceptionCounter];
A = [A;vars];
filename = [run.id, 'info.csv'];
filedir = ['C:\Users\GDD\Documents\MATLAB\Picosec\Michaela\SPEAnalysisPolyaMeanFinalVersion\', filename];
%filedir = ['C:\Users\Michaela\Documents\4th Year\CERN\MATLAB\SinglePEAnalysisPolyaMean\',filename];
writematrix(A,filename);
fileattrib(filedir,'+w','a');

mean
strcat('Upper Bound = ',num2str(mean+delta_arr(2)))
strcat('Lower Bound = ',num2str(mean-delta_arr(1)))


%% save data to mat file and store all m files used for processing

%save([store_folder '\Run' run.id '_' run.oscilloscope '_plots.mat'], 'run', 'pad', 'line', 'area', 'small', 'g2', 'twalk'); % save interesting data for plots

%% store tabular data
txt_tab.names = {'No','EventID','REF_MCP_timestamp','DUT_MCP_timestamp','REF_MCP_maxY', 'DUT_MCP_maxY','REF_MCP_e_peak', 'DUT_MCP_e_peak', 'TrackerX', 'TrackerY'};
fid = fopen([store_folder '/Run' run.id '_tabular.txt'],'w');
% print header with names
for i=1:length(txt_tab.names)
    fprintf(fid,'%s\t',txt_tab.names{i});
end
fprintf(fid,'\n');
for k=1:length(MCP_data)
    fprintf(fid,'%d\t%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n', k, MM_data(k).event_id, MCP_data(k).sigmoid.timepoint, MM_data(k).sigmoid.timepoint, MCP_data(k).sig.max.y,...
        MM_data(k).sig.max.y, MCP_data(k).sig.charge.e_peak, MM_data(k).sig.charge.e_peak, trackerX(k),trackerY(k));
end
fclose(fid);


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
    save([store_folder '\Run' run.id '-' run.oscilloscope '.mat'], 'run', 'MM_data', 'MCP_data', 'time_diff', 'time_diff_sigmoid', 'MCP_maxy', 'MM_maxy', 'trackerX', 'trackerY');
    MM_data_samplingCut = MM_data(cut_sampling);
    MCP_data_samplingCut = MCP_data(cut_sampling);
    time_diff_samplingCut = time_diff(cut_sampling);
    time_diff_sigmoid_samplingCut = time_diff_sigmoid(cut_sampling);
    MCP_maxy_samplingCut = MCP_maxy(cut_sampling);
    MM_maxy_samplingCut = MM_maxy(cut_sampling);
    trackerX_samplingCut = trackerX(cut_sampling);
    trackerY_samplingCut = trackerY(cut_sampling);
    save([store_folder '\Run' run.id '-' run.oscilloscope 'SamplingCut.mat'], 'run', 'MM_data_samplingCut', 'MCP_data_samplingCut', 'time_diff_samplingCut', 'time_diff_sigmoid_samplingCut', 'MCP_maxy_samplingCut', 'MM_maxy_samplingCut', 'trackerX_samplingCut', 'trackerY_samplingCut');
end