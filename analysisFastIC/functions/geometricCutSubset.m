clear all
close all

in_file = 'C:\Users\GDD\Documents\Picosec\May22\Analysed\Run103-GDD.mat';

alignToDUT = false; %if false, alignment to REF MCP

small.r = 5; %radius to select events

%pad size (visualisation)
pad.size = 11;     % MCP size in mm - Hamamatsu
%pad.size = 12;     % MCP size in mm - Planacon, 2x2 central pads
pad.cut_len = 8;   % move +/-8mm from median to find center
pad.isCircle = 1;   %plotting round if circle, otherwise square

%sampling area for 2D maps
%% calculate
area.step = 0.25;   % set grid resolution in mm
%area.size = 8;      % set half-size of the observed square area, mm - for
%pads and small MCP
area.size = 20;      % set half-size of the observed square area, mm - for large MCP
area.radius = 1;    % set radius of averaging circular window cut, mm


load(in_file);

runID = str2num(run.id);

%finish config
pad.curvature = [1,1]; %circle
if pad.isCircle == 0
    pad.curvature = [0,0]; %rect
end

%% define folder for storing data
store_folder = ['analyzed/Run' run.id '-' run.oscilloscope];
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
    %if(exitflag_MCP(i)+exitflag_MM(i)==0)
end

%% do some initial resolution measurement
% make simple cut with respect to median time and minimum amplitude magic
% numbers :-)
% MEDIAN/MEAN
time_avg_raw = median(time_diff_sigmoid);  % assume mu from the median vaule
%time_avg_raw = mean(time_diff_sigmoid);  % assume mu from the median vaule
time_min = time_avg_raw - 0.3;     % predicted resolution 100ps cut 3 sigma
time_max = time_avg_raw + 0.3;     % left and 3 sigma right from median

% find events within the cut and cut with respect to amplirudes and
% existance of the tracker data (GLOBAL CUT)

%with tracker
glbl_cut = time_diff_sigmoid>time_min & time_diff_sigmoid<time_max & MCP_maxy>0.01*max(MCP_maxy) & MM_maxy>0.01*max(MM_maxy)& MM_maxy<0.95*max(MM_maxy)& MCP_maxy<0.95*max(MCP_maxy) & trackerX~=0 & trackerY~=0;

%without tracker
%glbl_cut = time_diff_sigmoid>time_min & time_diff_sigmoid<time_max & MCP_maxy>0.01*max(MCP_maxy) & MM_maxy>0.01*max(MM_maxy)& MM_maxy<0.95*max(MM_maxy)& MCP_maxy<0.5*max(MCP_maxy);

% make cut on time difference
time_diff_cut = time_diff_sigmoid(glbl_cut);
mean(time_diff_cut)
std(time_diff_cut)


%% calculate DUT MCP centre if alignToDUT == false

if alignToDUT == true
    
    % calculate pad center
    % MEDIAN/MEAN
    pad.xc_med = median(trackerX(glbl_cut));
    pad.yc_med = median(trackerY(glbl_cut));
    %pad.xc_med = mean(trackerX(glbl_cut));
    %pad.yc_med = mean(trackerY(glbl_cut));
    pad.x_idx_cut = trackerX > (pad.xc_med-pad.cut_len) & trackerX < (pad.xc_med+pad.cut_len) & glbl_cut;
    pad.xc_n = mean(trackerX(pad.x_idx_cut)); % naive mean for x
    pad.y_idx_cut = trackerY > (pad.yc_med-pad.cut_len) & trackerY < (pad.yc_med+pad.cut_len) & glbl_cut;
    pad.yc_n = mean(trackerY(pad.y_idx_cut)); % naive mean for y
    
    % Plot x and y hit histograms
    figure;
    title(['Hits on DUT detector (X, Y projections)' run.name]);
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
    
    %saveas(gcf,[store_folder '\Run' run.id '_hits_DUT-Subset.png'])
    
    % find x center from electron peak mean charge
    for i = 1:length(pad.h_x.Values)
        tmp_cut = trackerX>pad.h_x.BinEdges(i) & trackerX<pad.h_x.BinEdges(i)+ pad.h_x.BinWidth;
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
    title(['E-peak mean over x-axis ' run.name]);
    grid on
    %saveas(gcf,[store_folder '/Run' run.id '_DUT_epeak_X_projSubset.png'])
    
    
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
    title(['E-peak mean over y-axis ' run.name]);
    grid on
    %saveas(gcf,[store_folder '\Run' run.id '_DUT_epeak_Y_projSubset.png'])
    
    %pad.xc = pad.xc_n;
    %pad.yc = pad.yc_n;
end

%% calculate REF MCP centre

if alignToDUT == false
    
    % calculate pad center
    % MEDIAN/MEAN
    pad.xc_med = median(trackerX(glbl_cut));
    pad.yc_med = median(trackerY(glbl_cut));
    %pad.xc_med = mean(trackerX(glbl_cut));
    %pad.yc_med = mean(trackerY(glbl_cut));
    pad.x_idx_cut = trackerX > (pad.xc_med-pad.cut_len) & trackerX < (pad.xc_med+pad.cut_len) & glbl_cut;
    pad.xc_n = mean(trackerX(pad.x_idx_cut)); % naive mean for x
    pad.y_idx_cut = trackerY > (pad.yc_med-pad.cut_len) & trackerY < (pad.yc_med+pad.cut_len) & glbl_cut;
    pad.yc_n = mean(trackerY(pad.y_idx_cut)); % naive mean for y
    
    % Plot x and y hit histograms
    figure;
    title(['Hits on REF detector (X, Y projections)' run.name]);
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
    
    %saveas(gcf,[store_folder '\Run' run.id '_hits_REF-Subset.png'])
    
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
    pad.xc_err = err(1);
    % plot to see how parabolic fit looks like
    figure
    bar(fit_data(1,:),fit_data(2,:));
    hold on
    plot(fit_data(1,:),parabola4_minuit(p, fit_data(1,:)),'LineWidth',2);
    xlabel('x-axis, mm');
    ylabel('Charge, pC');
    legend('RAW', 'Fit')
    title(['E-peak REF MCP mean over x-axis ' run.name]);
    grid on
    %saveas(gcf,[store_folder '/Run' run.id '_REF_epeak_X_proj-Subset.png'])
    
    
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
    pad.yc_err = err(1);
    % plot to see how parabolic fit looks like
    figure
    bar(fit_data(1,:),fit_data(2,:));
    hold on
    plot(fit_data(1,:),parabola4_minuit(p, fit_data(1,:)),'LineWidth',2);
    xlabel('y-axis, mm');
    ylabel('Charge, pC');
    legend('RAW', 'Fit')
    title(['E-peak REF MCP mean over y-axis ' run.name]);
    grid on
   % saveas(gcf,[store_folder '\Run' run.id '_REF_epeak_Y_proj-Subset.png'])
    
    %pad.xc = pad.xc_n;
    %pad.yc = pad.yc_n;
    
end

%circular cut_small -> select events to preserve
cut_small = (((trackerX - pad.xc).^2 + (trackerY - pad.yc).^2) < small.r^2) & glbl_cut;

% plot hits
figure;
scatter(trackerX, trackerY,'.')   % plot all hits
hold on;
scatter(trackerX(cut_small), trackerY(cut_small),'.'); % plot hits with global cut
%plot(line.x, line.y, 'LineWidth', 2);              % plot observation line
scatter(pad.xc,pad.yc);                            % plot center point
rectangle('Position',[pad.xc-pad.size/2 pad.yc-pad.size/2 pad.size pad.size], 'Curvature',pad.curvature, 'LineWidth',1);    %plot active area
rectangle('Position',[pad.xc-(2*small.r)/2 pad.yc-(2*small.r)/2 (2*small.r) (2*small.r)], 'Curvature',[1,1], 'LineWidth',1,'EdgeColor','green');    %plot selected area
axis equal
xlabel('x-axis, mm');
ylabel('y-axis, mm');
grid on
str_title = sprintf('%s tracker data', run.name);
title(str_title);
legend('Trigger','Hits','Observation path','DUT MCP center','Circle for resolution measurement');
xlim([min(trackerX(glbl_cut)) max(trackerX(glbl_cut))]);
ylim([min(trackerY(glbl_cut)) max(trackerY(glbl_cut))]);
saveas(gcf,[store_folder '\Run' run.id '_hits_overlay-Subset.png'])



%ring-shaped cut small
%cut_small = (((trackerX - pad.xc).^2 + (trackerY - pad.yc).^2) > small.r^2) & (((trackerX - pad.xc).^2 + (trackerY - pad.yc).^2) < medium.r^2) & glbl_cut;

mean(time_diff_sigmoid(cut_small))
std(time_diff_sigmoid(cut_small))

MM_data = MM_data(cut_small);
MCP_data = MCP_data(cut_small);
time_diff = time_diff(cut_small);
time_diff_sigmoid = time_diff_sigmoid(cut_small);
MCP_maxy = MCP_maxy(cut_small);
MM_maxy = MM_maxy(cut_small);
trackerX = trackerX(cut_small);
trackerY = trackerY(cut_small);


save(['C:\Users\GDD\Documents\Picosec\May22\Analysed\Run' run.id '-' run.oscilloscope 'GeoSubset.mat'], 'run', 'MM_data', 'MCP_data', 'time_diff', 'time_diff_sigmoid', 'MCP_maxy', 'MM_maxy', 'trackerX', 'trackerY');