close all

shouldLoad = true;



%runs for geometric scan of photek
%Run73 - center, 4000V ok
%Run76 - bottom left, 2cm, 4000V ok
%Run78 - bottom left, 1cm, 4000V ok
%Run81 - top right, 1.5cm, 4000V ok
%Run82 - top left, 1.5cm, 4000V ok
%Run83 - bottom rigt, 1.5cm, 4000V ok


%open multiple files, implement cuts and

;

if shouldLoad
    clear all
 
    
    
filesList = ["Run073-GDDGeoSubset.mat" "Run076-GDDGeoSubset.mat" "Run078-GDDGeoSubset.mat" "Run081-GDDGeoSubset.mat" "Run082-GDDGeoSubset.mat" "Run083-GDDGeoSubset.mat"];

%planacon scan
filesList = ["\\eosproject-smb\eos\project\p\picosec\testbeam\2022_May_h4\Results\Run101-GDD\Run101-GDDSamplingCut.mat" 
    "\\eosproject-smb\eos\project\p\picosec\testbeam\2022_May_h4\Results\Run103-GDD\Run103-GDDSamplingCut.mat" 
    "\\eosproject-smb\eos\project\p\picosec\testbeam\2022_May_h4\Results\Run104-GDD\Run104-GDDSamplingCut.mat"  ];
%filesList = ["Run073-GDDGeoSubset.mat" "Run081-GDDGeoSubset.mat"];

    for (ff=1:length(filesList))
        %loop through files
        in_file = filesList(ff);
        fileData = load(in_file);
        
%         if ff==1
%             MM_data = fileData.MM_data;
%             MCP_data = fileData.MCP_data;
%             time_diff = fileData.time_diff;
%             time_diff_sigmoid = fileData.time_diff_sigmoid;
%             MCP_maxy = fileData.MCP_maxy;
%             MM_maxy = fileData.MM_maxy;
%             trackerX = fileData.trackerX;
%             trackerY = fileData.trackerY;
%         else
%             MM_data = [MM_data,fileData.MM_data];
%             MCP_data = [MCP_data,fileData.MCP_data];
%             time_diff = [time_diff,fileData.time_diff];
%             time_diff_sigmoid = [time_diff_sigmoid,fileData.time_diff_sigmoid];
%             MCP_maxy = [MCP_maxy,fileData.MCP_maxy];
%             MM_maxy = [MM_maxy,fileData.MM_maxy];
%             trackerX = [trackerX,fileData.trackerX];
%             trackerY = [trackerY,fileData.trackerY];
%             
%         end

        
                if ff==1
            MM_data = fileData.MM_data_samplingCut;
            MCP_data = fileData.MCP_data_samplingCut;
            time_diff = fileData.time_diff_samplingCut;
            time_diff_sigmoid = fileData.time_diff_sigmoid_samplingCut;
            MCP_maxy = fileData.MCP_maxy_samplingCut;
            MM_maxy = fileData.MM_maxy_samplingCut;
            trackerX = fileData.trackerX_samplingCut;
            trackerY = fileData.trackerY_samplingCut;
        else
            MM_data = [MM_data,fileData.MM_data_samplingCut];
            MCP_data = [MCP_data,fileData.MCP_data_samplingCut];
            time_diff = [time_diff,fileData.time_diff_samplingCut];
            time_diff_sigmoid = [time_diff_sigmoid,fileData.time_diff_sigmoid_samplingCut];
            MCP_maxy = [MCP_maxy,fileData.MCP_maxy_samplingCut];
            MM_maxy = [MM_maxy,fileData.MM_maxy_samplingCut];
            trackerX = [trackerX,fileData.trackerX_samplingCut];
            trackerY = [trackerY,fileData.trackerY_samplingCut];
            
        end

        %clear MM_data MCP_data time_diff time_diff_sigmoid MCP_maxy MMM_maxy trackerX trackerY;
    end
    
end

%hard coded center of device (determned from amplitude profile)
centerX=12.7912; %from mean ampl of Run86
centerY=30.6735

centerX=16; %from circle on mean by eye
centerY=28

run.id = 'PlanaconScan';
run.oscilloscope = 'GDD';
run.year = '2022/05 ';
run.name = ['BEAM ' run.year ' RUN ' run.id];

%%configuration

% center resolution circular - radius
small.r = 20; %used for only circle
medium.r = 2.5; %used for ring only

alignToDUT = false; %if false, alignment to REF MCP

%pad size (visualisation)
pad.size = 11;     % MCP size in mm - Hamamatsu
%pad.size = 12;     % MCP size in mm - Planacon, 2x2 central pads
pad.size = 40;     % MCP size in mm - Planacon, 2x2 central pads
pad.cut_len = 8;   % move +/-8mm from median to find center
pad.isCircle = 1;   %plotting round if circle, otherwise square

%sampling area for 2D maps
%% calculate
area.step = 0.25;   % set grid resolution in mm
%area.size = 8;      % set half-size of the observed square area, mm - for
%pads and small MCP
area.size = 30;      % set half-size of the observed square area, mm - for large MCP
area.radius = 1;    % set radius of averaging circular window cut, mm


%load(dataFileLocalPath);


%finish config
pad.curvature = [1,1]; %circle
if pad.isCircle == 0
    pad.curvature = [0,0]; %rect
end




%% define folder for storing data
store_folder = ['\\eosproject-smb\eos\project\p\picosec\testbeam\2022_May_h4\Results\Run' run.id '-' run.oscilloscope '-Merged-'];
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
time_avg_raw = median(time_diff_sigmoid);  % assume mu from the median vaule
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

% Time walk analysis (this needs to be improoved)
twalk.en = 0;                                   % enable timewalk correction
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
grid
saveas(gcf,[store_folder '\Run' run.id '_timewalk-Merged.png'])

% plot resolution vs. e-charge
figure
errorbar(twalk.e_peak,twalk.rms*1000,[],[],twalk.e_peak_err_n,twalk.e_peak_err_p,'o');
hold on
xlabel('Electron peak charge, pC')
ylabel('Resolution, ps')
grid
saveas(gcf,[store_folder '\Run' run.id '_res_vs_charge-Merged.png'])

% make time walk correction
if(twalk.en == 1)
    time_diff_sigmoid = time_diff_sigmoid - twalk_fn_minuit(twalk.p, e_peak_MM);
end


%% calculate DUT MCP centre if alignToDUT == false

if alignToDUT == true
    
    % calculate pad center
    pad.xc_med = median(trackerX(glbl_cut));
    pad.yc_med = median(trackerY(glbl_cut));
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
    %xlim([0 50]);
    %ylim([0 50]);
    
    saveas(gcf,[store_folder '\Run' run.id '_hits_DUT-Merged.png'])
    
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
    saveas(gcf,[store_folder '/Run' run.id '_DUT_epeak_X_proj-Merged.png'])
    
    
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
    saveas(gcf,[store_folder '\Run' run.id '_DUT_epeak_Y_proj-Merged.png'])
    
    %pad.xc = pad.xc_n;
    %pad.yc = pad.yc_n;
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
    %xlim([0 50]);
    %ylim([0 50]);
    
    saveas(gcf,[store_folder '\Run' run.id '_hits_REF-Merged.png'])
    
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
    saveas(gcf,[store_folder '/Run' run.id '_REF_epeak_X_proj-Merged.png'])
    
    
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
    saveas(gcf,[store_folder '\Run' run.id '_REF_epeak_Y_proj-Merged.png'])
    
    %pad.xc = pad.xc_n;
    %pad.yc = pad.yc_n;
    
end

pad.xc = centerX;
pad.yc = centerY;

%% generate observation line and shift to MCP center position
line.angle = 45;     % line angle
line.step = 0.5;     % spacing in mm
line.half_len = 5.5;   % half length of the line
line.pts = -line.half_len:line.step:line.half_len;
line.x = cos(line.angle*pi/180)*line.pts + pad.xc;
line.y = sin(line.angle*pi/180)*line.pts + pad.yc;
line.n = length(line.pts);

% plot hits
figure;
scatter(trackerX, trackerY,'.')   % plot all hits
hold on;
scatter(trackerX(glbl_cut), trackerY(glbl_cut),'.'); % plot hits with global cut
plot(line.x, line.y, 'LineWidth', 2);              % plot observation line
scatter(pad.xc,pad.yc);                            % plot center point
%rectangle('Position',[pad.xc-pad.size/2 pad.yc-pad.size/2 pad.size pad.size], 'Curvature',pad.curvature, 'LineWidth',1);    %plot active area
%rectangle('Position',[pad.xc-(2*small.r)/2 pad.yc-(2*small.r)/2 (2*small.r) (2*small.r)], 'Curvature',[1,1], 'LineWidth',1,'EdgeColor','green');    %plot selected area
axis equal
xlabel('x-axis, mm');
ylabel('y-axis, mm');
grid on
str_title = sprintf('%s tracker data', run.name);
title(str_title);
legend('Trigger','Hits','Observation path','DUT MCP center','Circle for resolution measurement');
%xlim([min(trackerX(glbl_cut)) max(trackerX(glbl_cut))]);
%ylim([min(trackerY(glbl_cut)) max(trackerY(glbl_cut))]);
xlim([0 50]);
ylim([0 50]);
saveas(gcf,[store_folder '\Run' run.id '_hits_overlay-Merged.png'])



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
str_title = sprintf('%s:\n SAT line at %2.1f°',run.name, line.angle);
title(str_title);
xlabel('Distance from pad centre, mm');
ylabel('SAT, ps');
grid on;
saveas(gcf,[store_folder '\Run' run.id '_SAT_over_' num2str(line.angle) 'deg_line-Merged.png'])


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
set(h, 'EdgeColor', 'none');
%set(gca, 'XDir','reverse')
axis equal
colorbar
xlabel('x-axis, mm');
ylabel('y-axis, mm');
zlim([0 10]);

h = colorbar;
h.Label.String = 'SAT, ns';
h.Label.Position(1) = 3;
str_title = sprintf('%s:\n SAT over the PAD (\\phi_{avg} = %2.1f mm)', run.name, 2*area.radius);
title(str_title);
saveas(gcf,[store_folder '\Run' run.id '_SAT_map-Merged.png'])

%% plot time resolution over area
figure
h=pcolor(area.xx,area.yy,area.rms*1000);
hold on
%rectangle('Position',[pad.xc-pad.size/2 pad.yc-pad.size/2 pad.size pad.size], 'Curvature',pad.curvature, 'LineWidth',1);

%overlay nominal size
rectangle('Position',[centerX-pad.size/2 centerY-pad.size/2 pad.size pad.size], 'Curvature',pad.curvature, 'LineWidth',1);
set(h, 'EdgeColor', 'none');
axis equal
xlabel('x-axis, mm');
ylabel('y-axis, mm');
caxis([0 10]);

h = colorbar;
h.Label.String = 'Time resolution, ps';
h.Label.Position(1) = 3;
str_title = sprintf('%s:\n Time resolution over the PAD (\\phi_{avg} = %2.1f mm)', run.name, 2*area.radius);
title(str_title);
saveas(gcf,[store_folder '\Run' run.id '_RES_map-Merged.png'])

%% plot SAT mean errors
figure
h=pcolor(area.xx,area.yy,area.sat_err*1000);
hold on
rectangle('Position',[pad.xc-pad.size/2 pad.yc-pad.size/2 pad.size pad.size],  'Curvature',pad.curvature,'LineWidth',1);
set(h, 'EdgeColor', 'none');
axis equal
colorbar
xlabel('x-axis, mm');
ylabel('y-axis, mm');
h = colorbar;
h.Label.String = 'Error, ps';
h.Label.Position(1) = 3;
str_title = sprintf('%s:\n Mean error over the PAD (\\phi_{avg} = %2.1f mm)', run.name, 2*area.radius);
title(str_title);
saveas(gcf,[store_folder '\Run' run.id '_SAT_err-Merged.png'])


%% plot DUT MCP amplitude over area
figure
h=pcolor(area.xx,area.yy,area.amp);
hold on
rectangle('Position',[pad.xc-pad.size/2 pad.yc-pad.size/2 pad.size pad.size],  'Curvature',pad.curvature,'LineWidth',2);
set(h, 'EdgeColor', 'none');
axis equal
colorbar
xlabel('x-axis, mm');
ylabel('y-axis, mm');
h = colorbar;
h.Label.String = 'Amplitude, V';
h.Label.Position(1) = 3;
str_title = sprintf('%s:\n DUT MCP amplitude over the PAD (\\phi_{avg} = %2.1f mm)', run.name, 2*area.radius);
title(str_title);
saveas(gcf,[store_folder '\Run' run.id '_DUT_MCP_amp-Merged.png'])

%% plot DUT MCP amplitude over area as the 3D plot
figure
h=surf(area.xx,area.yy,area.amp);
hold on
rectangle('Position',[pad.xc-pad.size/2 pad.yc-pad.size/2 pad.size pad.size],  'Curvature',pad.curvature,'LineWidth',2);
set(h, 'EdgeColor', 'none');
%axis equal
colorbar
xlabel('x-axis, mm');
ylabel('y-axis, mm');
h = colorbar;
h.Label.String = 'Amplitude, V';
h.Label.Position(1) = 3;
str_title = sprintf('%s:\n DUT MCP amplitude over the PAD (\\phi_{avg} = %2.1f mm)', run.name, 2*area.radius);
title(str_title);
saveas(gcf,[store_folder '\Run' run.id '_DUT_MCP_amp_3D-Merged.png'])

%% plot E-peak charge over area
figure
h=pcolor(area.xx,area.yy,area.e_peak);
hold on
rectangle('Position',[pad.xc-pad.size/2 pad.yc-pad.size/2 pad.size pad.size],  'Curvature',pad.curvature, 'LineWidth',2);
set(h, 'EdgeColor', 'none');
axis equal
colorbar
xlabel('x-axis, mm');
ylabel('y-axis, mm');
h = colorbar;
h.Label.String = 'Charge, pC';
h.Label.Position(1) = 3;
str_title = sprintf('%s:\n DUT MCP E-peak charge (\\phi_{avg} = %2.1f mm)', run.name, 2*area.radius);
title(str_title);
saveas(gcf,[store_folder '\Run' run.id '_charge_map-Merged.png'])

%% plot REF MCP Amplitude over area
figure
h = pcolor(area.xx,area.yy,area.amp_MCP);
hold on
rectangle('Position',[pad.xc-pad.size/2 pad.yc-pad.size/2 pad.size pad.size], 'Curvature',pad.curvature, 'LineWidth',2);
set(h, 'EdgeColor', 'none');
axis equal
colorbar
xlabel('x-axis, mm');
ylabel('y-axis, mm');
h = colorbar;
h.Label.String = 'Amplitude, V';
h.Label.Position(1) = 3;
str_title = sprintf('%s:\n REF MCP amplitude over the PAD (\\phi_{avg} = %2.1f mm)', run.name, 2*area.radius);
title(str_title);
saveas(gcf,[store_folder '\Run' run.id '_REF_MCP_map-Merged.png'])

%% plot REF MCP Amplitude over area 3D
figure
h = surf(area.xx,area.yy,area.amp_MCP);
hold on
rectangle('Position',[pad.xc-pad.size/2 pad.yc-pad.size/2 pad.size pad.size], 'Curvature',pad.curvature,'LineWidth',2);
set(h, 'EdgeColor', 'none');
%axis equal
colorbar
xlabel('x-axis, mm');
ylabel('y-axis, mm');
h = colorbar;
h.Label.String = 'Amplitude, V';
h.Label.Position(1) = 3;
str_title = sprintf('%s:\n REF MCP amplitude over the PAD (\\phi_{avg} = %2.1f mm)', run.name, 2*area.radius);
title(str_title);
saveas(gcf,[store_folder '\Run' run.id '_REF_MCP_map_3D-Merged.png'])


%circular cut_small
%cut_small = (((trackerX - pad.xc).^2 + (trackerY - pad.yc).^2) < small.r^2) & glbl_cut;

%sample time at fixed position
cut_small = (((trackerX - centerX).^2 + (trackerY - centerY).^2) < small.r^2) & glbl_cut;

%array definining radii for sampling rings
ringCutArray = [[0 5];[5 10];[10 15];[15 20];[20 25];[25 30]];
ringRMSArray = [];
meanRadiusArray = [];

%loop through ring cuts
for (pos=1:length(ringCutArray))

    innerRadius = ringCutArray(pos,1);
    outerRadius = ringCutArray(pos,2);
    meanRadius = (innerRadius+outerRadius)/2;
    meanRadiusArray = [meanRadiusArray;meanRadius];
    
    radiusInfo = append("R",int2str(innerRadius),"-",int2str(outerRadius));
    
%ring-shaped cut small
cut_small = (((trackerX - pad.xc).^2 + (trackerY - pad.yc).^2) > innerRadius^2) & (((trackerX - pad.xc).^2 + (trackerY - pad.yc).^2) < outerRadius^2) & glbl_cut;

mean(time_diff_sigmoid(cut_small))
std(time_diff_sigmoid(cut_small))

%fractionEventsSelectedInPad = length(time_diff_sigmoid(cut_small))/length(time_diff_sigmoid);
fractionEventsSelectedForTiming = length(time_diff_sigmoid(cut_small))/length(time_diff_sigmoid);

% plot dual gauss histogram
figure
%rng(abs(pd.mu));
h = histogram(time_diff_sigmoid(cut_small),100);
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
    p0(2) = std(time_diff_sigmoid(cut_small));                       % sigma1
    p0(3) = mean(time_diff_sigmoid(cut_small));                      % mean1
    p0(4) = comb(i);                                  % combination factotr
    p0(5) = std(time_diff_sigmoid(cut_small));                       % sigma2
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
message = [message sprintf('  RMS_{tot} = %2.1f ps ',1000*std(time_diff_sigmoid(cut_small)))];

xlabel('Time difference, ns');
ylabel('events');
legend('RAW hist','Gauss combined','Gauss 1', 'Gauss 2');
grid
title_str = sprintf('2Gauss %s (\\phi = %2.1f mm)', run.name, 2*small.r );
title(title_str)
y_pos=get(gca,'ylim');
x_pos=get(gca,'xlim');
zlim([0 10]);
text(x_pos(1),0.75*y_pos(2),message)
 saveas(gcf,[store_folder '\Run' run.id '_center_res-Merged-' convertStringsToChars(radiusInfo) '.png'])
 
 rmsResolution = 1000*std(time_diff_sigmoid(cut_small));
 
    ringRMSArray = [ringRMSArray;rmsResolution];
end

%% plot of RMS vs rings
figure
plot(meanRadiusArray,ringRMSArray,'o');
str_title = sprintf('%s:\n Time resolution as function of distance from center (rings)', run.name);
xlabel('Distance from center (ring center) (mm)');
ylabel('Time resolution RMS (ps)');
title(str_title);
xlim([0 30]);
ylim([0 50]);
grid on
saveas(gcf,[store_folder '\Run' run.id '_TimeRes-Rings-Merged.png'])

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
str_title = sprintf('%s:\n DUT MCP amplitude cross-sections ',run.name);
title(str_title);
legend(line_amp.legend, 'Location', 'Best');
xlabel('Distance from pad centre, mm');
ylabel('Amplitude, V');
grid on;
saveas(gcf,[store_folder '\Run' run.id '_DUT_amp_cross-Merged.png'])


%% save data to mat file and store all m files used for processing

save([store_folder '/Run' run.id '_' run.oscilloscope '_plots.mat'], 'run', 'pad', 'line', 'area', 'small', 'g2', 'twalk'); % save interesting data for plots

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
uploadCommand = append('pscp -r -pw ','Win_Admin',' C:\Users\GDD\Documents\MATLAB\Picosec\Analysis\analyzed\Run',run.id,'-',run.oscilloscope,'Merged ','gdd','@lxplus.cern.ch:/eos/project/p/picosec/testbeam/2022_May_h4/Results');
uploadResult = system(convertStringsToChars(uploadCommand))


%zip([store_folder '/Run' run.id '_' run.oscilloscope '_mfiles.zip'],{dir('*.m').name}); % save all m files used for procesing as a zip file
%zip([store_folder '/Run' run.id '_' run.oscilloscope '_input.zip'],in_file); % ssave input file as zip



