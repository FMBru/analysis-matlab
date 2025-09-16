%clear all
close all

%load input file
%in_file = 'C:\Users\GDD\Documents\Picosec\May22\Analysed\Run0123-GDD.mat';

shouldSaveMat = false;

shouldUseCoarseMedianEstimate = true; %prefilter data with max counts for median time diff estiamation

%%configuration
twalk.en = 1;                                   % enable timewalk correction

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
    pad.isCircle = 1;
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

dataFilePath = append('/eos/project/p/picosec/testbeam/2024_September_h4/info/MatlabRunInfo');
runDataFileName=append('RunData',string(runID),'.mat');

%finish config
pad.curvature = [1,1]; %circle
if pad.isCircle == 0
    pad.curvature = [0,0]; %rect
end


%% define folder for storing data
store_folder = ['D:\asena\measurements' run.id '-' run.oscilloscope '-dut' opts_MM.chID '-ref' opts_MCP.chID '-' runInfoString];
mkdir(store_folder);



%% extract from structurearrays
k=1;
for i=1:length(MCP_data)
    %for i=1:length(MCP_data)
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

end




%% do some initial resolution measurement


[counts,centers] = hist(time_diff_sigmoid,1000);
[maxCounts,maxIdx] = max(counts);
coarseMedianEstimate = centers(maxIdx);
time_minCoarse = coarseMedianEstimate - 3;     % predicted resolution 100ps cut 3 sigma
time_maxCoarse = coarseMedianEstimate + 3;     % left and 3 sigma right from median

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
time_min = time_avg_raw - 0.3;     % predicted resolution 100ps cut 3 sigma
time_max = time_avg_raw + 0.3;     % left and 3 sigma right from median


% find events within the cut and cut with respect to amplirudes and
% existance of the tracker data (GLOBAL CUT)
%with tracker
glbl_cut = time_diff_sigmoid>time_min & time_diff_sigmoid<time_max & ...
    MCP_maxy>0.001*max(MCP_maxy)&...
    MCP_maxy<0.999*max(MCP_maxy) 

time_diff_cut = time_diff_sigmoid(glbl_cut);
mean(time_diff_cut)
std(time_diff_cut)


% % make time walk correction
% if(twalk.en == 1)
%     time_diff_sigmoidTWCorr = time_diff_sigmoid - twalk_fn_minuit(twalk.p, e_peak_MM);
% else
%     time_diff_sigmoidTWCorr = time_diff_sigmoid;
% end


% plot single gauss histogram
h = histogram(time_diff_cut,100);
hold on
fit_data = [];
fit_data(1,:) = h.BinEdges(1:end-1) + h.BinWidth/2;
fit_data(2,:) = h.Values;
% fit_data(3,:) = yerr;


p0=[];
p0(1) = sum(fit_data(2,:))*h.BinWidth/std(h.Data)/sqrt(2*pi);   % total normalization factor
p0(2) = std(time_diff_cut);                       % sigma1
p0(3) = mean(time_diff_cut);     


step = [1 2 3];
l_b = [p0(1)*0.8 0.1*p0(2) p0(3)-0.1];
u_b = [p0(1)*1.2 10*p0(2) p0(3)+0.1];
StepBounds= [step; l_b; u_b]';
MinuitCommands = 'min; ret;';
chi_idx = fit_data(2,:) > 0;  % take into account only non-empty bins

[p, err, chi2_min] = fminuit('gauss1_minuit', p0,fit_data(:,chi_idx), '-c', MinuitCommands,'-s', StepBounds)

plot(fit_data(1,:),gauss1_minuit(p,fit_data(1,:)),'Linewidth',2);
message =  sprintf('  \\mu = %2.3f ns \\pm %2.3f ps\n',p(3),1000*err(3));
message = [message sprintf('  \\sigma_{fit} = %2.1f ps \\pm %2.3f ps\n',1000*p(2),1000*err(2))];
message = [message sprintf('  RMS_{hist} = %2.1f ps ',1000*std(time_diff_cut))];

xlabel('Time difference: PICOSEC vs Reference, ns');
ylabel('Events');
%xlim ([-0.2 0.2])
legend('RAW hist','Gauss fit');
grid
title_str = sprintf('Time resolution - Single Gauss fit');
title(title_str)
y_pos=get(gca,'ylim');
x_pos=get(gca,'xlim');
text(x_pos(1),0.88*y_pos(2),message)
movegui(gcf,'southeast');
saveas(gcf,[store_folder '\Run' run.id '_time_res-samplingArea-singleGauss.png'])


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
legend('REF','MM','MM fit');
title_str = sprintf('e-peak amplitude \\mu = %4.4f V U_{max} = %4.4f V', e_peak_amp.mean,xbins(e_peak_amp.max_idx));
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
title_str = sprintf('e-peak amplitude \\mu = %4.4f V U_{max} = %4.4f V', e_peak_amp.mean,xbins(e_peak_amp.max_idx));
%title_str = sprintf('e-peak amplitude \\mu = %4.4f V U_{max} = %4.4f V, (files %2.0f - %2.0f)', e_peak_amp.mean,xbins(e_peak_amp.max_idx),file_start,file_stop);
title(title_str)
saveas(gcf,[store_folder '\Run' run.id '_signalAmplitude_Hist_MM.png'])

file_start_str = sprintf(' %2.0f',file_start);
file_stop_str =sprintf(' %2.0f',file_stop);
saveas(gcf,[store_folder '\Run' run.id 'signalAmplitude_Hist_MM_from_' file_start_str '_to_' file_stop_str '.png'])


riseTimeCut = abs(riseTime)<10;


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
[~, e_peak_amp.max_idx] = max(polya_cnt_fit);

xlabel('Electron lead charge, pC');
ylabel('Events');
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


xlabel('Electron peak charge, pC');
ylabel('Events');
grid on
legend('Hist','Polya fit');
title_str = sprintf('%s \n e-peak charge \\mu = %4.4f pC Q_{max} = %4.4f pC',runTitleString, e_peak_amp.mean,xbins(e_peak_amp.max_idx));
title(title_str)

saveas(gcf,[store_folder '\Run' run.id '_echarge_hist.png'])


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




% %plot resolution vs. e-charge
figure
errorbar(twalk.e_peak,twalk.rms*1000,[],[],twalk.e_peak_err_n,twalk.e_peak_err_p,'o');
hold on
xlabel('Electron peak charge, pC')
ylabel('Resolution, ps')
title_str = sprintf('%s \n Resolution vs. Electron peak charge',runTitleString);
title(title_str)

grid
saveas(gcf,[store_folder '\Run' run.id '_res_vs_charge.png'])


%signals




%histo1

timeDiffCutLower = -0.7;
timeDiffCutHigher = -0.50;

timeDiffCutLowerStr = sprintf(' %4.2f',timeDiffCutLower);
timeDiffCutHigherStr = sprintf( '%4.2f', timeDiffCutHigher);

timeDiffMask = time_diff_sigmoid>timeDiffCutLower & time_diff_sigmoid<timeDiffCutHigher;

MM_dataSelection = MM_data(timeDiffMask);

lengthMMDataSelection = length(MM_dataSelection);

numPoints = length(MM_dataSelection(1).waveformY);

Y_all = NaN(lengthMMDataSelection, numPoints);


% figure
for i = 1:lengthMMDataSelection
    plot(MM_dataSelection(i).waveformX*10^9,MM_dataSelection(i).waveformY);
    Y_all(i,:) = MM_dataSelection(i).waveformY(:)';
    hold on;
end

xlabel('Time, ns')
ylabel('Amplitude, V')
title_str = sprintf('Signals for time difference: PICOSEC vs Reference from %4.2f ns to %4.2f ns',timeDiffCutLower,timeDiffCutHigher);
title(title_str)

saveas(gcf,[store_folder '\Run' run.id '_signals_from_' timeDiffCutLowerStr 'ns_to_' timeDiffCutHigherStr 'ns.png'])

figure
Y_avg1 = mean(Y_all,1,'omitnan');
plot(MM_dataSelection(1).waveformX*10^9,Y_avg1);

%xlim([-0.25 2]);
xlabel('Time, ns')
ylabel('Amplitude, V')
title_str = sprintf('Mean of signals for time difference: PICOSEC vs Reference from %4.2f ns to %4.2f ns',timeDiffCutLower,timeDiffCutHigher);
title(title_str)


saveas(gcf,[store_folder '\Run' run.id '_signals_from_' timeDiffCutLowerStr 'ns_to_' timeDiffCutHigherStr 'ns_mean1.png'])




%histo2

timeDiffCutLower = -0.47;
timeDiffCutHigher = -0.4;

timeDiffCutLowerStr = sprintf(' %4.2f',timeDiffCutLower);
timeDiffCutHigherStr = sprintf( '%4.2f', timeDiffCutHigher);

timeDiffMask = time_diff_sigmoid>timeDiffCutLower & time_diff_sigmoid<timeDiffCutHigher;

MM_dataSelection = MM_data(timeDiffMask);

lengthMMDataSelection = length(MM_dataSelection);

numPoints = length(MM_dataSelection(1).waveformY);

Y_all = NaN(lengthMMDataSelection, numPoints);



for i = 1:lengthMMDataSelection
    Y_all(i,:) = MM_dataSelection(i).waveformY(:)';
end


Y_avg2 = mean(Y_all,1,'omitnan');


%

timeDiffCutLower = -0.3;
timeDiffCutHigher = -0.1;

timeDiffCutLowerStr = sprintf(' %4.2f',timeDiffCutLower);
timeDiffCutHigherStr = sprintf( '%4.2f', timeDiffCutHigher);

timeDiffMask = time_diff_sigmoid>timeDiffCutLower & time_diff_sigmoid<timeDiffCutHigher;

MM_dataSelection = MM_data(timeDiffMask);

lengthMMDataSelection = length(MM_dataSelection);

numPoints = length(MM_dataSelection(1).waveformY);

Y_all = NaN(lengthMMDataSelection, numPoints);


for i = 1:lengthMMDataSelection
    Y_all(i,:) = MM_dataSelection(i).waveformY(:)';
end

Y_avg3 = mean(Y_all,1,'omitnan');


%

% timeDiffCutLower = -0.27;
% timeDiffCutHigher = -0.23;
% 
% timeDiffCutLowerStr = sprintf(' %4.2f',timeDiffCutLower);
% timeDiffCutHigherStr = sprintf( '%4.2f', timeDiffCutHigher);
% 
% timeDiffMask = time_diff_sigmoid>timeDiffCutLower & time_diff_sigmoid<timeDiffCutHigher;
% 
% MM_dataSelection = MM_data(timeDiffMask);
% 
% lengthMMDataSelection = length(MM_dataSelection);
% 
% numPoints = length(MM_dataSelection(1).waveformY);
% 
% Y_all = NaN(lengthMMDataSelection, numPoints);
% 
% 
% for i = 1:lengthMMDataSelection
%     Y_all(i,:) = MM_dataSelection(i).waveformY(:)';
% end
% 
% Y_avg4 = mean(Y_all,1,'omitnan');



figure
plot(MM_dataSelection(1).waveformX*10^9,Y_avg1);
hold on;
plot(MM_dataSelection(1).waveformX*10^9,Y_avg2);
hold on;
plot(MM_dataSelection(1).waveformX*10^9,Y_avg3);
% hold on;
% plot(MM_dataSelection(1).waveformX*10^9,Y_avg4);


xlabel('Time, ns')
ylabel('Amplitude, V')
title_str = sprintf('Mean PICOSEC signals');
title(title_str)

saveas(gcf,[store_folder '\Run' run.id '_signals_from_' timeDiffCutLowerStr 'ns_to_' timeDiffCutHigherStr 'ns_meanAll.png'])





%% plot peak amplitude histogram cut included - Only MM in these regions
figure
h=histogram(MM_maxy(timeDiffMask),25);
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
 title_str = sprintf('e-peak amplitude \\mu = %4.4f V U_{max} = %4.4f V', e_peak_amp.mean,xbins(e_peak_amp.max_idx));
 title(title_str)
saveas(gcf,[store_folder '\Run' run.id '_amplitude_from_' timeDiffCutLowerStr 'ns_to_' timeDiffCutHigherStr 'ns.png'])





%  %% For Marta 
% 
% SATlist = time_diff_sigmoidTWCorr(cut_sampling); %for further analysis in ROOT
% writematrix(SATlist(:), [store_folder '/Run' run.id '_' runInfoString '- SATlist.txt']);


if shouldSaveMat
    save([store_folder '\Run' run.id '-' run.oscilloscope '.mat'], 'run', 'MM_data', 'MCP_data', 'time_diff', 'time_diff_sigmoidTWCorr', 'MCP_maxy', 'MM_maxy', 'trackerX', 'trackerY');
    MM_data_samplingCut = MM_data(cut_sampling);
    MCP_data_samplingCut = MCP_data(cut_sampling);
    time_diff_samplingCut = time_diff(cut_sampling);
    time_diff_sigmoidTWCorr_samplingCut = time_diff_sigmoidTWCorr(cut_sampling);
    MCP_maxy_samplingCut = MCP_maxy(cut_sampling);
    MM_maxy_samplingCut = MM_maxy(cut_sampling);
    trackerX_samplingCut = trackerX(cut_sampling);
    trackerY_samplingCut = trackerY(cut_sampling);
    save([store_folder '\Run' run.id '-' run.oscilloscope 'SamplingCut.mat'], 'run', 'MM_data_samplingCut', 'MCP_data_samplingCut', 'time_diff_samplingCut', 'time_diff_sigmoidTWCorr_samplingCut', 'MCP_maxy_samplingCut', 'MM_maxy_samplingCut', 'trackerX_samplingCut', 'trackerY_samplingCut');
end