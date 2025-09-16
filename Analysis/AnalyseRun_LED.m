%clear all
close all

%load input file
%in_file = 'C:\Users\GDD\Documents\Picosec\May22\Analysed\Run0123-GDD.mat';

shouldSaveMat = false;

shouldUseCoarseMedianEstimate = true; %prefilter data with max counts for median time diff estiamation

%%configuration
twalk.en = 1;                                   % enable timewalk correction


     min_v = 0.01;
     max_v = 0.07;

num_bins = 100;

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
runTitleString = [run.name ' ' run.oscilloscope ' DUT:C' opts_MM.chID '_' runInfoString ' '];

dataFilePath = append('/eos/project/p/picosec/testbeam/2024_September_h4/info/MatlabRunInfo');
runDataFileName=append('RunData',string(runID),'.mat');

%finish config
pad.curvature = [1,1]; %circle
if pad.isCircle == 0
    pad.curvature = [0,0]; %rect
end


%% define folder for storing data
store_folder = ['C:\Users\gdd.CERN\Desktop\a\measurements2' run.id '-' run.oscilloscope '-dut' opts_MM.chID '-' runInfoString];
mkdir(store_folder);



%% extract from structurearrays
k=1;
for i=1:length(MM_data)
    %for i=1:length(MCP_data)
    %time_MCP(i) = MCP_data(i).sigmoid.timepoint;
    time_MM(i) = MM_data(i).sigmoid.timepoint;
    blavg_MM(i) = MM_data(i).sig.blavg;   % baseline average
    %blavg_MCP(i) = MCP_data(i).sig.blavg; % baseline average MCP
    blrms_MM(i) = MM_data(i).sig.blrms;   % baseline RMS MM (noise)
    %blrms_MCP(i) = MCP_data(i).sig.blrms; % baseline RMS MCP (noise)
    %ymax_MCP(i)=MCP_data(i).sig.max.y;    % take maximum of MCP
    ymax_MM(i)=MM_data(i).sig.max.y;      % take maximum of MM
    ymin_MM(i)=MM_data(i).sig.min.y;      % take minimum of MM
    e_peak_MM(i)=MM_data(i).sig.charge.e_peak;
    %e_peak_MCP(i)=MCP_data(i).sig.charge.e_peak;% extract electrom peak charge
    riseTime(i)=MM_data(i).sigmoid.timepoint90-MM_data(i).sigmoid.timepoint10;% extract electrom peak charge

end




% do some initial resolution measurement


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

saveas(gcf,[store_folder '\Run' run.id '_RAW_signals.png'])



%% plot peak amplitude histogram cut included - Only MM

%cut_idx = 5;
%cutUp_idx = 90;

figure
h=histogram(MM_maxy,num_bins);




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

cut_idx = min_cut_idx;



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
[p, err, chi] = fminuit('polya_minuit',p0,fit_data(:,cut_idx:cutUp_idx),'-b','-c',cmd);
polya_cnt_fit = polya_minuit(p,fit_data(1,cut_idx:cutUp_idx));
plot(xbins(cut_idx:cutUp_idx),polya_cnt_fit,'Linewidth',2);
e_peak_amp.mean = sum(polya_cnt_fit.*xbins(cut_idx:cutUp_idx))/sum(polya_cnt_fit);
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

%%%


lengthMMData = length(MM_data);

numPoints = length(MM_data(1).waveformY);

Y_all = NaN(lengthMMData, numPoints);


figure
for i = 1:lengthMMData
    plot(MM_data(i).waveformX*10^9,MM_data(i).waveformY);
    Y_all(i,:) = MM_data(i).waveformY(:)';
%    hold on;
end


Y_avg1 = mean(Y_all,1,'omitnan');

    plot(MM_data(1).waveformX*10^9,Y_avg1);

xlim([0 100]);
xlabel('Time, ns')
ylabel('Amplitude, V')
title_str = sprintf('Mean of signals');
title(title_str)


saveas(gcf,[store_folder '\Run' run.id '_mean_signals.png'])



% % % % % % % % % % % % % % % % %%%

% % % % % % % % % % % % % % % 

%zip([store_folder '/Run' run.id '_' run.oscilloscope '_mfiles.zip'],{dir('*.m').name}); % save all m files used for procesing as a zip file
%zip([store_folder '/Run' run.id '_' run.oscilloscope '_input.zip'],in_file); % ssave input file as zip
if shouldSaveMat
    save([store_folder '\Run' run.id '-' run.oscilloscope '.mat'], 'run', 'MM_data', 'MCP_data', 'time_diff', 'time_diff_sigmoidTWCorr', 'MCP_maxy', 'MM_maxy', 'trackerX', 'trackerY');
    MM_data_samplingCut = MM_data(cut_sampling);
  %  MCP_data_samplingCut = MCP_data(cut_sampling);
    time_diff_samplingCut = time_diff(cut_sampling);
    time_diff_sigmoidTWCorr_samplingCut = time_diff_sigmoidTWCorr(cut_sampling);
  %  MCP_maxy_samplingCut = MCP_maxy(cut_sampling);
    MM_maxy_samplingCut = MM_maxy(cut_sampling);
   % trackerX_samplingCut = trackerX(cut_sampling);
   % trackerY_samplingCut = trackerY(cut_sampling);
    save([store_folder '\Run' run.id '-' run.oscilloscope 'SamplingCut.mat'], 'run', 'MM_data_samplingCut', 'MCP_data_samplingCut', 'time_diff_samplingCut', 'time_diff_sigmoidTWCorr_samplingCut', 'MCP_maxy_samplingCut', 'MM_maxy_samplingCut', 'trackerX_samplingCut', 'trackerY_samplingCut');
end