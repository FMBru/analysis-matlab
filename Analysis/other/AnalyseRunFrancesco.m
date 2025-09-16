%VERSION DATE: 30/06/2025 - FRANCESCO (summer student)
%code used to analyze single photoelectron response of PICOSEC using UV LED
% there's no timing info, used only to check if it possible to extract the
% number of p.e. emitted knowing the single p.e. response when using the same DUT 

close all

shouldSaveMat = false;
use_fixed_binning = true;

%% define folder for storing data
% add_str = '';
% if opts_MM.en_noiseRejection
%     add_str = [add_str '_rej'];
% end
% if opts_MM.en_filter
%     add_str = [add_str '_filt'];
% end

store_folder = ['\\eosproject-smb\eos\project\p\picosec\lab\Marta\PhotoDetectorResTi\Results\PhotoDetector_A275V_C' run.id 'V\AnalysisPlots\' runInfoString '\'];% add_str];
if exist(store_folder, 'dir') == 7
    rmdir(store_folder, 's');
end

mkdir(store_folder);

% mat_folder = ['\\eosproject-smb\eos\project\p\picosec\lab\Francesco\UV_LED_labtest\metallic_singlechannel\MatFile\Ti_C100_trigger\highStat\50mV-20ns\filt'];% add_str];
% % if exist(mat_folder, 'dir') == 7
% %     rmdir(mat_folder, 's');
% % end
% 
% if exist(mat_folder, 'dir') == 0
%     mkdir(mat_folder);
% end


%% extract from structurearrays
k=1;
for i=1:length(MM_data)
%     time_MM(i) = MM_data(i).sigmoid.timepoint;
%     blavg_MM(i) = MM_data(i).sig.blavg;   % baseline average
%     blrms_MM(i) = MM_data(i).sig.blrms;   % baseline RMS MM (noise)
%     ymax_MM(i)=MM_data(i).sig.max.y;      % take maximum of MM
%     ymin_MM(i)=MM_data(i).sig.min.y;      % take minimum of MM
    e_peak_MM(i)=MM_data(i).sig.charge.e_peak; %electron peak charge
    e_peak_width(i)=MM_data(i).sig.e_peak_width; %electron peak width = epeak endpoint - epeak startpoint
    riseTime(i)=MM_data(i).sigmoid.timepoint90-MM_data(i).sigmoid.timepoint10;% extract rise time
    all_charge(i) = MM_data(i).sig.charge.all;
    
    if opts_MM.en_filter
        %ion_tail(i) = MM_data(i).sig.ionTailAvg;    %ion tail average value
        isMultipeak(i) = MM_data(i).sig.is_multipeak; %logical array to know if it has more than one peak
    end
end

eventsValidDUT = logical(eventsValidDUT);
bin_width_amp = (max(MM_maxy) - min(MM_maxy))/100;
bin_width_epeak = (max(e_peak_MM) - min(e_peak_MM))/100;
bin_width_all = (max(all_charge) - min(all_charge))/100;

%% plot peak amplitude histogram cut included - Only MM
figure
if (opts_MM.en_filter)
    isMultipeak = logical(isMultipeak);
%     h=histogram(MM_maxy(~isMultipeak),100); 
end

if use_fixed_binning
    histogram(MM_maxy,'BinWidth', bin_width_amp, 'FaceColor','r'); 
    hold on
    h=histogram(MM_maxy(eventsValidDUT),'BinWidth', bin_width_amp, 'FaceColor','b'); 
else
    h=histogram(MM_maxy(eventsValidDUT),100, 'FaceColor','b'); 
    hold on
end

initial_bin = 5;


xbins = h.BinEdges(1:end-1)+h.BinWidth/2;
fit_data = [];
fit_data(1,:)=initial_bin:length(xbins);
fit_data(2,:) = h.Values(fit_data(1, :));
% fit_data(3,:) = yerr;
p0=[];
p0(1) = sum(h.Values)*h.BinWidth;   % normalization factor
p0(2) = 1;
p0(3) = 0.1;
cmd='min; ret';
[p, err, chi] = fminuit('polya_minuit',p0,fit_data(:,1:end),'-b','-c',cmd);
polya_cnt_fit = polya_minuit(p,fit_data(1,:));
plot(xbins(fit_data(1,:)),polya_cnt_fit,'Linewidth',2);
e_peak_amp.mean = sum(polya_cnt_fit.*xbins(fit_data(1,:)))/sum(polya_cnt_fit);
[dummy, e_peak_amp.max_idx] = max(polya_cnt_fit);

xlabel('Signal amplitude, V')
ylabel('Events');
ylim([0 1500]);
grid on
legend('wNoise', 'MM','MM fit');
title_str = sprintf('e-peak amplitude \\mu = %4.4f V U_{max} = %4.4f V', e_peak_amp.mean,xbins(e_peak_amp.max_idx));
%title_str = sprintf('e-peak amplitude \\mu = %4.4f V U_{max} = %4.4f V, (files %2.0f - %2.0f)', e_peak_amp.mean,xbins(e_peak_amp.max_idx),file_start,file_stop);
title(title_str)
saveas(gcf,[store_folder '\Run' run.id '_signalAmplitude_Hist_MM.png'])

%% plot e_peak charge histogram cut included - Only MM
figure
% if (opts_MM.en_filter)
%     hCharge=histogram(e_peak_MM(~isMultipeak),100); 
% else

if use_fixed_binning
    histogram(e_peak_MM,'BinWidth', bin_width_epeak, 'FaceColor','r'); 
    hold on
    hCharge=histogram(e_peak_MM(eventsValidDUT),'BinWidth', bin_width_epeak, 'FaceColor','b');
else
    hCharge=histogram(e_peak_MM(eventsValidDUT),100, 'FaceColor','b');
    hold on
end
 
xbins = hCharge.BinEdges(1:end-1)+hCharge.BinWidth/2;
fit_data = [];
fit_data(1,:)=initial_bin:length(xbins);
fit_data(2,:) = hCharge.Values(fit_data(1,:));
% fit_data(3,:) = yerr;
p0=[];
p0(1) = sum(hCharge.Values)*hCharge.BinWidth;   % normalization factor
p0(2) = 1;
p0(3) = 0.1;
cmd='min; ret';
[p, err, chi] = fminuit('polya_minuit',p0,fit_data(:,1:end),'-b','-c',cmd);
polya_cnt_fit = polya_minuit(p,fit_data(1,:));
plot(xbins(fit_data(1,:)),polya_cnt_fit,'Linewidth',2);
e_peak_charge.mean = sum(polya_cnt_fit.*xbins(fit_data(1,:)))/sum(polya_cnt_fit);
[dummy, e_peak_charge.max_idx] = max(polya_cnt_fit);

xlabel('Peak charge, pC')
ylabel('Events');
ylim([0 1500]);
grid on
legend('wNoise', 'MM','MM fit');
title_str = sprintf('e-peak charge \\mu = %4.4f pC Q_{max} = %4.4f pC', e_peak_charge.mean, xbins(e_peak_charge.max_idx));
title(title_str);
saveas(gcf,[store_folder '\Run' run.id '_signalCharge_Hist_MM.png'])
hold off
%close all;

%% plot all charge histogram cut included - Only MM
if opts_MM.wholeIonTail && opts_MM.en_filter
    figure
    
    if use_fixed_binning
        histogram(all_charge,'BinWidth', bin_width_all, 'FaceColor','r'); 
        hold on
        hAll=histogram(all_charge(eventsValidDUT),'BinWidth', bin_width_all, 'FaceColor','b');
    else
        hAll=histogram(all_charge(eventsValidDUT), 100, 'FaceColor','b'); 
    end

    hold on
    xbins = hAll.BinEdges(1:end-1)+hAll.BinWidth/2;
    fit_data = [];
    fit_data(1,:)=initial_bin:length(xbins);
    fit_data(2,:) = hAll.Values(fit_data(1,:));
    % fit_data(3,:) = yerr;
    p0=[];
    p0(1) = sum(hAll.Values)*hAll.BinWidth;   % normalization factor
    p0(2) = 1;
    p0(3) = 0.1;
    cmd='min; ret';
    [p, err, chi] = fminuit('polya_minuit',p0,fit_data(:,1:end),'-b','-c',cmd);
    polya_cnt_fit = polya_minuit(p,fit_data(1,:));
    plot(xbins(fit_data(1,:)),polya_cnt_fit,'Linewidth',2);
    allCharge.mean = sum(polya_cnt_fit.*xbins(fit_data(1,:)))/sum(polya_cnt_fit);
    [dummy, allCharge.max_idx] = max(polya_cnt_fit);

    xlabel('Peak charge, pC')
    ylabel('Events');
    ylim([0 1500]);
%     ylim([0 1000]);
%     xlim([0 4]);
    grid on
    legend('wNoise', 'MM','MM fit');
    title_str = sprintf('All charge \\mu = %4.4f pC Q_{max} = %4.4f pC', allCharge.mean, xbins(allCharge.max_idx));
    title(title_str);
    saveas(gcf,[store_folder '\Run' run.id '_signalAll_Hist_MM.png'])
    hold off
    %close all;
    
%     parfor k=1:length(MM_data)
%         if (all_charge(k) < 4.5)% && all_charge(k) > 32)
%             figure
%             plot(-1*MM_data(k).waveformY);
%             hold on
%             plot(MM_data(k).sig.smoothed.y, 'r-');
%             scatter(MM_data(k).sig.startpoint.idx, MM_data(k).sig.startpoint.y, 'g', 'filled');
%             scatter(MM_data(k).sig.endpoint.idx, MM_data(k).sig.endpoint.y, 'b', 'filled');
%             scatter(MM_data(k).sig.e_peak_end.idx, MM_data(k).sig.e_peak_end.y, 'r', 'filled');
%             saveas(gcf,[store_folder '\Noise' int2str(k) '.png']);
%             hold off
%             close all;
%         end
%     end
end

%% rise time analysis
riseTimeCut = abs(riseTime)<3;

figure
h=histogram(riseTime(eventsValidDUT&riseTimeCut),100);
hold on
xbins = h.BinEdges(1:end-1)+h.BinWidth/2;

xlabel('Rise time (ns)')
ylabel('Events');
grid on
xlim([0 2]);
%legend('MM','MM fit');
title_str = sprintf('Signal rise time mean = %4.4f ns', mean(riseTime(eventsValidDUT&riseTimeCut)));
title(title_str)
saveas(gcf,[store_folder '\Run' run.id '_riseTimeHist.png'])


%% calculate mean signal 
% lengthMMData = length(MM_data);
% 
% numPoints = length(MM_data(1).waveformY);
% 
% Y_all = NaN(lengthMMData, numPoints);
% 
% %in waveformX there's the virtual time vector (in seconds)
% figure
% for i = 1:lengthMMData
%     plot(MM_data(i).waveformX*10^9,MM_data(i).waveformY); %needs to be corrected for trigger offset??
%     Y_all(i,:) = MM_data(i).waveformY(:)';
% %    hold on;
% end
% 
% 
% Y_avg1 = mean(Y_all,1,'omitnan');
% 
% plot(MM_data(1).waveformX*10^9,Y_avg1);
% 
% xlim([0 200]);
% xlabel('Time, ns')
% ylabel('Amplitude, V')
% title_str = sprintf('Mean of signals');
% title(title_str)
% 
% 
% saveas(gcf,[store_folder '\Run' run.id '_mean_signals.png'])
% close all;

% %% correlations between variables
% % correlation between e_peak_charge and peak amplitude
% figure;
% histogram2(e_peak_MM, MM_maxy, 'DisplayStyle', 'tile', 'Normalization', 'count', 'NumBins', [100, 100]);
% xlabel('e\_peak\_MM');
% ylabel('MM\_maxy');
% title('2D Histogram of e\_peak\_MM vs MM\_maxy');
% xlim([0, 20]);
% ylim([0, 0.15]);
% colorbar;
% 
% saveas(gcf,[store_folder '\Run' run.id '_charge_ampli_2d.png'])
% close all;
% 
% % correlation between e_peak_charge and peak width
% figure;
% histogram2(e_peak_MM, e_peak_width, 'DisplayStyle', 'tile', 'Normalization', 'count', 'NumBins', [100, 100]);
% xlabel('e\_peak\_MM');
% ylabel('e\_peak\_width');
% title('2D Histogram of e\_peak\_MM vs e\_peak\_width');
% colorbar;
% xlim([0, 20]);
% ylim([0, 15]);
% saveas(gcf,[store_folder '\Run' run.id '_charge_peakWidth_2d.png'])
% close all;
% 
% if opts_MM.en_filter
%     % correlation between e_peak_charge and ion tail
%     figure;
%     histogram2(e_peak_MM, ion_tail, 'DisplayStyle', 'tile', 'Normalization', 'count', 'NumBins', [100, 100]);
%     xlabel('e\_peak\_MM');
%     ylabel('ionTail');
%     title('2D Histogram of e\_peak\_MM vs ionTail');
%     colorbar;
%     xlim([0, 20]);
%     ylim([0, 0.015]);
%     
%     saveas(gcf,[store_folder '\Run' run.id '_charge_ionTail_2d.png'])
%     close all;
%     
%     % correlation between ion tail and peak amplitude
%     figure;
%     histogram2(ion_tail, MM_maxy, 'DisplayStyle', 'tile', 'Normalization', 'count', 'NumBins', [100, 100]);
%     xlabel('ionTail');
%     ylabel('MM\_maxy');
%     title('2D Histogram of ionTail vs MM\_maxy');
%     colorbar;
%     xlim([0, 0.015]);
%     ylim([0, 0.15]);
%     saveas(gcf,[store_folder '\Run' run.id '_ionTail_ampli_2d.png'])
%     close all;
% 
%     % correlation between e_peak_charge and peak amplitude with multipeak excluded
%     figure;
%     histogram2(e_peak_MM(~isMultipeak), MM_maxy(~isMultipeak), 'DisplayStyle', 'tile', 'Normalization', 'count', 'NumBins', [100, 100]);
%     xlabel('e\_peak\_MM');
%     ylabel('MM\_maxy');
%     title('2D Histogram of e\_peak\_MM vs MM\_maxy (multipeak excluded)');
%     colorbar;
%     xlim([0, 20]);
%     ylim([0, 0.15]);
%     saveas(gcf,[store_folder '\Run' run.id '_charge_ampli_2d_selected.png'])
%     close all;
% 
%     % correlation between e_peak_charge and ion tail (no multipeak)
%     figure;
%     histogram2(e_peak_MM(~isMultipeak), ion_tail(~isMultipeak), 'DisplayStyle', 'tile', 'Normalization', 'count', 'NumBins', [100, 100]);
%     xlabel('e\_peak\_MM');
%     ylabel('ionTail');
%     title('2D Histogram of e\_peak\_MM vs ionTail (multipeak excluded)');
%     colorbar;
%     xlim([0, 20]);
%     ylim([0, 0.015]);
%     
%     saveas(gcf,[store_folder '\Run' run.id '_charge_ionTail_2d_selected.png'])
%     close all;
% 
%     % correlation between e_peak_charge and ion tail (only multipeak)
%     figure;
%     histogram2(e_peak_MM(isMultipeak), ion_tail(isMultipeak), 'DisplayStyle', 'tile', 'Normalization', 'count', 'NumBins', [100, 100]);
%     xlabel('e\_peak\_MM');
%     ylabel('ionTail');
%     title('2D Histogram of e\_peak\_MM vs ionTail (only multipeak)');
%     colorbar;
%     xlim([0, 20]);
%     ylim([0, 0.015]);
% 
%     saveas(gcf,[store_folder '\Run' run.id '_charge_ionTail_2d_multipeak.png'])
%     close all;
% 
%     % correlation between e_peak_charge and peak width
%     figure;
%     histogram2(e_peak_MM(~isMultipeak), e_peak_width(~isMultipeak), 'DisplayStyle', 'tile', 'Normalization', 'count', 'NumBins', [100, 100]);
%     xlabel('e\_peak\_MM');
%     ylabel('e\_peak\_width');
%     title('2D Histogram of e\_peak\_MM vs e\_peak\_width (no multipeak)');
%     colorbar;
%     xlim([0, 20]);
%     ylim([0, 15]);
%     saveas(gcf,[store_folder '\Run' run.id '_charge_peakWidth_2d_selected.png'])
%     close all;
%     
% 
%     % correlation between e_peak_charge and peak width
%     figure;
%     histogram2(e_peak_MM(isMultipeak), e_peak_width(isMultipeak), 'DisplayStyle', 'tile', 'Normalization', 'count', 'NumBins', [100, 100]);
%     xlabel('e\_peak\_MM');
%     ylabel('e\_peak\_width');
%     title('2D Histogram of e\_peak\_MM vs e\_peak\_width (only multipeak)');
%     colorbar;
%     xlim([0, 20]);
%     ylim([0, 15]);
% 
%     saveas(gcf,[store_folder '\Run' run.id '_charge_peakWidth_2d_multipeak.png'])
%     close all;
% 
%     % correlation between ion tail and peak amplitude
%     figure;
%     histogram2(ion_tail(~isMultipeak), MM_maxy(~isMultipeak), 'DisplayStyle', 'tile', 'Normalization', 'count', 'NumBins', [100, 100]);
%     xlabel('ionTail');
%     ylabel('MM\_maxy');
%     title('2D Histogram of ionTail vs MM\_maxy (no multipeak)');
%     colorbar;
%     xlim([0, 0.015]);
%     ylim([0, 0.15]);
%     saveas(gcf,[store_folder '\Run' run.id '_ionTail_ampli_2d_selected.png'])
%     close all;
% 
%     % correlation between ion tail and peak amplitude
%     figure;
%     histogram2(ion_tail(isMultipeak), MM_maxy(isMultipeak), 'DisplayStyle', 'tile', 'Normalization', 'count', 'NumBins', [100, 100]);
%     xlabel('ionTail');
%     ylabel('MM\_maxy');
%     title('2D Histogram of ionTail vs MM\_maxy (only multipeak)');
%     colorbar;
%     xlim([0, 0.015]);
%     ylim([0, 0.15]);
%     saveas(gcf,[store_folder '\Run' run.id '_ionTail_ampli_2d_multipeak.png'])
%     close all;

% end

% % % % % % % % % % % % % % % % %%%

% % % % % % % % % % % % % % % 

if shouldSaveMat
    if opts_MM.en_filter
            save([mat_folder '\Run' run.id '-' run.oscilloscope '.mat'], 'MM_data', 'MM_maxy', 'e_peak_MM', 'isMultipeak');
    else
            save([mat_folder '\Run' run.id '-' run.oscilloscope '.mat'], 'MM_data', 'MM_maxy', 'e_peak_MM');
    end
end