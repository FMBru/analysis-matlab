% created by Francesco - 06/07/25

close all;
enableAmplitudePlot = 1;
enableChargePlot = 0;

if exist('dataFolder','var') == 0
    dataFolder = '\\eosproject-smb\eos\project\p\picosec\lab\Francesco\UV_LED_labtest\metallic_singlechannel\sealedMode\MatFile';
    strBoundaryCond = '';
    strSubBoundary = '';
    strTypeOfAnalysis = '';
    strSavePictureName = 'comparisonSealedMode4days';
    tic;
    data1 = load([dataFolder strBoundaryCond strSubBoundary strTypeOfAnalysis '\' 'Run0807_Ti_C100_A275_C470_newCSA_spectrum_highStat']);
    data1.name = '1st day';
    files{1} = data1;

    data2 = load([dataFolder strBoundaryCond strSubBoundary strTypeOfAnalysis '\' 'Run0907_Ti_C100_A275_C470_newCSA_spectrum_highStat2']);
    data2.name = '2nd day';
    files{2} = data2;

    data3 = load([dataFolder strBoundaryCond strSubBoundary strTypeOfAnalysis '\' 'Run1007_Ti_C100_A275_C470_newCSA_spectrum_highStat']);
    data3.name = '3rd day';
    files{3} = data3;

    data4 = load([dataFolder strBoundaryCond strSubBoundary strTypeOfAnalysis '\' 'Run1107_Ti_C100_A275_C470_newCSA_spectrum_highStat']);
    data4.name = '4th day';
    files{4} = data4;
    
    toc

    mat_folder = ['\\eosproject-smb\eos\project\p\picosec\lab\Francesco\UV_LED_labtest\metallic_singlechannel\sealedMode\Results' strBoundaryCond strSubBoundary strTypeOfAnalysis];% add_str];

    if exist(mat_folder, 'dir') == 0
        mkdir(mat_folder);
    end
end

colors = lines(length(files));  % distinct colors for each file

if enableChargePlot
    figure(1); hold on; title('Histograms Electron Peak Charge');
    xlabel('e peak charge (pC)'); ylabel('Relative frequency');
end

if enableAmplitudePlot
    figure(2); hold on; title('Histograms Peak Amplitude');
    xlabel('Amplitude (V)'); ylabel('Relative frequency');
end

legendEntries1 = {};
legendEntries2 = {};
bin_width_e_peak = 0.4;
bin_width_maxy = 0.002;

for k = 1:length(files)
    
    % Existence check
    if isfield(files{k}, 'e_peak_MM') && enableChargePlot
        % Extract data
        e_peak = files{k}.e_peak_MM;
        
        % to normalize histogram e_peak
        figure(3)
%         bin_width_e_peak = 0.2 * k;
        edges_e_peak = 0:bin_width_e_peak:(max(e_peak)+bin_width_e_peak);
        counts_e_peak = histcounts(e_peak, edges_e_peak);
        norm_e_peak = 100 / max(counts_e_peak(5:end)); %normalization relative to max bin (after noise)
        close;
        
        % Plot e_peak_MM
        figure(1);
        counts_e_peak_scaled = counts_e_peak * norm_e_peak;
        histogram('BinEdges', edges_e_peak, 'BinCounts', counts_e_peak_scaled, 'DisplayStyle', 'stairs', 'EdgeColor', colors(k, :));

        legendEntries1{end+1} = files{k}.name;
    else
        fprintf('ATTENZIONE: %s non contiene e_peak_MM\n', files{k}.name);
    end
    
    % Existence check
    if enableAmplitudePlot && isfield(files{k}, 'MM_maxy')
        % Extract data
        maxy = files{k}.MM_maxy;
        
        % to normalize histogram amplitude
        figure(4)
        edges_maxy = 0:bin_width_maxy:(max(maxy)+bin_width_maxy);
        counts_maxy = histcounts(maxy, edges_maxy);
        norm_maxy = 100 / max(counts_maxy(3:end)); %normalization relative to max bin (after noise)
        close;
        
        % Plot MM_maxy
        figure(2);
        counts_maxy_scaled = counts_maxy * norm_maxy;
        histogram('BinEdges', edges_maxy, 'BinCounts', counts_maxy_scaled, 'DisplayStyle', 'stairs', 'EdgeColor', colors(k, :));

        legendEntries2{end+1} = files{k}.name;
    else
        fprintf('ATTENZIONE: %s non contiene MM_maxy\n', files{k}.name);
    end
end

if enableChargePlot
    figure(1); legend(legendEntries1, 'Interpreter', 'none'); xlim([-1.5 30]); ylim([0 120]);
    saveas(figure(1),[mat_folder '\' strSavePictureName '_Charge.png']);
end

if enableAmplitudePlot
    figure(2); legend(legendEntries2, 'Interpreter', 'none'); ylim([0 120]);
    saveas(figure(2),[mat_folder '\' strSavePictureName '_Amplitude.png']);
end

close all;
