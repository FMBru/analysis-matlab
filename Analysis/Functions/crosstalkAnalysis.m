function [crosstalk_vector] = crosstalkAnalysis(allpads_table, padArray, pad_centers_array2, geom, store_folder)
%CROSSTALKANALYSIS Summary of this function goes here


full_geom.x_centroid = mean(pad_centers_array2(:,1));
full_geom.y_centroid = mean(pad_centers_array2(:,2));
full_geom.x_window_length = ceil(length(padArray)/2) * geom.pad_window_length;
full_geom.y_window_length = ceil(length(padArray)/2) * geom.pad_window_length;
full_geom.pad_window_length = geom.pad_window_length;

MCP_ampCut = allpads_table.MCP_amp >= 0.15 & allpads_table.MCP_amp > 0.01*max(allpads_table.MCP_amp) & allpads_table.MCP_amp < 0.99*max(allpads_table.MCP_amp);
SATCut = abs(allpads_table.(['SAT' num2str(padArray(1))]) - median(allpads_table.(['SAT' num2str(padArray(1))]), 'omitnan')) < 10;
riseTimeCut = abs(allpads_table.(['riseTime' num2str(padArray(1))]) - median(allpads_table.(['riseTime' num2str(padArray(1))]), 'omitnan')) < 10;

for i=2:length(padArray)
    SATCut = SATCut & (abs(allpads_table.(['SAT' num2str(padArray(i))]) - median(allpads_table.(['SAT' num2str(padArray(i))]), 'omitnan')) < 10);
    riseTimeCut = riseTimeCut & abs(allpads_table.(['riseTime' num2str(padArray(i))]) - median(allpads_table.(['riseTime' num2str(padArray(i))]), 'omitnan')) < 10;
end

subset = allpads_table(MCP_ampCut & SATCut & riseTimeCut, :);
[full_frame, full_table] = geometricalFrame(subset, full_geom, 1, 'center');

% Compute counts map
[~,~,countMap] = twoDimensionsMap(full_frame, ones(height(full_table),1));
plotMap(full_frame, countMap);
title('Count map full');
close;

for i=1:length(padArray)
    padExcluded = padArray;
    padExcluded(i) = [];

    %amplitude cuts on one pad to see what's on the others
    mmAmpVec = full_table.(['MM_amp' num2str(padArray(i))]);

    highAmpCut = mmAmpVec > 0.2 * max(mmAmpVec, [], 'omitnan');
    %lowAmpCut = mmAmpVec < 0.2 * max(mmAmpVec, [], 'omitnan');

    %geometrical cut on the pad
    pad_x_center = pad_centers_array2(i,1);
    pad_y_center = pad_centers_array2(i,2);

    pad_window_left = pad_x_center - 0.3*full_geom.pad_window_length/2;
    pad_window_right = pad_x_center + 0.3*full_geom.pad_window_length/2;
    pad_window_up = pad_y_center + 0.3*full_geom.pad_window_length/2;
    pad_window_bottom = pad_y_center - 0.3*full_geom.pad_window_length/2;

    padCut = full_table.X > pad_window_left & full_table.X < pad_window_right & full_table.Y > pad_window_bottom & full_table.Y < pad_window_up;

    cross_high_table = full_table(highAmpCut & padCut,:);
    %cross_low_table = full_table(lowAmpCut & padCut,:);
   
    crossCopy = cross_high_table;

    for j=1:length(padExcluded)
        cross_high_table = crossCopy;
        satvar = cross_high_table.(['SAT' num2str(padExcluded(j))]);
        SATfiner = abs(satvar - mean(satvar, 'omitnan')) < 1.5;
        cross_high_table = cross_high_table(SATfiner,:);
        
        %ratio between real signal and crosstalk
        ratio_high = cross_high_table.(['MM_amp' num2str(padArray(i))]) ./ cross_high_table.(['MM_amp' num2str(padExcluded(j))]);
        %ratio_low = cross_low_table.(['MM_amp' num2str(padArray(i))]) ./ cross_low_table.(['MM_amp' num2str(padExcluded(j))]);

        distribution(cross_high_table.(['SAT' num2str(padExcluded(j))]),100,['SAT-highCutOn' num2str(padArray(i))],'ns', store_folder, padExcluded(j), true); %SAT distribution (high Amp)
        %distribution(cross_low_table.(['SAT' num2str(padExcluded(j))]),100,['SAT-lowCutOn' num2str(padArray(i))],'ns', store_folder, padExcluded(j)); %SAT distribution (low Amp)
        distribution(cross_high_table.(['MM_amp' num2str(padExcluded(j))]),100,['Amplitude-highCutOn' num2str(padArray(i))],'V', store_folder, padExcluded(j), false); %amplitude distribution (high Amp)
        %distribution(cross_low_table.(['MM_amp' num2str(padExcluded(j))]),100,['Amplitude-lowCutOn' num2str(padArray(i))],'V', store_folder, padExcluded(j)); %amplitude distribution (low Amp)
        distribution(cross_high_table.(['riseTime' num2str(padExcluded(j))]),100,['riseTime-highCutOn' num2str(padArray(i))],'V', store_folder, padExcluded(j),true); %riseTime distribution (high Amp)
        %distribution(cross_low_table.(['riseTime' num2str(padExcluded(j))]),100,['riseTime-lowCutOn' num2str(padArray(i))],'V', store_folder, padExcluded(j)); %riseTime distribution (low Amp)

        correlationPlots(cross_high_table.(['SAT' num2str(padExcluded(j))]), ratio_high, ['SAT-highCutOn' num2str(padArray(i))], 'ns', ['Ratio-pad' num2str(padArray(i)) 'by-pad' num2str(padExcluded(j))], 'au', store_folder, padExcluded(j), true);
        %correlationPlots( cross_low_table.(['SAT' num2str(padExcluded(j))]), ratio_low,  ['SAT-lowCutOn' num2str(padArray(i))], 'ns', ['Ratio-pad' num2str(padArray(i)) 'by-pad' num2str(padExcluded(j))], 'au', store_folder, padExcluded(j));
    end
end

end

