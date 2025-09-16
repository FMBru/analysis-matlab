%A strong filtering of noise is applied to have a clean signal, so that
%is easier to define e-peak, identify multipeak signal and redefine startpoint and endpoint

function sig = filteredAnalysis(x, y,sig, enIon)
    smoothed.y = sgolayfilt(y, 3, 85);
    plot(x,smoothed.y,'r');
    [smoothed.max.y, smoothed.max.idx] = max(smoothed.y);

    %to define the e-peak end the ion tail avg value is calculated and
    %an intersection between the identified the negative slope e-peak
    %line and the ion-tail avg line is made

    %To find if this signal is a multipeak (multiavalanche) signal and
    %to the define the starting point of the ion-tail (rough estimate
    %of e-peak-end) local max and min are searched
    
    sig.fail = 0;

    if (smoothed.max.idx+300 > sig.length || sig.max.idx+300 > sig.length)
        sig.fail = 7;
        return
    elseif (smoothed.max.idx+700 > sig.length && enIon)
        sig.fail = 7;
        return    
    end
    
    if (smoothed.max.y < 2*sig.blrms)
        sig.fail = 8;
        return
    end

    %if another local max is found to be >thr_value probably it's a multipeak signal
    % in this case the local minium is searched after the last peak 

%     if enIon
%         rms_filtered = std(smoothed.y(smoothed.max.idx+(200:700)));
%         mean_filtered = mean(smoothed.y(smoothed.max.idx+(200:700)));
%     else
%         rms_filtered = std(smoothed.y(smoothed.max.idx+200:sig.length));
%         mean_filtered = mean(smoothed.y(smoothed.max.idx+200:sig.length));
%     end
    
    rms_filtered = std(smoothed.y(smoothed.max.idx+(200:700)));
    mean_filtered = mean(smoothed.y(smoothed.max.idx+(200:700)));
    sig.is_multipeak = 0;


    peaks_interval = sig.startpoint.idx:sig.max.idx+300; %300 is 15% of signal length

    if(smoothed.max.y < mean_filtered + 3*rms_filtered)
        sig.fail = 9;
        return
    end

    if(sum(smoothed.y(peaks_interval) > mean_filtered + 3*rms_filtered) < 1)
        sig.fail = 9;
        return
    end

    [ smoothed.peaks, smoothed.peakidx ] = findpeaks(smoothed.y(peaks_interval), 'MinPeakDistance', 15, 'MinPeakHeight', mean_filtered + 3*rms_filtered);
    
    %this could mean that there's no significant peak -> noise
    if (length(smoothed.peaks) < 1)
        sig.fail = 9;
%         plot(x,y, 'b');
%         hold on
%         plot(x,smoothed.y,'r');
%         hold off
%         close all;
        return
    end
    smoothed.peakidx = peaks_interval(1) + smoothed.peakidx; 
    %scatter(x(smoothed.peakidx), smoothed.peaks, 'g', 'filled');
    
    %threshold value is higher if the signal is lower because it could be
    %that noise peaks are confused for real ones

    thr_value = 0.4*smoothed.max.y;
    if (smoothed.max.y < 0.01)
        thr_value = 0.6*smoothed.max.y;
    end 

    multipeak_boolean = smoothed.peaks > thr_value;
    if (sum(multipeak_boolean)>1)
        sig.is_multipeak = 1;
    end

    if (sig.is_multipeak)
        temp_min = findLocalMin(smoothed.y, smoothed.peakidx(end)+5, smoothed.peakidx(end)+100); %search for the first local min (5 and 100 are magical numbers)
        temp_max.idx = smoothed.peakidx(end);
    else
        temp_min = findLocalMin(smoothed.y, smoothed.max.idx+5, smoothed.max.idx+100); %search for the first local min (5 and 100 are magical numbers)
        temp_max.idx = smoothed.max.idx;
    end
    
    
    if enIon
        sig.ionTailAvg = mean(smoothed.y(temp_min.idx+(1:500))); %I'm using 500 that is almost 50 ns after
    else
        % here I'm assuming (as observed) that the ion-tail continues till the end of the signal
        sig.ionTailAvg = mean(smoothed.y(temp_min.idx:sig.length - 5));
    end
    

    %searching for a line describing the negative slope part of e_peak
    %(between 20% and 70% points of the interval between max and min)
    line_start = round(temp_max.idx + (temp_min.idx-temp_max.idx)*0.2);
    line_end = round(temp_max.idx + (temp_min.idx-temp_max.idx)*0.7);
    line_fit_range = [line_start line_end];

    line_coeff = polyfit(line_fit_range, smoothed.y(line_fit_range), 1);
    e_peak_line = polyval(line_coeff, line_start:line_end);
%     plot(x(line_start:line_end), e_peak_line, 'g-');
%     plot(x(temp_min.idx:sig.length), sig.ionTailAvg*ones(size(temp_min.idx:sig.length)), 'b-');
    e_peak_slope = line_coeff(1);
    e_peak_interc = line_coeff(2);

    sig.e_peak_end.idx = round((sig.ionTailAvg - e_peak_interc) / e_peak_slope);
    if(sig.e_peak_end.idx < sig.max.idx)
        sig.e_peak_end.idx = temp_min.idx;
    end
    
    if  (sig.e_peak_end.idx > sig.length)
        sig.fail = 10;
        sig.e_peak_end.idx = sig.length;
        return
    end
    sig.e_peak_end.x = x(sig.e_peak_end.idx);
    sig.e_peak_end.y = y(sig.e_peak_end.idx);
    sig.smoothed = smoothed;
    sig.multipeak = multipeak_boolean;

    %% find new start point
    sig.startpoint.idx = 1;
    for i=smoothed.peakidx(1):-1:1
        if(smoothed.y(i)-rms_filtered < 0)   
            sig.startpoint.idx = i;
            break;
        end
    end
    sig.startpoint.x = x(i);
    sig.startpoint.y = y(i);

    %% find new end point
    sig.endpoint.idx = sig.length;  %I added this as condition when the endpoint is not found
    for i=sig.length:-1:smoothed.peakidx(end)
        if(smoothed.y(i) - rms_filtered > 0)
            sig.endpoint.idx = i;
            break;
        end
    end
    
%     sig.endpoint.idx = sig.length;  %I added this as condition when the endpoint is not found
%     for i=smoothed.peakidx(end):sig.length
%         if(smoothed.y(i) - rms_filtered < 0)
%             sig.endpoint.idx = i;
%             break;
%         end
%     end
    
    sig.endpoint.y = y(i);
    sig.endpoint.x = x(i);
    return
end

%% finds the first local min in a certain range
% points = y vector points
% rangeStart (End) = index of range starting (ending) point
% local_min = structure with a .y and .idx


function local_min = findLocalMin(points, rangeStart, rangeEnd)
    local_min.y = points(rangeStart);
    local_min.idx = rangeStart;
    %from the first point it searches for a local min and it breaks when it finds the first one 
    for i=rangeStart:rangeEnd
        %if a certain point is greater than the currently identified
        %local_min, break the loop
        if (points(i) > 1.1*local_min.y) %1.1 to avoid little bumps problems
            break
        elseif (points(i) < local_min.y)
            local_min.y = points(i);
            local_min.idx = i;
        end
    end
    return
end