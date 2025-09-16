%updated by Francesco (summer student) in 30/06/25

function out = process_signal_minuit(x,y,options,run,shouldSave,store_folder,eventID)
% function for processing the waveform to identify timing parameters

    if(options.invert)
        y =-y;        % invert signal from MM 
    end
    
    %magical interval where there's only noise (probably) -> precursor
    sig.blrms = std(y(10:50));
    sig.blavg = mean(y(10:50));

    %% correct for DC offset
    y = y - sig.blavg;
    
    x = x*1e9;    % convert time to nS

    sig.length = length(x);
    
    if(options.en_plot || shouldSave)
        subplot(2,2,[1 2]);
        plot(x,y);
        hold on
    end

    out.fail=0;
    %% find basic stats from the signal
    [sig.max.y, sig.max.idx] = max(y);
    [sig.min.y, sig.min.idx] = min(y);
    
    if (sig.max.y >= options.saturation_value)
        out.fail = 1; %if the max is at the saturation value it could be "truncated"
        if shouldSave
            saveImage(store_folder, eventID);
        end 
        return
    end

    %% calculate mean value and noise
    %t_prerms is a constant to be used to select a time range for evaluating the noise rms and avg
    %This time range starts with t_dead (first sampling point) and ends...
    %...with a point just before the e-peak: to find this ending point (sig.max.idx - t_prerms) is used... 
    %... or this ending point is evaluated as sig.startpoint (see later) with a raw estimation of sig.blrms (now I'm using this last one)

    if(sig.max.idx > options.t_prerms)

        %                   ---- sig.max.idx - t_prerms mode -----
        if options.en_noiseRejection == 0
            sig.blrms = std(y(options.t_dead:sig.max.idx-options.t_prerms)); 
            sig.blavg = mean(y(options.t_dead:sig.max.idx-options.t_prerms));
            if(options.en_plot || shouldSave)
                scatter(x([options.t_dead sig.max.idx-options.t_prerms]),y([options.t_dead sig.max.idx-options.t_prerms]));
            end
        else
            %                   ---- sig.upperrms.idx mode -----
            sig.upperrms.idx = 1;
            for i=sig.max.idx:-1:1
                if(y(i)-sig.blrms < 0)   
                    sig.upperrms.idx = i;
                    break;
                end
            end
            sig.blrms = std(y(options.t_dead:sig.upperrms.idx)); 
            sig.blavg = mean(y(options.t_dead:sig.upperrms.idx));
            if(options.en_plot || shouldSave)
                scatter(x([options.t_dead sig.upperrms.idx]),y([options.t_dead sig.upperrms.idx]));
            end
        end
        
    else
        out.fail = 1; %if the max is before t_prerms probably it's noise
        if shouldSave
            saveImage(store_folder, eventID);
        end 
        return
    end
    %% dont process signals with statistically unsignificant peaks
    if (sig.max.y < 3*sig.blrms)
        out.fail = 2;
        if shouldSave
            saveImage(store_folder, eventID); %save image before returning the function
        end
        return
    end
    %% correct for DC offset
    y = y - sig.blavg;
    %sig.max.y should be updated

    %% find start point
    sig.startpoint.idx = 1;
    for i=sig.max.idx:-1:1
        if(y(i)-sig.blrms < 0)   
            sig.startpoint.idx = i;
            break;
        end
    end
    sig.startpoint.x = x(i);
    sig.startpoint.y = y(i);

    %% find end point
    sig.endpoint.idx = sig.length;  %I added this as condition when the endpoint is not found
    for i=sig.max.idx:sig.length
        if(y(i) - sig.blrms < 0)
            sig.endpoint.idx = i;
            break;
        end
    end
     
    sig.endpoint.y = y(i);
    sig.endpoint.x = x(i);

    %this startpoint and endpoint are raw estimation of them

    
    %to discriminate signal from noise, the average value between the
    %starpoint and the endpoint is evaluated: since between those points
    %there should be the e-peak followed by the ion-tail the expected avg
    %is considerably more than 0. The threshold (empirically found) is 2.5 rms

    if out.fail == 0 && options.en_noiseRejection
        yavg = mean(y(sig.startpoint.idx:sig.endpoint.idx));
        if (yavg < 2.5*sig.blrms)
            disp("Noise Rejection");
            out.fail = 8;
            if shouldSave
                saveImage(store_folder, eventID);
            else
                plot(x,y);
                mkdir([store_folder '\rejected']);
                saveImage([store_folder '\rejected'], eventID);
            end
            return
        end 
    end


    %-------------------------------------------------------------------------------
    %A strong filtering of noise is applied to have a clean signal, so that
    %is easier to define e-peak, identify multipeak signal and redefine startpoint and endpoint
    
    sig.is_multipeak = 0;
    if options.en_filter
        sig = filteredAnalysis(x, y, sig);
        out.fail = sig.fail;
        if (out.fail > 0)
            if shouldSave
            saveImage(store_folder, eventID);
            end
            return
        end
    else

        % find electron peak point
        sig.e_peak_end.y = 11111111111.;
        sig.e_peak_end.idx = sig.max.idx;
    
        j = sig.max.idx + 180; % set to 60 in original code
        if( j > sig.length - 1 ) 
            j = sig.length - 2; 
        end
    
        % algorithm taken from c code
        for i = sig.max.idx+3:j   
            if(y(i) < sig.e_peak_end.y)            
                sig.e_peak_end.y = y(i);
                sig.e_peak_end.idx = i;
                if( i+1 > sig.length)
                    out.fail=3;
                    if shouldSave
                        saveImage(store_folder, eventID);
                    end
                    return
                else
                    if(y(i+1) < 0)
                        break;
                    end
                end
            end
        end
        sig.e_peak_end.x = x(sig.e_peak_end.idx);
        sig.e_peak_end.y = y(sig.e_peak_end.idx);
    end
    
    if(options.en_plot || shouldSave)
        scatter(sig.startpoint.x, sig.startpoint.y);
        scatter(sig.endpoint.x, sig.endpoint.y);
    end


    % finish first subplot
    if(options.en_plot || shouldSave)
        scatter(sig.e_peak_end.x, sig.e_peak_end.y,'blue', 'filled');        
        if out.fail == 0
            if sig.is_multipeak
                status_str = 'Multipeak';
            else
                status_str = 'Good';
            end
        else
            status_str = 'NOISE';
        end
        title_str = ['RAW signal - ', status_str];
        title(title_str);
        % title('RAW signal');
        ylabel('Voltage, V');
        xlabel('Time, ns');
        grid on;
        if options.en_filter
            legend('RAW', 'RMS range', 'Filtered', 'locals','startpoint', 'endpoint', 'e-peak end');
        else
            legend('RAW', 'RMS range', 'startpoint', 'endpoint', 'e-peak end');
        end
        hold off;
    end
    
    %to evaluate if there's a way to discriminate signal with the peak width
    sig.e_peak_width = sig.e_peak_end.x-sig.startpoint.x;

    if sig.e_peak_width < 0
        out.fail=11;
        if shouldSave
            saveImage(store_folder, eventID);
        end
        return
    end

    %if a peak is too wide, it could be noise or multipeak
    % the value 7 was found looking at the correlation between e_peak
    % charge and width for single peaks and for multipeaks
    
    if options.shouldExcludeWidePeaks && sig.e_peak_width > 7
        out.fail=12;
        if shouldSave
            saveImage(store_folder, eventID);
        end
        return      
        %sig.is_multipeak = 1;
    end
    
    %% Calculate charges
    sig.k_V2C = options.Ts*1e12/options.Rin;  % convert to picoColoumbs

    sig.charge.lead_edge = sig.k_V2C * sum(y(sig.startpoint.idx:sig.max.idx));
    sig.charge.e_peak = sig.k_V2C * sum(y(sig.startpoint.idx:sig.e_peak_end.idx));
    sig.charge.all = sig.k_V2C * sum(y(sig.startpoint.idx:sig.endpoint.idx));
    %sig.charge.full_int = sig.k_V2C * sum(y(sig.startpoint.idx:sig.endpoint.idx));

    %% calculate CFD timing
    cfd.wav = y(1:end-3)-0.8*y(4:end);      % superimpose with shifted signal
    [cfd.max.val, cfd.max.idx] = max(cfd.wav);
    [cfd.min.val, cfd.min.idx] = min(cfd.wav);
    cfd.zci = find(diff(sign(cfd.wav(cfd.min.idx:cfd.max.idx))));

%     if(isempty(cfd.zci))
%         out.fail=4;
%         return;
%     else
%         cfd.x(1) = x(cfd.zci(1)+cfd.min.idx-1);
%         cfd.x(2) = x(cfd.zci(1)+cfd.min.idx);
%         cfd.y(1) = cfd.wav(cfd.zci(1)+cfd.min.idx-1);
%         cfd.y(2) = cfd.wav(cfd.zci(1)+cfd.min.idx);
%         cfd.time = lin_interp(cfd.x,cfd.y,0);
%     end
       cfd.time = 0; %not using CFD, assign dummy value

%     if(options.en_plot)
%         subplot(2,2,3);
%         plot(x(1:end-3),cfd.wav)
%         hold on
%         plot(x,y)
%         scatter(x(cfd.zci+cfd.min.idx),cfd.wav(cfd.zci+cfd.min.idx));
%         xlim([x(sig.startpoint.idx)-1, x(sig.startpoint.idx)+20]);
%         title('CFD result');
%         ylabel('Voltage, V');
%         xlabel('Time, ns');
%         grid on;
%         hold off;
%     end
    
    cfd.wav=[];  % dont save waveform to
    
    %% calculate sigmoid fit to the leading edge
    sigmoid.start = sig.startpoint.idx - 5;  % 5 samples in original code
    if (sigmoid.start < 1)
        sigmoid.start = 1;
        out.fail = 5;
        if shouldSave
            saveImage(store_folder, eventID);
        end
        return
    end
    
    %this part can be improved
    if options.en_filter && sig.is_multipeak
        sigmoid.end = sig.smoothed.peakidx(1) + 5; %if there's more than one peak, try to take the timing info from the first one
    else
        sigmoid.end = sig.max.idx + 5;
    end

    if (sigmoid.end > sig.length)   
        sigmoid.end = sig.length;
        out.fail = 6;
        if shouldSave
            saveImage(store_folder, eventID);
        end        
        return
    end
    sigmoid.npoints = sigmoid.end-sigmoid.start+1;
    if(sigmoid.npoints>100 || sigmoid.npoints<=0)
        out.fail = 7;
        if shouldSave
            saveImage(store_folder, eventID);
        end
        return
    end
    
    xfit = x(sigmoid.start:sigmoid.end);
    yfit = y(sigmoid.start:sigmoid.end);
    yerr = ones(size(xfit))*sig.blrms;
    
    % fill initial parameters
    p0(1) = sig.max.y;
    p0(2) = 0.5*(x(sigmoid.end)+x(sigmoid.start));
    p0(3) = 5.0/(x(sigmoid.end)-x(sigmoid.start));
    p0(4) = 1.0;
    p0(5) = 0.0;
    
    data(1,:) = xfit;
    data(2,:) = yfit;
    data(3,:) = yerr;
    
    cmd='min; ret';
    
%    evalc("[p, err, chi] = fminuit('ferdirac5_minuit',p0,data(:,1:end-5),'-b','-c',cmd);");

    data(2,end-4:end)=data(2,end-5); %improved sigmoid - last 5 points fixed
    %p is where the parameters of the sigmoid are stored
    evalc("[p, err, chi] = fminuit('ferdirac5_minuit',p0,data(:,1:end),'-b','-c',cmd);");

    % call unconstrained optimization
    %optimopts = optimset('Display','off','MaxFunEvals',1000);
    %[p,fval,exitflag,outputs] = fminsearch(@(p)sigmoid_fit(xfit(1:end-5),yfit(1:end-5),p),p0,optimopts);
    sigmoid.p = p;
    sigmoid.err = err;
    sigmoid.chi=chi;
    %sigmoid.exitflag = exitflag;
    %sigmoid.niter = outputs.funcCount;
    %sigmoid.ndf = length(sigmoid.p)-1;
    %sigmoid.chisquare = sum(((yfit(1:end-5)-ferdirac5(xfit(1:end-5),p)).^2)./yfit(1:end-5))./sigmoid.ndf;

    
    % check fit
    if(options.en_plot || shouldSave) 
        subplot(2,2,[3 4]);
        scatter(xfit,yfit);
        hold on;
        plot(xfit,ferdirac5_minuit(p,xfit'));
        title('Sigmoid fit result');
        ylabel('Voltage, V');
        xlabel('Time, ns');
        grid on
        hold on;
        drawnow;
        p;
        %pause(1);
    end
    
    % calculate timing from sigmoid
    sigmoid.time20=0;
    sigmoid.sevals=0; %number of evaluted points
    sigmoid.timepoint=xfit(1);
    
    sigmoid.time10=0;
    sigmoid.sevals10=0;
    sigmoid.timepoint10=xfit(1);

    sigmoid.time90=0;
    sigmoid.sevals90=0;
    sigmoid.timepoint90=xfit(1);

    
    % 20% val search algorithm
    % first go forward
    while(sigmoid.time20 < 0.2*sig.max.y)
        sigmoid.sevals = sigmoid.sevals + 1;
        sigmoid.timepoint = sigmoid.timepoint + 0.01;          % 10ps stepping forward
        sigmoid.time20 = ferdirac5_minuit(p,sigmoid.timepoint);
        if(sigmoid.timepoint>x(sig.max.idx))
            sigmoid.timepoint=0;  % maybe add fail here
            break;
        end   
    end
    % refine backwards
    while(sigmoid.time20 > 0.2*sig.max.y)
        sigmoid.sevals = sigmoid.sevals + 1;
        sigmoid.timepoint = sigmoid.timepoint - 0.0001;        % 0.1ps stepping backwards
        sigmoid.time20 = ferdirac5_minuit(p,sigmoid.timepoint);
        if(sigmoid.timepoint<xfit(1))
            sigmoid.timepoint=0;  % maybe add fail here
            break;
        end  
    end

    % 10% val search algorithm
    % first go forward
    while(sigmoid.time10 < 0.1*sig.max.y)
        sigmoid.sevals10 = sigmoid.sevals10 + 1;
        sigmoid.timepoint10 = sigmoid.timepoint10 + 0.01;          % 10ps stepping forward
        sigmoid.time10 = ferdirac5_minuit(p,sigmoid.timepoint10);
        if(sigmoid.timepoint10>x(sig.max.idx))
            sigmoid.timepoint10=0;  % maybe add fail here
            break;
        end   
    end
    % refine backwards
    while(sigmoid.time10 > 0.1*sig.max.y)
        sigmoid.sevals10 = sigmoid.sevals10 + 1;
        sigmoid.timepoint10 = sigmoid.timepoint10 - 0.0001;        % 0.1ps stepping backwards
        sigmoid.time10 = ferdirac5_minuit(p,sigmoid.timepoint10);
        if(sigmoid.timepoint10<xfit(1))
            sigmoid.timepoint10=0;  % maybe add fail here
            break;
        end  
    end

        % 90% val search algorithm
    % first go forward
    while(sigmoid.time90 < 0.9*sig.max.y)
        sigmoid.sevals90 = sigmoid.sevals90 + 1;
        sigmoid.timepoint90 = sigmoid.timepoint90 + 0.01;          % 10ps stepping forward
        sigmoid.time90 = ferdirac5_minuit(p,sigmoid.timepoint90);
        if(sigmoid.timepoint90>x(sig.max.idx))
            sigmoid.timepoint90=0;  % maybe add fail here
            break;
        end   
    end
    % refine backwards
    while(sigmoid.time90 > 0.9*sig.max.y)
        sigmoid.sevals90 = sigmoid.sevals90 + 1;
        sigmoid.timepoint90 = sigmoid.timepoint90 - 0.0001;        % 0.1ps stepping backwards
        sigmoid.time90 = ferdirac5_minuit(p,sigmoid.timepoint90);
        if(sigmoid.timepoint90<xfit(1))
            sigmoid.timepoint90=0;  % maybe add fail here
            break;
        end  
    end

    if shouldSave

        if (sig.is_multipeak)
            saveas(gcf,[store_folder '\multipeak' int2str(eventID) '.png']);
        else
            saveas(gcf,[store_folder '\signal' int2str(eventID) '.png']);
        end
        %pause(1);
        hold off
        close all;
    
        hold off 
    end

    %% wrap outputs
    out.cfd = cfd;
    out.sig = sig;
    out.sigmoid = sigmoid;

return
end

%% AUX functions below

%% linear interpolation function
function x=lin_interp(xq,yq,y)
    b = (yq(2)-yq(1))/(xq(2)-xq(1));
    a = yq(2) - b*xq(2);
    x = (y-a)/b;
end

%% save images before returning
function saveImage(folder, ID)
    saveas(gcf,[folder '\noise' int2str(ID) '.png']);
    %pause(1);
    hold off
    close all;

    hold off
end


