function out = process_signal_minuit(x,y,options,run,shouldSave,store_folder)
% function for processing the waveform to identify timing parameters

    if(options.invert)
        y =-y;        % invert signal from MM 
    end
    
    
    %-- sig.blrms = std(y(10:100));
    %-- sig.blavg = mean(y(10:100));

    %for short acq
        %sig.blrms = std(y(2:20));
    %sig.blavg = mean(y(2:20));

    %% correct for DC offset
    %-- y = y - sig.blavg;
  
    
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

    %% calculate mean value and noise
    if(sig.max.idx > options.t_prerms)
        sig.blrms = std(y(options.t_dead:sig.max.idx-options.t_prerms));
        sig.blavg = mean(y(options.t_dead:sig.max.idx-options.t_prerms));
        if(options.en_plot || shouldSave)
            scatter(x([options.t_dead sig.max.idx-options.t_prerms]),y([options.t_dead sig.max.idx-options.t_prerms]));
        end
    else
        out.fail = 1;
        return 
    end
    %% dont process signals with statistically unsignificant peaks
    if (sig.max.y < 3*sig.blrms)
        out.fail = 2;
        return
    end
    %% correct for DC offset
    y = y - sig.blavg;

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
    sig.endpoint.idx = sig.max.idx:sig.length;
    for i=sig.max.idx:sig.length
        if(y(i) - sig.blrms < 0)
            sig.endpoint.idx = i;
            break;
        end
    end
    sig.endpoint.y = y(i);
    sig.endpoint.x = x(i);

    if(options.en_plot || shouldSave)
        scatter(sig.startpoint.x, sig.startpoint.y);
        scatter(sig.endpoint.x, sig.endpoint.y);
    end

    %% find electron peak point
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
    % finish first subplot
    if(options.en_plot || shouldSave)
        scatter(sig.e_peak_end.x, sig.e_peak_end.y);
        title('RAW signal');
        ylabel('Voltage, V');
        xlabel('Time, ns');
        grid on;
        hold off;
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

    if(isempty(cfd.zci))
        %out.fail=4;
        %return;
        cfd.time = 9999.9;
    else
        cfd.x(1) = x(cfd.zci(1)+cfd.min.idx-1);
        cfd.x(2) = x(cfd.zci(1)+cfd.min.idx);
        cfd.y(1) = cfd.wav(cfd.zci(1)+cfd.min.idx-1);
        cfd.y(2) = cfd.wav(cfd.zci(1)+cfd.min.idx);
        cfd.time = lin_interp(cfd.x,cfd.y,0);
    end


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
        return;
    end
    
    sigmoid.end = sig.max.idx + 5;
    if (sigmoid.end > sig.length)   
        sigmoid.end = sig.length;
        out.fail = 6;
        return;
    end
    sigmoid.npoints = sigmoid.end-sigmoid.start+1;
    if(sigmoid.npoints>100 || sigmoid.npoints<=0)
        out.fail = 7;
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
    
    evalc("[p, err, chi] = fminuit('ferdirac5_minuit',p0,data(:,1:end-5),'-b','-c',cmd);");
    
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
        hold off;
        drawnow;
        p;
        pause(0.5);
    end
    
    % calculate timing from sigmoid
    sigmoid.time20=0;
    sigmoid.sevals=0;
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
        saveas(gcf,[store_folder '\signal' int2str(run.savedSignals) '.png']);
            pause(1);
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




    
