function out = process_signal_sampic(x,y,options,baseline)

%{

close all
options.t_dead = 1;      % blanking time from start of the recording samples
options.t_prerms = 10  % blanking time before global maximum samples
%opts_MM.Ts = 1/20e9;     % sampling speed
options.Rin=50;          % input impedance 
options.invert = 1;      % is inverted (1-inverted, 0-non-inverted)
options.type=0;          % detector type (0-DUT, 1-REF) for future use
%opts_MM.ch_name = ['C2' run.lecroy_name]; % LeCroy file name format
options.en_plot = 1;     % enable debugging plots
x = ch_mm(:,1);
y = ch_mm(:,2);

%}
Ts = x(2)-x(1);
y= y-baseline;

    if(options.invert)
        y =-y;        % invert signal from MM 
    end

    %x = x - x(1); % start from t=0;
    x = x*1e9;    % convert time to nS    
    sig.length = length(x);
    
    if(options.en_plot)
        close all
        subplot(2,2,[1 2]);
        plot(x,y,'.-');
        hold on
    end

    out.fail=0;  % no failure at begining
    %% find basic stats from the signal
    [sig.max.y, sig.max.idx] = max(y);
    [sig.min.y, sig.min.idx] = min(y);

    %% calculate mean value and noise
    if(sig.max.idx>options.t_prerms)
        sig.blrms = std(y(options.t_dead:sig.max.idx-options.t_prerms));
        sig.blavg = mean(y(options.t_dead:sig.max.idx-options.t_prerms));
        if(options.en_plot)
            scatter(x([options.t_dead sig.max.idx-options.t_prerms]),y([options.t_dead sig.max.idx-options.t_prerms]));
        end
    else
        out.fail = 1;
        return 
    end
    %% dont process signals with statistically unsignificant peaks
    if (sig.max.y < 4*sig.blrms)
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

    if(options.en_plot)
        scatter(sig.startpoint.x, sig.startpoint.y);
        scatter(sig.endpoint.x, sig.endpoint.y);
    end

    %% find electron peak point
    sig.e_peak_end.y = 11111111111.;
    sig.e_peak_end.idx = sig.max.idx;

    j = sig.max.idx + 60;
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
    if(options.en_plot)
        scatter(sig.e_peak_end.x, sig.e_peak_end.y);
        title('RAW signal');
        ylabel('Voltage, V');
        xlabel('Time, ns');
        grid on;
        hold off;
    end

    %% Calculate charges
    sig.k_V2C = Ts*1e12/options.Rin;  % convert to picoColoumbs

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
        out.fail=4;
        return;
    else
        cfd.x(1) = x(cfd.zci(1)+cfd.min.idx-1);
        cfd.x(2) = x(cfd.zci(1)+cfd.min.idx);
        cfd.y(1) = cfd.wav(cfd.zci(1)+cfd.min.idx-1);
        cfd.y(2) = cfd.wav(cfd.zci(1)+cfd.min.idx);
        cfd.time = lin_interp(cfd.x,cfd.y,0);
    end


    if(options.en_plot)
        subplot(2,2,3);
        plot(x(1:end-3),cfd.wav)
        hold on
        plot(x,y)
        scatter(x(cfd.zci+cfd.min.idx),cfd.wav(cfd.zci+cfd.min.idx));
        xlim([x(sig.startpoint.idx)-1, x(sig.startpoint.idx)+20]);
        title('CFD result');
        ylabel('Voltage, V');
        xlabel('Time, ns');
        grid on;
        hold off;
    end
    cfd.wav=[];  % dont save waveform to
    
    %% calculate sigmoid fit to the leading edge
    sigmoid.start = sig.startpoint.idx - 5;  % 5 samples in original code (why?)
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
    
    % cut the leading edge
    xfit = x(sigmoid.start:sigmoid.end);
    yfit = y(sigmoid.start:sigmoid.end);
    
    % fill initial parameters
    p0(1) = sig.max.y;
    p0(2) = 0.5*(x(sigmoid.end)+x(sigmoid.start));
    p0(3) = 5.0/(x(sigmoid.end)-x(sigmoid.start));
    p0(4) = 2.0;   % changed from 1.0 to 2.0
    p0(5) = 0.0;
    sigmoid.p0 = p0; % store p0
  
    % call unconstrained optimization (1e-3 for TolFun and Tolx seems like
    % good compromise) - note that fit is done from 1:end-5!!!
    optimopts = optimset('Display','off','TolFun',1e-3,'TolX', 1e-3);
    [p,fval,exitflag,outputs] = fminsearch(@(p)sigmoid_fit(xfit(1:end-5),yfit(1:end-5),p),p0,optimopts);
    sigmoid.p = p;
    sigmoid.fval = fval;
    sigmoid.exitflag = exitflag;
    sigmoid.niter = outputs.funcCount;
    
    % check fit if ploting enabled
    if(options.en_plot)
        
        subplot(2,2,4);
        scatter(xfit,yfit);
        hold on;
        plot(xfit,ferdirac5(xfit,p));
        title('Sigmoid fit result');
        ylabel('Voltage, V');
        xlabel('Time, ns');
        grid on
        hold off;
        drawnow;
        p
        sigmoid.niter
        pause(1);
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
    end    %% wrap outputs
    out.cfd = cfd;
    out.sig = sig;
    out.sigmoid = sigmoid;

return
end    %% AUX functions below

%% linear interpolation function
function x=lin_interp(xq,yq,y)
    b = (yq(2)-yq(1))/(xq(2)-xq(1));
    a = yq(2) - b*xq(2);
    x = (y-a)/b;
end

%% sigmoid fittness function
function out=sigmoid_fit(xq,yq,p)
    ycalc = ferdirac5(xq,p);
    out = (ycalc-yq).'*(ycalc-yq); % calculate  chi-square fittnes
end

%% 5 parameter fermi dirac function
function f=ferdirac5(x,p)
    f=p(1)./power((1+exp(-(x-p(2))*p(3))),p(4))+p(5);
end



