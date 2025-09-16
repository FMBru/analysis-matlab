function tr = calc_sig_rise_time(p,start, max)

    % calculate timing from sigmoid
    sigmoid.time10=0;
    sigmoid.time90=0;
    sigmoid.sevals=0;
    sigmoid.timepoint10=start(1);
    
    % 20% val search algorithm
    % first go forward
    while(sigmoid.time10 < 0.1*max)
        sigmoid.sevals = sigmoid.sevals + 1;
        sigmoid.timepoint10 = sigmoid.timepoint10 + 0.01;          % 10ps stepping forward
        sigmoid.time10 = ferdirac5_minuit(p,sigmoid.timepoint10);
        if(sigmoid.sevals>1000)
            tr=0;  % maybe add fail here
            break;
        end   
    end
    % refine backwards
    while(sigmoid.time10 > 0.1*max)
        sigmoid.sevals = sigmoid.sevals + 1;
        sigmoid.timepoint10 = sigmoid.timepoint10 - 0.0001;        % 0.1ps stepping backwards
        sigmoid.time10 = ferdirac5_minuit(p,sigmoid.timepoint10);
        if(sigmoid.sevals>1000)
            tr=0;  % maybe add fail here
            break;
        end  
    end
    
    sigmoid.timepoint90 = sigmoid.timepoint10;
    
    % 90 % val search algorithm
    % first go forward
    while(sigmoid.time90 < 0.9*max)
        sigmoid.sevals = sigmoid.sevals + 1;
        sigmoid.timepoint90 = sigmoid.timepoint90 + 0.01;          % 10ps stepping forward
        sigmoid.time90 = ferdirac5_minuit(p,sigmoid.timepoint90);
        if(sigmoid.sevals>1000)
            tr=0;  % maybe add fail here
            break;
        end   
    end
    % refine backwards
    while(sigmoid.time90 > 0.9*max)
        sigmoid.sevals = sigmoid.sevals + 1;
        sigmoid.timepoint90 = sigmoid.timepoint90 - 0.0001;        % 0.1ps stepping backwards
        sigmoid.time90 = ferdirac5_minuit(p,sigmoid.timepoint90);
        if(sigmoid.sevals>1000)
            tr=0;  % maybe add fail here
            break;
        end   
    end
    tr = sigmoid.timepoint90-sigmoid.timepoint10;

end

