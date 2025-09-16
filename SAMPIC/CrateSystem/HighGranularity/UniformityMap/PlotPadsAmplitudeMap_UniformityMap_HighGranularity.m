%take padResults vector and visualise

%config for fitting
addpath 'C:\Users\GDD\Documents\MATLAB\Picosec\SPEAnalysisBeam'

%% Parameters for chaning

% cutoffs\
min_v = 30e-3;
max_v = 0.3;

% bins
n_bins = 0;



% polya starting parameters for fitting

clear MM_maxyDUT

        padIDsArray = [1:37];
          channelToAnalyse = [1:37];
        for i=1:length(padIDsArray)
            channelToAnalyse (i) = getChannelForPadNumberHighGranularity(padIDsArray(i));
        end


%for chPos = 1:length(channelToAnalyse)
for chPos = 1:length(channelToAnalyse) %loop through hannels, starrt from other channel if fit fails

    clearvars -except padIDsArray padCounts2DArray padMeanArray padMeanFitArray padCountsArray run store_folder dutChannelArray MM_maxy chPos padCHArray padMaxYArray padMaxYFitArray padCountsArray min_v max_v n_bins channelToAnalyse
        x0 = [1 1 0.2];  %parameters 3LEDs test ortec

    if chPos==101
        
    else
        ch = channelToAnalyse(chPos);
        
        channelMask = dutChannelArray==ch;
        MM_maxyDUT = MM_maxy(channelMask);
        meanAmp = mean(MM_maxyDUT);
        counts = length(MM_maxyDUT);
        
        if counts>10
        
        padCHArray(chPos) = ch;
        padMaxYArray(chPos) = 1000*meanAmp;
        padCountsArray(chPos) = counts;
        end
        %fitting specturm
        
        mm_max_y = MM_maxyDUT;
        
        glbl_cut = mm_max_y>0.00 & mm_max_y<0.95*max(mm_max_y)
        
        num_bins = 250;
        
        figure
        x_fig=10;
        y_fig=10;
        width=1500;
        height=700;
        set(gcf,'position',[x_fig,y_fig,width,height]);
        h=histogram(mm_max_y(glbl_cut),num_bins)
        set(gca, 'YScale', 'log') %plot on log scale
        hold on
        xlabel('Signal amplitude, V');
        ylabel('Events');
        ylim([0.5 max(h.Values)+1000]);
        ax = gca;
        ax.FontSize = 20;
        xlim([-0.005 0.5]);
        grid on
        title_str = sprintf('PICOSEC uniformity test - Run %s - Max e-peak amplitude',run.id);
        title(title_str)
        
        % set the upper and lower bounds for the lower cutoff of the polya fit that
        % we would like to try
        min_lower = 10e-3;
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
        minCutIndex = min_cut_idx+n_bins;
        if minCutIndex>length(h.BinEdges)
            minCutIndex = length(h.BinEdges);
        end
        % array of starting voltages is the first n_bins starting with
        min_arr = h.BinEdges(min_cut_idx:minCutIndex);
        
        % initialize arrays
        mean_arr = [];
        chi2_arr = [];
        
        % take upper cut as very end of array
        %cutUp_idx = length(h.Values);
        
        % count exceptions
        exceptionCounter = 0;
        
        % initialize error array
        mean_err_arr = [];
        
        % colors for plotting
        colors = distinguishable_colors(n_bins+1);
        c = 1;
        P = ["N","theta","nBar"];
        
        %for i=1:length(min_arr)
        for i=1:1
            cut_idx = min_cut_idx+i-1;
            % define polya
            fitfun = fittype( @(N,theta,nBar,x) (N./nBar).*((theta+1).^(theta+1)).*((x./nBar).^theta).*exp(-(theta+1).*x./nBar)./gamma(theta+1));
            fitfun
            options = fitoptions(fitfun);
            options.Upper = [max(h.Values) 5 1];
            options.Lower = [0 0 0];
            options.StartPoint = x0;
            
            try
                [fitted_curve,gof] = fit(h.BinEdges(cut_idx:cutUp_idx)',h.Values(cut_idx:cutUp_idx)',fitfun,options)
            catch exception
                exceptionCounter = exceptionCounter+1;
                mean_arr(i) = 0;
                chi2_arr(i) = 0;
                mean_err_arr(i,:) = [0;0];
                continue
            end
            
            P = [P;fitted_curve.N,fitted_curve.theta,fitted_curve.nBar];
            
            
            % plot the fitted curve
            if (fitted_curve.theta > 0)
                hold on
                plot(h.BinEdges(cut_idx:cutUp_idx),fitted_curve(h.BinEdges(cut_idx:cutUp_idx)),'LineWidth',3,'Color',[colors(c,:)])
                % used color, increment color counter
                
                % calculate chi squared and degrees of freedom
                ch2 = sum(((h.Values(cut_idx:cutUp_idx)-fitted_curve(h.BinEdges(cut_idx:cutUp_idx))').^2)./fitted_curve(h.BinEdges(cut_idx:cutUp_idx))');
                dof = size(h.Values(cut_idx:cutUp_idx),2)-3;
                nch2 = ch2/dof
                np = 1-chi2cdf(ch2,dof) % P(\chi^2>ch2)
                
                mean_arr(i) = fitted_curve.nBar;
                chi2_arr(i) = nch2;
                mean_err = confint(cfit(fitted_curve),0.68);
                mean_err_arr(i,:) = mean_err(:,3);
                
                c = c+1;
            else
                exceptionCounter = exceptionCounter+1;
                mean_arr(i) = 0;
                chi2_arr(i) = 0;
                mean_err_arr(i,:) = [0;0];
            end
            
        end
        
        % print out file with params
        % filename = [run.id, 'optParams.csv'];
        % filedir = ['C:\Users\GDD\Documents\MATLAB\Picosec\SPEAnalysisPolyaMeanFinalVersion2\',filename];
        % writematrix(P,filename);
        % fileattrib(filedir,'+w','a');
        % annotate and label plot
        %text(min(h.BinEdges)+0.02*max(h.BinEdges),max(h.Values)+2,'\chi^2/dof = '+string(round(ch2)) + '/' + string(size(h.Values(cut_idx:cutUp_idx),2)) + ' = ' + string(nch2),'FontSize',14)
        %text();
        
        %store_folder = ['\\eosproject-smb\eos\project\p\picosec\testbeam\2022_July_h4\Results\TestMarta\'];
        %mkdir(store_folder);
        
        
        %% calculate the final mean from all means
        format shortE
        % mean
        % remove 0s which correspond to fits that threw exceptions, do not want 0s
        % biasing the mean
        filter = mean_arr ~= 0;
        mean_arr = mean_arr(filter);
        chi2_arr = chi2_arr(filter);
        mean_err_arr_temp = mean_err_arr(filter,:);
        mean_err_arr = mean_err_arr_temp;
        min_arr = min_arr(filter);
        meanSum = sum(mean_arr)/length(mean_arr);
        
        %error propagation
        delta_arr(1) = (1/2)*(sum((mean_err_arr(:,1)'-mean_arr).^2))^(1/2); % lower
        delta_arr(2) = (1/2)*(sum((mean_err_arr(:,2)'-mean_arr).^2))^(1/2); % upper
        
        % print values
        %str = sprintf('Mean: %f \n Upper bound: %f \n Lower bound: %f',mean,mean+delta_arr(2),mean-delta_arr(1));
        %display(str);
        
        str1 = sprintf('CH: %d - Mean = %0.5f +/- %0.5f',ch,meanSum,delta_arr(2));
        annotation('textbox',[.50 .6 .35 .3], 'String',str1,'FitBoxToText','on','FontSize',18);
        %    saveas(figure(1),[store_folder 'RUN' run.id ' - Channel' channel.id '- E-PeakAmplitude.png'])
        
        
                if counts>10

        padMaxYFitArray(chPos) = 1000*meanSum;
                end
        % hist of amp only
%         
%         figure;
%         hist(MM_maxyDUT,1000); hold on
%         grid on
%         xlabel('Signal amplitude, V')
%         ylabel('Counts');
%         title_str = sprintf('CH %d \n Amplitude ',ch);
%         title(title_str)
%         xlim([0 0.3]);
%         vline(meanAmp,'g','Mean');
%         vline(meanSum,'r','MeanFit');
%         saveas(gcf,[store_folder '\ch' num2str(ch) '-amp.png'])
%         
        chPos
        pause(1);
        
        [row,col] = find(padIDsArray==ch)
        padMeanArray(row,col) = 1000*meanAmp;
        padCounts2DArray(row,col) = counts;
        padMeanFitArray(row,col) = 1000*meanSum;

        
                padMaxYArray(chPos) = 1000*meanAmp;
        


        padCountsArray(chPos) = counts;

        
        close all
    end
end

VisualiseHighGranularity(channelsEnabled,chCountsArray,'Hitmap Multipad Picosec Relative',[store_folder '\Run' run.id '_HitmapRelative.png'],'Hits','%0.1f %',0,0);



    
figure (5);
plot(padCHArray,padCountsArray,'.');
title('Hits per channel');
xlabel('Channel');
ylabel('Hits');
ax = gca;
grid on;
saveas(figure (5),[store_folder '\RUN' run.id  '- Hits.png'])

figure (6);
plot(padCHArray,padMaxYArray,'.');
title('Mean Amp');
xlabel('Channel');
ylabel('Mean amp');
ax = gca;
grid on;
saveas(figure (6),[store_folder '\RUN' run.id  '- MaxAmp.png'])

figure (7);
plot(padCHArray,padMaxYFitArray,'.');
title('Mean amp from fit');
xlabel('Channel');
ylabel('Mean amp');
ax = gca;
grid on;
saveas(figure (7),[store_folder '\RUN' run.id  '- MaxAmpFromFit.png'])



