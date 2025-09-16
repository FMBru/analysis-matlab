close all

%% histogram and polya fit
format long; % print 15 decimals intead of 4delat
%clear ch_tr ch_mmclose

%geometric cut
%radius = 2; % 2 mm radius

spacial_cut = 1;
%% Parameters for changing

shouldOutputAmplitudesTxtFile = true;

useNaiveMeanPos = true; %if median pos extract fails, use simple mean of tracker values instead

% start / end voltages for fitting
% if trackerExist == 1
%     min_v = 30e-3;
%     max_v = 0.22;
%     % bins
%     binOffset = 0;
% else
    min_v = 0;
    max_v =4e-3;
    % bins
    binOffset = 0;

%end

% polya starting parameters for fitting
%x0 = [1 1 0.2];  %parameters 3LEDs test ortec
%x0 = [2 2 0.5];

figureWidth=800;
figureHeight=500;


num_bins = 200;

% if trackerExist == 1
% 
%     %sampling area for 2D maps
%     %% calculate
%     area.step = 0.25;   % set grid resolution in mm
%     area.size = 8;      % set half-size of the observed square area, mm - for
%     %pads and small MCP
%     %area.size = 20;      % set half-size of the observed square area, mm - for large MCP
%     area.radius = 1;    % set radius of averaging circular window cut, mm
% end


%% Global cut
% if trackerExist == 1
%     glbl_cut = mm_max_y>0.00 & mm_max_y<0.95*max(mm_max_y) & trackerX' ~= 0 & trackerY' ~= 0; %remove saturated datapoints
% else
    glbl_cut = mm_max_y>0.00 & mm_max_y<0.95*max(mm_max_y)
    
%end



figure(1)
x_fig=10;
y_fig=10;
set(gcf,'position',[x_fig,y_fig,figureWidth,figureHeight]);
h=histogram(mm_max_y(glbl_cut),num_bins)
set(gca, 'YScale', 'log') %plot on log scale
movegui(gcf,'south');
hold on
xlabel('Signal amplitude, V');
ylabel('Events');
ylim([0.5 max(h.Values)+1000]);
ax = gca;
ax.FontSize = 20;
if trackerExist == 1
    % xlim([-0.005 0.3]);
else
    %xlim([-0.005 0.1]);
end
grid on
if trackerExist == 1
    title_str = sprintf('PICOSEC beam test - Run %s - Max e-peak amplitude',run.id);
else
    title_str = sprintf('PICOSEC LED test - Run %s - Max e-peak amplitude',run.id);
end
title(title_str)

%could delete following lines
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

% array of starting voltages is the first n_bins starting with
min_arr = h.BinEdges(min_cut_idx:min_cut_idx+binOffset);

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
colors = distinguishable_colors(binOffset+1);
c = 1;
%P = ["N","theta","nBar"];
disp("ciao")
for i=1:length(min_arr)
    cut_idx = min_cut_idx+i-1;
    % define polya
    %fitfun = fittype( @(N,theta,nBar,x) (N./nBar).*((theta+1).^(theta+1)).*((x./nBar).^theta).*exp(-(theta+1).*x./nBar)./gamma(theta+1));
    
    %define landau
    %fitfun = fittype(@(x) exp(-(x+exp(-x)/2))/sqrt(2*pi));
    
    %options = fitoptions(fitfun);
    %options.Upper = [100 5 1];
    %options.Lower = [0 0 0];
    %options.StartPoint = x0;
   
    %try
        %fitted_curve = fit((h.BinEdges(cut_idx:cutUp_idx)+h.BinWidth/2)',h.Values(cut_idx:cutUp_idx)','gauss1')
        fitted_curve = fit((h.BinEdges(1:length(h.Values))+h.BinWidth/2)',h.Values','gauss2')
   
    % catch exception
    %     disp("ciao3")
    %     exceptionCounter = exceptionCounter+1;
    %     mean_arr(i) = 0;
    %     chi2_arr(i) = 0;
    %     mean_err_arr(i,:) = [0;0];
    %     continue
    %     exceptionCounter
    % end
    
 %   P = [P;fitted_curve.N,fitted_curve.theta,fitted_curve.nBar];
    
    % plot the fitted curve
%     if (fitted_curve.theta > 0)
%         hold on
%          plot(h.BinEdges(cut_idx:cutUp_idx),fitted_curve(h.BinEdges(cut_idx:cutUp_idx)),'LineWidth',3,'Color',[colors(c,:)])
%          xlim([0 4e-3] );


          %%save fit params
 gaus_mean = [];
 gaus_sigma = [];

 gaus_mean = [gaus_mean; fitted_curve.b1];
 gaus_sigma = [gaus_sigma; fitted_curve.c1];

 %%plot the fitted curve
 plot(h.BinEdges(cut_idx:cutUp_idx),fitted_curve(h.BinEdges(cut_idx:cutUp_idx)),'LineWidth',3,'Color',[colors(c,:)])
 xlim([0 4e-3] );
 
 dlmwrite(['noiseGausseFit.txt'],[fitted_curve.b1; fitted_curve.c1] )

%%extract fit errors
cfit_curve = cfit(fitted_curve);
param_errors = confint(cfit_curve);
err_b1 = param_errors(2, 2)-fitted_curve.b1;
err_c1 = param_errors(2, 3)-fitted_curve.c1;


%%save fit
%store_folder = ['\\eosproject-smb\eos\project\p\picosec\testbeam\2023_August_h4\Results\SPE\Gauss+Polya\Run' run.id '\Noise\'];
%store_folder = ['\\eosproject-smb\eos\project\p\picosec\lab\Chiara\resistive\results\' run.id '\'];
% store_folder = ['\\eosproject-smb\eos\project\p\picosec\lab\Chiara\antonija\results\' run.id '\'];
% store_folder = ['\\eosproject-smb\eos\project\p\picosec\lab\Chiara\amplitudeVsTime\C103\results\' run.id '\'];
%mkdir(store_folder);

str1 = sprintf('Mean = %0.5f +/- %0.5f\n Sigma = %0.5f +/- %0.5f',fitted_curve.b1,err_b1,fitted_curve.c1,err_c1);
annotation('textbox',[.50 .6 .35 .3], 'String',str1,'FitBoxToText','on','FontSize',18);
saveas(figure(1),[store_folder 'RUN' run.id ' - Channel' channel.id '- E-PeakAmplitude-Noise.png'])
% saveas(figure(1),[store_folder ' Channel' channel.id '- E-PeakAmplitude.png'])

%         % used color, increment color counter
%     
%         % calculate chi squared and degrees of freedom
%         ch2 = sum(((h.Values(cut_idx:cutUp_idx)-fitted_curve(h.BinEdges(cut_idx:cutUp_idx))').^2)./fitted_curve(h.BinEdges(cut_idx:cutUp_idx))');
%         dof = size(h.Values(cut_idx:cutUp_idx),2)-3;
%         nch2 = ch2/dof
%         np = 1-chi2cdf(ch2,dof) % P(\chi^2>ch2)
%     
%         mean_arr(i) = fitted_curve.nBar;
%         chi2_arr(i) = nch2;
%         mean_err = confint(cfit(fitted_curve),0.68);
%     
%         c = c+1;
%     else
%         exceptionCounter = exceptionCounter+1;
%         mean_arr(i) = 0;
%         chi2_arr(i) = 0;
%     end
end
%ProcessRawFiles_PEAnalysis_noise_plus_sig



% close all
% 
% %% histogram and polya fit
% format long; % print 15 decimals intead of 4delat
% %clear ch_tr ch_mmclose
% 
% %geometric cut
% radius = 2; % 2 mm radius
% 
% spacial_cut = 1;
% %% Parameters for changing
% 
% shouldOutputAmplitudesTxtFile = true;
% 
% useNaiveMeanPos = true; %if median pos extract fails, use simple mean of tracker values instead
% 
% % start / end voltages for fitting
% if trackerExist == 1
%     min_v = 30e-3;
%     max_v = 0.22;
%     % bins
%     binOffset = 0;
% else
%     min_v = 0;
%     max_v =8e-3;
%     % bins
%     binOffset = 0;
% 
% end
% 
% % polya starting parameters for fitting
% %x0 = [1 1 0.2];  %parameters 3LEDs test ortec
% x0 = [2 2 0.5];
% 
% figureWidth=800;
% figureHeight=500;
% 
% 
% num_bins = 500;
% 
% if trackerExist == 1
% 
% %sampling area for 2D maps
% %% calculate
% area.step = 0.25;   % set grid resolution in mm
% area.size = 8;      % set half-size of the observed square area, mm - for
% %pads and small MCP
% %area.size = 20;      % set half-size of the observed square area, mm - for large MCP
% area.radius = 1;    % set radius of averaging circular window cut, mm
% end
% 
% 
% %% Global cut
% if trackerExist == 1
%     glbl_cut = mm_max_y>0.00 & mm_max_y<0.95*max(mm_max_y) & trackerX' ~= 0 & trackerY' ~= 0; %remove saturated datapoints
% else
%     glbl_cut = mm_max_y>0.00 & mm_max_y<0.95*max(mm_max_y)
%     
% end
% 
% figure(1)
% x_fig=10;
% y_fig=10;
% set(gcf,'position',[x_fig,y_fig,figureWidth,figureHeight]);
% h=histogram(mm_max_y(glbl_cut),num_bins)
% set(gca, 'YScale', 'log') %plot on log scale
% movegui(gcf,'south');
% hold on
% xlabel('Signal amplitude, V');
% ylabel('Events');
% ylim([0.5 max(h.Values)+1000]);
% ax = gca;
% ax.FontSize = 20;
% grid on
% if trackerExist == 1
%     title_str = sprintf('PICOSEC beam test - Run %s - Max e-peak amplitude',run.id);
% else
%     title_str = sprintf('PICOSEC LED test - Run %s - Max e-peak amplitude',run.id);
% end
% title(title_str)
% 
% %could delete following lines
% % set the upper and lower bounds for the lower cutoff of the polya fit that
% % we would like to try
% min_lower = 3e-3;
% %min_upper = 3.5e-3;
% 
% % If min and max for fit are not within dataset range pick edges of dataset
% if (min(h.BinEdges) > min_lower)
%     min_lower = min(h.BinEdges);
% end
% if (max(h.BinEdges(1:num_bins)) < max_v)
%     max_v = max(h.BinEdges(1:num_bins)) % cut the last bin edge because it's on the right side of the bin to make array lengths match
% end
% %max_v = max(h.BinEdges(1:num_bins));
% 
% 
% % min cut value should be the closest value to min_lower that is greater than
% % min_v
% [min_cut_val,min_cut_idx] = min(abs(h.BinEdges-min_v));
% if (min_cut_val >= min_v)
% else
%     min_cut_idx=min_cut_idx+1;
% end
% [diff,cutUp_idx] = min(abs(h.BinEdges-max_v));
% cutUp_val = h.BinEdges(cutUp_idx);
% if (cutUp_val <= max_v)
% else
%     cutUp_idx=cutUp_idx-1;
% end
% 
% % array of starting voltages is the first n_bins starting with
% min_arr = h.BinEdges(min_cut_idx:min_cut_idx+binOffset);
% 
% % initialize arrays
% mean_arr = [];
% chi2_arr = [];
% 
% % take upper cut as very end of array
% %cutUp_idx = length(h.Values);
% 
% % count exceptions
% exceptionCounter = 0;
% 
% % initialize error array
% mean_err_arr = [];
% 
% % colors for plotting
% colors = distinguishable_colors(binOffset+1);
% c = 1;
% 
% % for i=1:length(min_arr)
%      cut_idx = min_cut_idx+i-1;
%     
%     try
%         %fitted_curve = fit((h.BinEdges(cut_idx:cutUp_idx)+h.BinWidth/2)',h.Values(cut_idx:cutUp_idx)','gauss1')
%         fitted_curve = fit((h.BinEdges(1:length(h.Values))+h.BinWidth/2)',h.Values','gauss2')
%     catch exception
%         exceptionCounter = exceptionCounter+1
%         mean_arr(i) = 0;
%         chi2_arr(i) = 0;
%         mean_err_arr(i,:) = [0;0];
%        % continue
%     end
%  
%  %%save fit params
%  gaus_mean = [];
%  gaus_sigma = [];
% 
%  gaus_mean = [gaus_mean; fitted_curve.b1];
%  gaus_sigma = [gaus_sigma; fitted_curve.c1];
% 
%  %%plot the fitted curve
%  plot(h.BinEdges(cut_idx:cutUp_idx),fitted_curve(h.BinEdges(cut_idx:cutUp_idx)),'LineWidth',3,'Color',[colors(c,:)])
%  xlim([0 4e-3] );
%  
%  dlmwrite(['noiseGausseFit.txt'],[fitted_curve.b1; fitted_curve.c1] )
% 
% %%extract fit errors
% cfit_curve = cfit(fitted_curve);
% param_errors = confint(cfit_curve);
% err_b1 = param_errors(2, 2);
% err_c1 = param_errors(2, 3);
% 
% 
% %%save fit
% store_folder = ['\\eosproject-smb\eos\project\p\picosec\testbeam\2023_July_h4\Results\SPE\Gauss+Polya\Run' run.id '\Noise\'];
% mkdir(store_folder);
% 
% str1 = sprintf('Mean = %0.5f +/- %0.5f\n Sigma = %0.5f +/- %0.5f',fitted_curve.b1,err_b1,fitted_curve.c1,err_c1);
% annotation('textbox',[.50 .6 .35 .3], 'String',str1,'FitBoxToText','on','FontSize',18);
% saveas(figure(1),[store_folder 'RUN' run.id ' - Channel' channel.id '- E-PeakAmplitude.png'])
% 
% %         % used color, increment color counter
% %     
% %         % calculate chi squared and degrees of freedom
% %         ch2 = sum(((h.Values(cut_idx:cutUp_idx)-fitted_curve(h.BinEdges(cut_idx:cutUp_idx))').^2)./fitted_curve(h.BinEdges(cut_idx:cutUp_idx))');
% %         dof = size(h.Values(cut_idx:cutUp_idx),2)-3;
% %         nch2 = ch2/dof
% %         np = 1-chi2cdf(ch2,dof) % P(\chi^2>ch2)
% %     
% %         mean_arr(i) = fitted_curve.nBar;
% %         chi2_arr(i) = nch2;
% %         mean_err = confint(cfit(fitted_curve),0.68);
% %     
% %         c = c+1;
% %     else
% %         exceptionCounter = exceptionCounter+1;
% %         mean_arr(i) = 0;
% %         chi2_arr(i) = 0;
% %     end
% 
% 
% % end
% 
% %ProcessRawFiles_PEAnalysis_noise_plus_sig
