clear all
close all

load Run045.mat

%% extract from structurearrays
k=1;
for i=1:length(MCP_data)
    time_MCP(i) = MCP_data(i).sigmoid.timepoint;
    time_MM(i) = MM_data(i).sigmoid.timepoint;
    blavg_MM(i) = MM_data(i).sig.blavg;
    blavg_MCP(i) = MCP_data(i).sig.blavg;
    blrms_MM(i) = MM_data(i).sig.blrms;
    blrms_MCP(i) = MCP_data(i).sig.blrms;
    ymax_MCP(i)=MCP_data(i).sig.max.y;
    ymax_MM(i)=MM_data(i).sig.max.y;
    e_peak_MM(i)=MM_data(i).sig.charge.e_peak;
    %if(exitflag_MCP(i)+exitflag_MM(i)==0)
end

%% do some initial resolution measurement
% make simple cut with respect to median time and minimum amplitude magic
% numbers :-)
time_avg_raw = median(time_diff_sigmoid);  % assume mu from the median vaule
time_min = time_avg_raw-0.3;     % predicted resolution 100ps cut 3 sigma
time_max = time_avg_raw+0.3;     % left and 3 sigma right from median

% find events within the cut and cut with respect to amplirudes and
% existance of the tracker data (GLOBAL CUT)
glbl_cut = time_diff_sigmoid>time_min & time_diff_sigmoid<time_max & MCP_maxy>0.01*max(MCP_maxy) & MM_maxy>0.01*max(MM_maxy)& MM_maxy<0.99*max(MM_maxy)& MCP_maxy<0.99*max(MCP_maxy);

% make cut on time difference
time_diff_cut = time_diff_sigmoid(glbl_cut);
mean(time_diff_cut)
std(time_diff_cut)

%% Time walk analysis (this needs to be improoved)
twalk.en = 0;                                   % enable timewalk correction
twalk.n_epk = 50;                               % number of e-peak bins
e_peak_srt = sort(e_peak_MM(glbl_cut));         % make temporary sort of e-peaks
% try to distrubute e_peak vector evenly over the amplitude range 
twalk.epeak_vec = e_peak_srt(1:round(length(e_peak_srt)/twalk.n_epk):end); % 50 bins
for i=1:length(twalk.epeak_vec)-1
    temp_cut = e_peak_MM > twalk.epeak_vec(i) & e_peak_MM < twalk.epeak_vec(i+1) & glbl_cut;
    twalk.mean_sat(i) = mean(time_diff_sigmoid(temp_cut));
    twalk.e_peak(i) =  mean(e_peak_MM(temp_cut));
    twalk.npts(i) = sum(temp_cut);
    twalk.rms(i) = std(time_diff_sigmoid(temp_cut));
    twalk.err(i) = std(time_diff_sigmoid(temp_cut))./sqrt(twalk.npts(i));
    twalk.e_peak_err_p(i) = twalk.e_peak(i)-(twalk.epeak_vec(i+1));
    twalk.e_peak_err_n(i) = -twalk.e_peak(i)+(twalk.epeak_vec(i));
end

% fit correction function using minuit
fit_data = [];
fit_data(1,:) = twalk.e_peak;
fit_data(2,:) = twalk.mean_sat;
fit_data(3,:) = twalk.err;
p0=[];
p0(1) = min(twalk.mean_sat);   
p0(2) = 1;        
p0(3) = 0.5;  
cmd='min; ret';
[p, err, chi] = fminuit('twalk_fn_minuit',p0,fit_data,'-b','-c',cmd);
twalk.p = p;

% plot time walk correction and fit
figure
errorbar(twalk.e_peak,twalk.mean_sat,twalk.err,twalk.err,twalk.e_peak_err_n,twalk.e_peak_err_p,'o');
hold on
plot(twalk.e_peak,twalk_fn_minuit(p,twalk.e_peak),'LineWidth',1.5);
xlabel('Electron peak charge, pC')
ylabel('SAT, ns')
grid

% plot resolution vs. e-charge
figure
errorbar(twalk.e_peak,twalk.rms*1000,[],[],twalk.e_peak_err_n,twalk.e_peak_err_p,'o');
hold on
xlabel('Electron peak charge, pC')
ylabel('Resolution, ps')
grid

% make time walk correction
if(twalk.en == 1)
    time_diff_sigmoid = time_diff_sigmoid - twalk_fn_minuit(twalk.p, e_peak_MM);
end

mean(time_diff_sigmoid);
std(time_diff_sigmoid);


%% plot peak amplitude histogram cut included
figure
histogram(MCP_maxy(glbl_cut),100)
hold on
h=histogram(MM_maxy(glbl_cut),100);
hold on
xbins = h.BinEdges(1:end-1)+h.BinWidth/2;
fit_data = [];
fit_data(1,:)=1:length(xbins);
fit_data(2,:) = h.Values;
% fit_data(3,:) = yerr;
p0=[];
p0(1) = sum(h.Values)*h.BinWidth;   % normalization factor
p0(2) = 1;         
p0(3) = 0.1;  
cmd='min; ret';
[p, err, chi] = fminuit('polya_minuit',p0,fit_data(:,1:end),'-b','-c',cmd);
polya_cnt_fit = polya_minuit(p,fit_data(1,:));
plot(xbins,polya_cnt_fit,'Linewidth',2);
e_peak_amp.mean = sum(polya_cnt_fit.*xbins)/sum(polya_cnt_fit);
[dummy, e_peak_amp.max_idx] = max(polya_cnt_fit);

xlabel('Signal amplitude, V')
ylabel('events');
grid on
legend('MCP','MM','MM fit');
title_str = sprintf('RUN %s - e-peak amplitude \\mu = %4.4f V U_{max} = %4.4f V',run.id, e_peak_amp.mean,xbins(e_peak_amp.max_idx));
title(title_str)

%% plot MM BG mean amplitude histogram cut included
figure
histogram(blavg_MM(glbl_cut),100);
blavg_MM_mean = mean(blavg_MM(glbl_cut));

xlabel('MM BG mean amplitude, V')
ylabel('events');
grid on
title_str = sprintf('RUN %s - MM BG mean amplitude = %4.4f mV',run.id, blavg_MM_mean*1000);
title(title_str);

%% plot MM BG mean RMS histogram cut included
figure
histogram(blrms_MM(glbl_cut),100);
blrms_MM_mean = mean(blrms_MM(glbl_cut));

xlabel('MM BG mean RMS, V')
ylabel('events');
grid on
title_str = sprintf('RUN %s - MM BG mean RMS = %4.4f mV',run.id, blrms_MM_mean*1000);
title(title_str);

%% plot dual gauss histogram
figure
%rng(abs(pd.mu));
h = histogram(time_diff_sigmoid(glbl_cut),100);
hold on
fit_data = [];
fit_data(1,:) = h.BinEdges(1:end-1) + h.BinWidth/2;
fit_data(2,:) = h.Values;
% fit_data(3,:) = yerr;

% try different combinations factors in order to find the min chi-square
comb = 0.05:0.05:0.95;
for i = 1:length(comb)
    p0=[];
    p0(1) = sum(fit_data(2,:))*h.BinWidth/std(h.Data)/sqrt(2*pi);   % total normalization factor
    p0(2) = std(time_diff_sigmoid(glbl_cut));                       % sigma1
    p0(3) = mean(time_diff_sigmoid(glbl_cut));                      % mean1
    p0(4) = comb(i);                                  % combination factotr
    p0(5) = std(time_diff_sigmoid(glbl_cut));                       % sigma2
    %p0(6) = mean(fit_data(1,:));                     % mean2 (not used) same as p3
    step = [1 2 3 4 5];
    l_b = [p0(1)*0.8 0.1*p0(2) p0(3)-0.1 0.5 0.1*p0(5)];
    u_b = [p0(1)*1.2 10*p0(2) p0(3)+0.1 0.999 10*p0(5)];
    StepBounds= [step; l_b; u_b]';
    MinuitCommands = 'min; ret;';
    chi_idx = fit_data(2,:) > 0;  % take into account only non-empty bins
    [p, err, chi2_min] = fminuit('gauss2_minuit', p0,fit_data(:,chi_idx), '-c', MinuitCommands,'-s', StepBounds);
    g2.chi2_vec(i) = chi2_min;
end
% find minimum chi-square
[chi2_min,idx_min_chi] = min(g2.chi2_vec);
p0(4) = comb(idx_min_chi);
[p, err, chi2_min] = fminuit('gauss2_minuit', p0,fit_data(:,chi_idx), '-c', MinuitCommands,'-s', StepBounds);

g2.chi2 = chi2_min;
g2.ndf = sum(chi_idx) - length(p);
plot(fit_data(1,:),gauss2_minuit(p,fit_data(1,:)),'Linewidth',2);
plot(fit_data(1,:),gauss1_minuit([p(1)*p(4) p(2) p(3)],fit_data(1,:)),'Linewidth',2);
plot(fit_data(1,:),gauss1_minuit([p(1)*(1-p(4)) p(5) p(3)],fit_data(1,:)),'Linewidth',2);
g2.p=p(4);
g2.mu = [p(3) p(3)];
g2.mu_err = [err(3) err(3)];
g2.sigma = [p(2) p(5)];
g2.sigma_err = [err(2) err(5)];
g2.sigma_all = sqrt(g2.p*g2.sigma(1)^2+(1-g2.p)*g2.sigma(2)^2+g2.p*(1-g2.p)*(g2.mu(1)-g2.mu(2))^2);
g2.mu_all = g2.p*g2.mu(1)+(1-g2.p)*g2.mu(2);

message = sprintf('  \\chi^2 / NDF = %2.1f / %d\n\n',g2.chi2,g2.ndf);
message = [message sprintf('  \\mu = %2.3f ns \\pm %2.3f ps\n',g2.mu_all,1000*g2.mu_err(1))];
message = [message sprintf('  \\sigma_1 = %2.1f ps \\pm %2.3f ps\n',1000*g2.sigma(1),1000*g2.sigma_err(1))];
message = [message sprintf('  \\sigma_2 = %2.1f ps \\pm %2.3f ps\n',1000*g2.sigma(2),1000*g2.sigma_err(2))];
message = [message sprintf('  \\sigma_{tot} = %2.1f ps \n',1000*g2.sigma_all)];
message = [message sprintf('  RMS_{tot} = %2.1f ps ',1000*std(time_diff_sigmoid(glbl_cut)))];

xlabel('Time difference, ns');
ylabel('events');
%xlim([-1 1]);
legend('RAW hist','Gauss combined','Gauss 1', 'Gauss 2');
grid
title_str = sprintf('2Gauss - %s',run.name);
title(title_str)
y_pos=get(gca,'ylim');
x_pos=get(gca,'xlim');
text(x_pos(1),0.75*y_pos(2),message)


%% save data to mat file
save(['Run' run.id '_plots.mat'], 'run', 'g2', 'twalk');

