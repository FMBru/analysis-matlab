clear all
close all

load Run019.mat

tbl_idx = 1;
table(tbl_idx) = str2num(run.id);
tbl_idx = tbl_idx+1;

k=1;
for i=1:length(MCP_data)
    %exitflag_MCP(i) = MCP_data(i).sigmoid.exitflag;
    %exitflag_MM(i) = MM_data(i).sigmoid.exitflag;
    %fval_MM(i) = MM_data(i).sigmoid.fval;
    %fval_MCP(i) = MCP_data(i).sigmoid.fval;
    time_MCP(i) = MCP_data(i).sigmoid.timepoint;
    time_MM(i) = MM_data(i).sigmoid.timepoint;
    blavg_MM(i) = MM_data(i).sig.blavg;
    blavg_MCP(i) = MCP_data(i).sig.blavg;
    blrms_MM(i) = MM_data(i).sig.blrms;
    blrms_MCP(i) = MCP_data(i).sig.blrms;
    ymax_MCP(i)=MCP_data(i).sig.max.y;
    ymax_MM(i)=MM_data(i).sig.max.y;
    %if(exitflag_MCP(i)+exitflag_MM(i)==0)
end

%% do some initial resolution measurement
% make simple cut with respect to median time and minimum amplitude magic
% numbers :-)
time_avg_raw = median(time_diff_sigmoid);  % assume mu from the median vaule
time_min = time_avg_raw-0.3;     % predicted resolution 100ps cut 3 sigma
time_max = time_avg_raw+0.3;     % left and 3 sigma right from median

% find events within the cut and cut with respect to amplirude ugly :-)
%idx_cut = time_diff_sigmoid>time_min & time_diff_sigmoid<time_max & MCP_maxy>0.1*max(MCP_maxy) & MM_maxy>0.05*max(MM_maxy);
idx_cut = time_diff_sigmoid>time_min & time_diff_sigmoid<time_max & MCP_maxy>0.1*max(MCP_maxy) & MM_maxy>0.01*max(MM_maxy)& MM_maxy<0.99*max(MM_maxy)& MCP_maxy<0.95*max(MCP_maxy);

time_diff_cut=time_diff_sigmoid(idx_cut);

%% plot time difference histogram
figure;
h=histogram(time_diff_cut,100);
hold on
fit_data = [];
fit_data(1,:) = h.BinEdges(1:end-1)+h.BinWidth/2;
fit_data(2,:) = h.Values;
% fit_data(3,:) = yerr;
p0=[];
p0(1) = sum(fit_data(2,:))*h.BinWidth;;   % normalization factor
p0(2) = std(h.Data);           % sigma
p0(3) = mean(fit_data(1,:));   % mean
cmd='min; ret';
[p, err, chi] = fminuit('gauss1_minuit',p0,fit_data(:,1:end),'-b','-c',cmd);
plot(fit_data(1,:),gauss1_minuit(p,fit_data(1,:)),'Linewidth',2);

grid on
ylabel('events')
xlabel('Time difference, ns');
title_str = sprintf('%s: \\mu = %4.4f ns \\sigma = %4.4f ns',run.name, p(3), p(2));
title(title_str)

table(tbl_idx) = p(3); 
tbl_idx = tbl_idx+1;
table(tbl_idx) = p(2); 
tbl_idx = tbl_idx+1;


%% plot dual gauss histogram
figure
%rng(abs(pd.mu));
h = histogram(time_diff_cut,100);
hold on
fit_data = [];
fit_data(1,:) = h.BinEdges(1:end-1)+h.BinWidth/2;
fit_data(2,:) = h.Values;
% fit_data(3,:) = yerr;
p0=[];
p0(1) = 0.2*sum(fit_data(2,:))*h.BinWidth/std(h.Data)/sqrt(2*pi);   % normalization of gauss
p0(2) = std(h.Data);           % sigma1
p0(3) = mean(fit_data(1,:));   % mean1
p0(4) = 0.8*sum(fit_data(2,:))*h.BinWidth/std(h.Data)/sqrt(2*pi);
p0(5) = std(h.Data);           % sigma2
p0(6) = mean(fit_data(1,:));   % mean1
cmd='min; ret';
[p, err, chi] = fminuit('gauss2_minuit',p0,fit_data(:,1:end),'-b','-c',cmd);
chi
plot(fit_data(1,:),gauss2_minuit(p,fit_data(1,:)),'Linewidth',2);
plot(fit_data(1,:),gauss1_minuit(p(1:3),fit_data(1,:)),'Linewidth',2);
plot(fit_data(1,:),gauss1_minuit(p(4:6),fit_data(1,:)),'Linewidth',2);
g2.p=p(1)/(p(1)+p(4));
g2.mu=[p(3) p(6)];
g2.sigma = [p(2) p(5)];
g2.sigma_all = sqrt(g2.p*g2.sigma(1)^2+(1-g2.p)*g2.sigma(2)^2+g2.p*(1-g2.p)*(g2.mu(1)-g2.mu(2))^2);
g2.mu_all = g2.p*g2.mu(1)+(1-g2.p)*g2.mu(2);
g2.chi = chi;
g2.err = err;

xlabel('Time difference, ns');
ylabel('events');
legend('RAW hist','Gauss combined','Gauss 1', 'Gauss 2');
grid
title_str = sprintf('2Gauss %s: \\mu = %4.4f ns \\sigma = %4.4f ns',run.name, g2.mu_all, g2.sigma_all);
title(title_str)

table(tbl_idx) = g2.mu_all; 
tbl_idx = tbl_idx+1;
table(tbl_idx) = g2.sigma_all; 
tbl_idx = tbl_idx+1;


%% plot peak amplitude histogram cut included
figure
histogram(MCP_maxy(idx_cut),100)
hold on
h=histogram(MM_maxy(idx_cut),100);
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
title_str = sprintf('%s - e-peak amplitude \\mu = %4.4f V U_{max} = %4.4f V',run.name, e_peak_amp.mean,xbins(e_peak_amp.max_idx));
title(title_str)

table(tbl_idx) = e_peak_amp.mean; 
tbl_idx = tbl_idx+1;
table(tbl_idx) = xbins(e_peak_amp.max_idx); 
tbl_idx = tbl_idx+1;



%% plot electron charge
for i=1:length(MM_data)
MM_e_charge(i)=MM_data(i).sig.charge.e_peak;
end
figure;
plot(MM_maxy, MM_e_charge,'.');
grid on
xlabel('Signal amplitude, V')
ylabel('Electron charge, pC');
title_str = sprintf('%s - Electron charge plot',run.name);
title(title_str)

%% plot electron lead charge histogram
MM_e_lead_charge = zeros(length(MM_data),1);
for i=1:length(MM_data)
MM_e_lead_charge(i)=MM_data(i).sig.charge.lead_edge;
end
figure;
h=histogram(MM_e_lead_charge(idx_cut),100);
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

xlabel('Electron charge, pC');
ylabel('events');
grid on
legend('Hist','Polya fit');
title_str = sprintf('%s - e-lead charge \\mu = %4.4f pC Q_{max} = %4.4f pC',run.name, e_peak_amp.mean,xbins(e_peak_amp.max_idx));
title(title_str)

table(tbl_idx) = e_peak_amp.mean; 
tbl_idx = tbl_idx+1;
table(tbl_idx) = xbins(e_peak_amp.max_idx); 
tbl_idx = tbl_idx+1;

%% plot electron peak charge histogram
MM_e_peak_charge = zeros(length(MM_data),1);
for i=1:length(MM_data)
MM_e_peak_charge(i)=MM_data(i).sig.charge.e_peak;
end
figure;
h=histogram(MM_e_peak_charge(idx_cut),100);
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

xlabel('Electron charge, pC');
ylabel('events');
grid on
legend('Hist','Polya fit');
title_str = sprintf('%s - e-peak charge \\mu = %4.4f pC  Q_{max} = %4.4f pC',run.name, e_peak_amp.mean,xbins(e_peak_amp.max_idx));
title(title_str)

table(tbl_idx) = e_peak_amp.mean; 
tbl_idx = tbl_idx+1;
table(tbl_idx) = xbins(e_peak_amp.max_idx); 
tbl_idx = tbl_idx+1;

%% plot 10-90% leadedge
MM_t_rise = zeros(length(MM_data),1);
for i=1:length(MM_t_rise)
    MM_t_rise(i)=calc_sig_rise_time(MM_data(i).sigmoid.p,MM_data(i).sig.startpoint.x-1,MM_data(i).sig.max.y);
end
figure;
h=histogram(MM_t_rise(idx_cut),100);
hold on
fit_data = [];
fit_data(1,:) = h.BinEdges(1:end-1)+h.BinWidth/2;
fit_data(2,:) = h.Values;
% fit_data(3,:) = yerr;
p0=[];
p0(1) = sum(fit_data(2,:))*h.BinWidth;   % normalization factor
p0(2) = std(h.Data);           % sigma
p0(3) = mean(fit_data(1,:));   % mean
cmd='min; ret';
[p, err, chi] = fminuit('gauss1_minuit',p0,fit_data(:,1:end),'-b','-c',cmd);
plot(fit_data(1,:),gauss1_minuit(p,fit_data(1,:)),'Linewidth',2);

grid on
ylabel('events')
xlabel('Risetime, ns');
title_str = sprintf('%s Rise time 10-90%% \\mu = %4.4f ns \\sigma = %4.4f ns',run.name, p(3), p(2));
title(title_str)

table(tbl_idx) =  p(3); 
tbl_idx = tbl_idx+1;
table(tbl_idx) =  p(2); 
tbl_idx = tbl_idx+1;

%% plot endpoint for integration of electron peak
for i=1:length(MM_data)
    MM_e_peak_end(i) = 1e9*100e-12*(MM_data(i).sig.e_peak_end.idx-MM_data(i).sig.max.idx);
end
figure;
plot(MM_maxy, MM_e_peak_end,'.');
grid on
xlabel('Signal amplitude, V')
ylabel('Time, ns');
title_str = sprintf('%s - Electron charge end time',run.name);
title(title_str)
% plot histogram
figure
h=histogram(MM_e_peak_end,177);
grid on
ylabel('events')
xlabel('Time, ns');
title_str = sprintf('%s - Electron charge end time',run.name);
title(title_str)

Y = fft(h.Values);
L = length(h.Values);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = 1000/h.BinWidth*(0:(L/2))/L;
[ares fres_idx] = max(P1(5:end));
fres = f(4+fres_idx);

figure
plot(f,P1)
title('Spectrum of oscillations of end time')
xlabel('f, MHz')
ylabel('Amplitude')
grid on

table(tbl_idx) =  fres; 
tbl_idx = tbl_idx+1




