
glbl_cut = mm_max_y>0.00 & mm_max_y<0.99*max(mm_max_y);

figure(1)
h=histogram(mm_max_y(glbl_cut),100)
hold on

xbins = h.BinEdges(1:end-1)+h.BinWidth/2;

% remove bins below noise threshold
% indices = find(abs(xbins)<3e-3);
% xbins(indices) = [];
% h1.Values(indices) = [];


% remove all bins and values in bins less than noise threshold 3e-3, make a
% new array with them for Gaussian
idx = find(abs(xbins)<3E-3);
gauss_x(idx) = xbins(idx);
xbins(idx) = [];
hvals = h.Values;
gauss_vals(idx) = hvals(idx);
hvals(idx) = [];

% for fit with <3e-3 noise subtracted
fit_data = [];
fit_data(1,:)=1:length(xbins);
fit_data(2,:) = hvals;

% for curve to subtract noise function
fit_data_n = [];
fit_data_n(1,:)=1:length(gauss_x);
fit_data_n(2,:) = gauss_vals;



% for i = 1:length(xbins)
%     if (xbins(i)<3e-3)
%         xbins(i)=[];
%     end
% end



% for curve fit with <3e-3 noise subtracted
p0=[];
p0(1) = sum(hvals)*h.BinWidth;   % normalization factor
p0(2) = 1;         
p0(3) = 0.1;  
cmd='min; ret';

[p, err, chi] = fminuit('polya_minuit',p0,fit_data(:,1:end),'-b','-c',cmd);
polya_cnt_fit = polya_minuit(p,fit_data(1,:));
plot(xbins,polya_cnt_fit,'Linewidth',2);
e_peak_amp.mean = sum(polya_cnt_fit.*xbins)/sum(polya_cnt_fit);
[dummy, e_peak_amp.max_idx] = max(polya_cnt_fit);

% repeat for noise *******************************************
% but use the same bins as for the real signal so that we can subtract
% easily later
noise_glbl_cut = noise_mm_max_y>0.00 & noise_mm_max_y<0.99*max(noise_mm_max_y);

%figure
%noise_h =histogram(noise_mm_max_y(noise_glbl_cut),100)
%hold on
noise_arr = [];
for ( i = 1:length(h.BinEdges)-1 )
    bin_arr = noise_mm_max_y(noise_mm_max_y>=h.BinEdges(i) & noise_mm_max_y<h.BinEdges(i+1));
    noise_arr = [noise_arr;length(bin_arr)];
end

%noise_xbins = noise_h.BinEdges(1:end-1)+noise_h.BinWidth/2;


noise_fit_data = [];
noise_fit_data(1,:)=1:length(xbins);
noise_fit_data(2,:) = noise_arr;

% for curve fit
noise_p0=[];
noise_p0(1) = sum(noise_arr)*h.BinWidth;   % normalization factor
noise_p0(2) = 1;         
noise_p0(3) = 0.1;  
cmd='min; ret';

[noise_p, noise_err, noise_chi] = fminuit('polya_minuit',noise_p0,noise_fit_data(:,1:end),'-b','-c',cmd);
noise_polya_cnt_fit = polya_minuit(noise_p,noise_fit_data(1,:));
%plot(xbins,noise_polya_cnt_fit,'Linewidth',2);
e_peak_amp.mean = sum(noise_polya_cnt_fit.*xbins)/sum(noise_polya_cnt_fit);
[noise_dummy, noise_e_peak_amp.max_idx] = max(noise_polya_cnt_fit);

% ************************************************************
%f = fit(gauss_x.',gauss_vals.','gauss2'); %the .' just transposes the vector
%plot(f,gauss_x,gauss_vals,'or');

% subtract noise from signal

xlabel('Signal amplitude, V');
ylabel('events');
%xlim([-0.5 1.5]);
grid on
title_str = sprintf('Picosec LED test - Run %s - e-peak amplitude, mean %s',run.id,e_peak_amp.mean);
title(title_str)

figure(2)
hold on;
plot(xbins,polya_cnt_fit,'Linewidth',2);
plot(xbins,noise_polya_cnt_fit,'Linewidth',2);
plot(xbins,polya_cnt_fit-noise_polya_cnt_fit,'Linewidth',2);
legend('Full data fit','Noise only fit','Full data fit - Noise only fit')
%% save data to mat file
save([ run.id '.mat'], 'run');
