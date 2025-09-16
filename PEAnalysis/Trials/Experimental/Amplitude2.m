
glbl_cut = mm_max_y>0.00 & mm_max_y<0.99*max(mm_max_y);

figure
h=histogram(mm_max_y(glbl_cut),100)
hold on

% make a new histogram cut without the noise
indices = find(abs(mm_max_y)<3E-3);
mm_max_y_new = mm_max_y;
mm_max_y_new(indices) = [];
glbl_cut_new = mm_max_y_new>0.00 & mm_max_y_new<0.99*max(mm_max_y_new);
h1 = histogram(mm_max_y_new(glbl_cut_new),100)

xbins = h1.BinEdges(1:end-1)+h1.BinWidth/2;

% remove bins below noise threshold
% indices = find(abs(xbins)<3e-3);
% xbins(indices) = [];
% h1.Values(indices) = [];

fit_data = [];
fit_data(1,:)=1:length(xbins);
fit_data(2,:) = h1.Values;



% for i = 1:length(xbins)
%     if (xbins(i)<3e-3)
%         xbins(i)=[];
%     end
% end

p0=[];
p0(1) = sum(h1.Values)*h1.BinWidth;   % normalization factor
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
%xlim([-0.5 1.5]);
grid on
title_str = sprintf('Picosec LED test - Run %s - e-peak amplitude, mean %s',run.id,e_peak_amp.mean);
title(title_str)


%% save data to mat file
save([ run.id '.mat'], 'run');
