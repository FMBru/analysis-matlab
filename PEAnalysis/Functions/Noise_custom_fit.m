%% histogram and polya fit

glbl_cut = mm_max_y>0.00 & mm_max_y<0.95*max(mm_max_y); %remove saturated datapoints
noise_glbl_cut = noise_mm_max_y>0.00 & noise_mm_max_y<0.99*max(noise_mm_max_y);


close all

num_bins = 500;
figure(1)
h=histogram(mm_max_y(glbl_cut),num_bins)
set(gca, 'YScale', 'log') %plot on log scale
hold on


min_v = 3.2e-3;
max_v = 100;

% If min and max for fit are not within dataset range pick edges of dataset
if (min(h.BinEdges) > min_v)
    min_v = min(h.BinEdges)
end
if (max(h.BinEdges(1:num_bins)) < max_v)
    max_v = max(h.BinEdges(1:num_bins)) % cut the last bin edge because it's on the right side of the bin to make array lengths match
end

% min cut value should be the closest value to min_v that is greater than
% min_v
[cut_val,cut_idx] = min(abs(h.BinEdges-min_v));
if (cut_val >= min_v)
else
    cut_idx=cut_idx+1;
end
% max cut value should be the closest value to max_v that is less than
% max_v
[cutUp_val,cutUp_idx] = min(abs(h.BinEdges-max_v));
if (cutUp_val <= max_v)
else
    cutUp_idx=cutUp_idx-1;
end

%cutUp = 200;


% define polya and set starting parameters for fitting
%x0 = [10 2 0.01]; 
x0 = [10 1 0.01];
fitfun = fittype( @(N,theta,nBar,x) (N./nBar).*((theta+1).^(theta+1)).*((x./nBar).^theta).*exp(-(theta+1).*x./nBar)./gamma(theta-1));
[fitted_curve,gof] = fit(h.BinEdges(cut_idx:cutUp_idx)',h.Values(cut_idx:cutUp_idx)',fitfun,'StartPoint',x0)

% plot the fitted curve
hold on
plot(h.BinEdges(cut_idx:cutUp_idx),fitted_curve(h.BinEdges(cut_idx:cutUp_idx)),'LineWidth',3)


% calculate chi squared and degrees of freedom
ch2 = sum(((h.Values(cut_idx:cutUp_idx)-fitted_curve(h.BinEdges(cut_idx:cutUp_idx))').^2)./fitted_curve(h.BinEdges(cut_idx:cutUp_idx))');
dof = size(h.Values(cut_idx:cutUp_idx),2);
nch2 = ch2/dof
np = 1-chi2cdf(ch2,dof) % P(\chi^2>ch2)

% annotate and label plot
%text(min(h.BinEdges)+0.02*max(h.BinEdges),max(h.Values)+2,'\chi^2/dof = '+string(round(ch2)) + '/' + string(size(h.Values(cut_idx:cutUp_idx),2)) + ' = ' + string(nch2),'FontSize',14)
%text();
title_str = sprintf('Amplitude spectrum - Mean = %s - chi^2/dof = %s',fitted_curve.nBar,string(nch2))
title(title_str)
xlabel('Max signal value [V]') 
ylabel('Entries [1]')

%%
