close all
cut = 20;
cutUp = 300;

h=histogram(mm_max_y(glbl_cut),500)

x0 = [10 2 0.01]; 
fitfun = fittype( @(N,theta,nBar,x) (N./nBar).*((theta+1).^(theta+1)).*((x./nBar).^theta).*exp(-(theta+1).*x./nBar)./gamma(theta-1));
[fitted_curve,gof] = fit(h.BinEdges(cut:end-1-cutUp)',h.Values(cut:end-cutUp)',fitfun,'StartPoint',x0)

hold on
plot(h.BinEdges(cut:end-1-cutUp),fitted_curve(h.BinEdges(cut:end-1-cutUp)),'LineWidth',3)
set(gca, 'YScale', 'log')

ch2 = sum(((h.Values(cut:end-cutUp)-fitted_curve(h.BinEdges(cut:end-1-cutUp))').^2)./fitted_curve(h.BinEdges(cut:end-1-cutUp))');
dof = size(h.Values(cut:end-cutUp),2);
nch2 = ch2/dof
np = 1-chi2cdf(ch2,dof) % P(\chi^2>ch2)

text(0.0,max(h.Values)+200,'\chi^2/dof = '+string(round(ch2)) + '/' + string(size(h.Values(cut:end-cutUp),2)) + ' = ' + string(nch2),'FontSize',14)

title('Amplitude spectrum')
xlabel('Max signal value [V]') 
ylabel('Entries [1]')