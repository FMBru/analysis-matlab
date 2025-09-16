function [meanPolya, stdPolya] = ampDistribution(valuesMCP,valuesMM, titleString, store_folder, shouldSave)
%AMPDISTRIBUTION to plot distributions of MM and MCP on the same graph and fit with a Polya 

figure;
histogram(valuesMCP, 100);
hold on
h = histogram(valuesMM, 100);
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
evalc("[p, err, chi] = fminuit('polya_minuit',p0,fit_data(:,1:end),'-b','-c',cmd);");
polya_cnt_fit = polya_minuit(p,fit_data(1,:));
plot(xbins,polya_cnt_fit,'Linewidth',2);
meanPolya = sum(polya_cnt_fit.*xbins)/sum(polya_cnt_fit);
stdPolya = sqrt(sum(polya_cnt_fit.*((xbins-meanPolya).^2)) / sum(polya_cnt_fit));
meanerrPolya = sqrt(sum(polya_cnt_fit.*((xbins-meanPolya).^2)))/sum(polya_cnt_fit);
[dummy, maxPolya] = max(polya_cnt_fit);
legend('MCP', 'MM', 'Polya fit');
xlabel('Amplitude (V)');
ylabel('Events');
title(['Amplitude distribution' titleString ' - Mean: ' num2str(meanPolya,2) ' \pm ' num2str(meanerrPolya,2)]);
hold off
if shouldSave
    saveas(gcf,[store_folder '\amplitudeDistribution' titleString '.png']);
end
close

end

