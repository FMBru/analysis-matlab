function [meanValue, stdValue] = distribution(values,nBins, nameQuantity, unityOfMeas, store_folder, titleString, shouldSave)
%DISTRIBUTION plot distribution of 1 quantity

figure;
histogram(values, nBins);
xlabel([nameQuantity ' (' unityOfMeas ')']);
ylabel('Events');
title([nameQuantity ' distribution' titleString]);

meanValue = mean(values(:), 'omitnan');
stdValue = std(values(:), 'omitnan');
message = sprintf('\\mu: %.3f\n', meanValue);
message = [message sprintf('RMS: %.3f', stdValue)];
y_pos=get(gca,'ylim');
x_pos=get(gca,'xlim');
text(x_pos(1)+0.1*(x_pos(2)-x_pos(1)),y_pos(1)+0.8*(y_pos(2)-y_pos(1)),message);
if shouldSave
    saveas(gcf,[store_folder '\' nameQuantity 'Distribution' titleString '.png']);
end
hold off
close

end

