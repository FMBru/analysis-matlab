close all
store_folder

    uniqueChIDs = [];
    for i = 1:length(padIDArray)
        if length(uniqueChIDs)==0 || any(ismember(uniqueChIDs,padIDArray(i))) == 0
            uniqueChIDs = [uniqueChIDs;padIDArray(i)];
        end
    end

for chPos = 1:length(uniqueChIDs)

close all;

trackerMask = trackerX>0 & trackerY>0 & padIDArray==uniqueChIDs(chPos) & MM_maxy>0.1;
trackerXMasked = trackerX(trackerMask);
trackerYMasked = trackerY(trackerMask);

numberPoints = length(trackerXMasked)

meanX = mean(trackerXMasked)
meanY = mean(trackerYMasked)

% hist(trackerXMasked)

% hist(trackerYMasked);

scatter(trackerXMasked,trackerYMasked);

% Normally distributed sample points:
x = trackerXMasked;
y = trackerYMasked;

% Bin the data:
pts = linspace(20, 50, 101);
N = histcounts2(y(:), x(:), pts, pts);

% Plot scattered data (for comparison):
%subplot(1, 2, 1);
%scatter(x, y, 'r.');
axis equal;
set(gca, 'XLim', pts([1 end]), 'YLim', pts([1 end]));

% Plot heatmap:
%subplot(1, 2, 2);
imagesc(pts, pts, N);
axis equal;
set(gca, 'XLim', pts([1 end]), 'YLim', pts([1 end]), 'YDir', 'normal');
hold on


padSide = 0.5*3.5;

colorMap = colormap;
str_title = sprintf('ch %d', uniqueChIDs(chPos));
title(str_title);


for r=1:19 %loop through pads
    padNumber = r;
    padCenter = padCenters(r,:);
    padX = padCenter(1);
    padY = padCenter(2);
    %check if value for hitmap exists
    channelID = padNumber;
    padNumber;
    channelID;
    hexagonOutlineRotated(padSide,padX,padY);
    hold on
end
pause(1);
saveas(gcf,[store_folder '\Run' run.id '-ch' num2str(uniqueChIDs(chPos)) '_Hits.png'])

end


