close all

store_folder

trackerMask = trackerX>0 & trackerY>0 & padIDArray==24 & MM_maxy>0.1;
trackerXMasked = trackerX(trackerMask);
trackerYMasked = trackerY(trackerMask);

numberPoints = length(trackerXMasked)

meanX = mean(trackerXMasked)
meanY = mean(trackerYMasked)

pad1Pos = [94.0599 97.2475];
padLength = 10;
padHeight=padLength;

padCenters = [];
for c=1:10
    for r=1:10
        padX = pad1Pos(1) - (r-1)*padLength;
        padY = pad1Pos(2) - (c-1)*padHeight;
        padNumber = (10-r+1)+(10-c)*10;
        %check if value for hitmap exists
        channelID = getChannelForPadNumber(padNumber,run.id);
        padCenters = [padCenters; padX padY];
        
    end
end

figure
hist(trackerXMasked)

xlabel('X position (mm)');
ylabel('Counts');
saveas(gcf,[store_folder '\Run' run.id '_trackerXHist.png'])
movegui(gcf,'northeast');

pause(1);
figure
hist(trackerYMasked);

xlabel('Y position (mm)');
ylabel('Counts');
saveas(gcf,[store_folder '\Run' run.id '_trackerYHist.png'])
movegui(gcf,'northeast');

figure
scatter(trackerXMasked,trackerYMasked);
saveas(gcf,[store_folder '\Run' run.id '_ScatterMap.png'])
movegui(gcf,'northeast');

% Normally distributed sample points:
x = trackerXMasked;
y = trackerYMasked;

% Bin the data:
pts = linspace(20, 100, 101);
N = histcounts2(y(:), x(:), pts, pts);

% Plot scattered data (for comparison):
%subplot(1, 2, 1);
%scatter(x, y, 'r.');
axis equal;
set(gca, 'XLim', pts([1 end]), 'YLim', pts([1 end]));

figure
% Plot heatmap:
%subplot(1, 2, 2);
imagesc(pts, pts, N);
axis equal;
set(gca, 'XLim', pts([1 end]), 'YLim', pts([1 end]), 'YDir', 'normal');
hold on

padSide = 10;

colorMap = colormap;

for r=1:length(padCenters) %loop through pads
    padNumber = r;
    padCenter = padCenters(r,:);
    padX = padCenter(1);
    padY = padCenter(2);
    %check if value for hitmap exists
    channelID = padNumber;
    padNumber;
    channelID;
    rectangle('Position',[(padX-padLength/2) (padY-padLength/2) padLength padHeight]);
    hold on
end

xlabel('X position (mm)');
ylabel('Y position (mm)');
saveas(gcf,[store_folder '\Run' run.id '_AlignmentMap.png'])
movegui(gcf,'northeast');


