close all

store_folder

trackerMask = trackerX>0 & trackerY>0 %& padIDArray==1% & MM_maxy>0.1;
trackerXMasked = trackerX(trackerMask);
trackerYMasked = trackerY(trackerMask);

numberPoints = length(trackerXMasked)

meanX = mean(trackerXMasked)
meanY = mean(trackerYMasked)
figure
hist(trackerXMasked)

xlabel('X position (mm)');
ylabel('Counts');
saveas(gcf,[store_folder '\Run' run.id '_trackerXHist.png'])
movegui(gcf,'northeast');

pause(1);
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

padSide = 0.5*3.5;

colorMap = colormap;

padCenters0= [
    0   0;
    0 3.464;
    3 1.732;
    3 -1.732;
    0 -3.464;
    -3 -1.732;
    -3 1.732;
    0 6.928;
    2.981 5.207;
    6 3.464;
    6 0;
    6 -3.464;
    2.981 -5.207;
    0 -6.928;
    -2.981 -5.207;
    -6 -3.464;
    -6 0;
    -6 3.464;
    -2.981 5.207;
    ];


 pad0Pos = [meanX meanY];
        
        for r=1:19 %loop through pads
            padNumber = r;
            padCenter = padCenters0(r,:);
            padX = padCenter(1) + pad0Pos(1);
            padY = padCenter(2) + pad0Pos(2);

            padCenters(r,1) = padX;
            padCenters(r,2) = padY;
        end
        
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
xlabel('X position (mm)');
ylabel('Y position (mm)');

saveas(gcf,[store_folder '\Run' run.id '_AlignmentMap.png'])
movegui(gcf,'northeast');


