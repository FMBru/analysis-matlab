%take padResults vector and visualise

%%load pad results from previously processed run
%run.id = '354';
run.oscilloscope='SAMPIC';
load(['\\eosproject-smb\eos\project\p\picosec\testbeam\2025_July_h4\Results\Run' run.id '-' run.oscilloscope '\padResults.mat']);


requireMinHitsToAccept = 0; %set to 0 if accepting all pads regardless of hit number

addpath '\\eosproject-smb\eos\project\p\picosec\testbeam\2025_July_h4\Matlab\SAMPIC\CrateSystem\Functions';

runTitleString = ['Run' run.id '-' run.oscilloscope] ;

padIDArray = [];
padCHArray = [];
padXArray = [];
padYArray = [];
RMSSamplingArray = [];
SATSamplingArray = [];
bgAvgSamplingArray = [];
bgRMSSamplingArray = [];
MCPAmpSamplingArray = [];
MMAmpSamplingArray = [];
NumberEntriesSamplingArray = [];
NumberEntriesFullPadArray = [];
NumberEntriesFullPadArrayGlblCut = [];
riseTimeSamplingArray = [];
xCenterNominal = [];
yCenterNominal = [];
areaMapsArray = [];
twCorrArray = [];

store_folder = ['\\eosproject-smb\eos\project\p\picosec\testbeam\2025_July_h4\Results\Run' run.id '-' run.oscilloscope];


%xOffset = 54;
%yOffset = 49;

xOffset = 0;
yOffset = 0;

%Calculate expected position for each pad

%xy coords for strarting point for coordinate of pad 1
startXCenter = 50;
startYCenter = 50;




for resultPos = 1:length(padResults)
    
    acceptPad = 0;
    
    if padResults(resultPos).padID>0 && padResults(resultPos).numberEntriesSampling>requireMinHitsToAccept && padResults(resultPos).rmsSampling>0 && padResults(resultPos).numberEntriesSampling>20
        acceptPad=1;
    end
    
    
    
    if acceptPad
        padIDArray = [padIDArray;padResults(resultPos).padID];
        padCHArray = [padCHArray;getChannelForPadNumber(padResults(resultPos).padID)];
        padXArray = [padXArray;100-padResults(resultPos).xc-xOffset];
        padYArray = [padYArray;100-padResults(resultPos).yc-yOffset];
        RMSSamplingArray = [RMSSamplingArray;padResults(resultPos).rmsSampling];
        SATSamplingArray = [SATSamplingArray;padResults(resultPos).SATMeanSampling];
        MCPAmpSamplingArray = [MCPAmpSamplingArray;padResults(resultPos).mcpAmpSampling];
        MMAmpSamplingArray = [MMAmpSamplingArray;padResults(resultPos).mmAmpSampling];
        NumberEntriesSamplingArray = [NumberEntriesSamplingArray;padResults(resultPos).numberEntriesSampling];
        NumberEntriesFullPadArray = [NumberEntriesFullPadArray;padResults(resultPos).numberEntriesFullPad];
        NumberEntriesFullPadArrayGlblCut = [NumberEntriesFullPadArrayGlblCut;padResults(resultPos).numberEntriesFullPadGlblCut];
        
        bgAvgSamplingArray = [bgAvgSamplingArray;padResults(resultPos).bgAvg];
        bgRMSSamplingArray = [bgRMSSamplingArray;padResults(resultPos).bgRMS];
        riseTimeSamplingArray = [riseTimeSamplingArray;padResults(resultPos).riseTime];
        areaMapsArray = [areaMapsArray;padResults(resultPos).areaMaps];
        twCorrArray = [twCorrArray;padResults(resultPos).twCorr];
        
        
    end
end

referencePadForStartingPoint = 33;
%determine starting point for residuals
for resultPos = 1:length(padIDArray)
    if padIDArray(resultPos) == referencePadForStartingPoint
        %save as starting point
        startXCenter = padXArray(resultPos);
        startYCenter = padYArray(resultPos);
    end
end

%find residuals for all pads
xResiduals = [];
yResiduals = [];
for pos = 1:length(padIDArray)
    %nominal offset from reference pad
    padIDArray(pos);
    colMain = mod(padIDArray(pos), 10);
    rowMain = (padIDArray(pos) - colMain) / 10;
    
    colTest = mod(referencePadForStartingPoint, 10);
    rowTest = (referencePadForStartingPoint - colTest) / 10;
    
    diffCol = colMain-colTest;
    diffRow = rowMain-rowTest;
    
    if colMain==0
        diffCol = 10-colTest-colMain;
        diffRow = diffRow-1;
    end
    
    estimatePosX = startXCenter+diffCol*10;
    estimatePosY = startYCenter+diffRow*10;
    
    xResiduals = [xResiduals; padXArray(pos)-estimatePosX];
    yResiduals = [yResiduals; padYArray(pos)-estimatePosY];
    
end

figure
hist(xResiduals,100)
xlabel('X position residual, mm');
ylabel('Events');
grid on
str_title = ['Xresidual distribution - Run' run.id '-' run.oscilloscope];
title(str_title);
saveas(gcf,[store_folder '\Run' run.id '_xResiduals.png'])

figure
hist(yResiduals,100)
xlabel('Y position residual, mm');
ylabel('Events');
grid on
str_title = ['Yresidual distribution - Run' run.id '-' run.oscilloscope];
title(str_title);
saveas(gcf,[store_folder '\Run' run.id '_yResiduals.png'])

xMeanResidual = median(xResiduals);
yMeanResidual = median(yResiduals);

refXCenter = startXCenter+xMeanResidual;
refYCenter = startYCenter+yMeanResidual;

%find residuals for all pads
xResiduals = [];
yResiduals = [];
for pos = 1:length(padIDArray)
    %nominal offset from reference pad
    padIDArray(pos);
    colMain = mod(padIDArray(pos), 10);
    rowMain = (padIDArray(pos) - colMain) / 10;
    
    colTest = mod(referencePadForStartingPoint, 10);
    rowTest = (referencePadForStartingPoint - colTest) / 10;
    
    diffCol = colMain-colTest;
    diffRow = rowMain-rowTest;
    
    if colMain==0
        diffCol = 10-colTest-colMain;
        diffRow = diffRow-1;
    end
    
    estimatePosX = refXCenter+diffCol*10;
    estimatePosY = refYCenter+diffRow*10;
    
    
    
    xResiduals = [xResiduals; padXArray(pos)-estimatePosX];
    yResiduals = [yResiduals; padYArray(pos)-estimatePosY];
    
    xCenterNominal = [xCenterNominal; estimatePosX];
    yCenterNominal = [yCenterNominal; estimatePosY];
    
end

figure
hist(xResiduals,100)
xlabel('X position residual correctedAlignment, mm');
ylabel('Events');
grid on
str_title = ['Xresidual distribution - Run' run.id '-' run.oscilloscope];
title(str_title);
saveas(gcf,[store_folder '\Run' run.id '_xResidualsCorr.png'])

figure
hist(yResiduals,100)
xlabel('Y position residual correctedAlignment, mm');
ylabel('Events');
grid on
str_title = ['Yresidual distribution - Run' run.id '-' run.oscilloscope];
title(str_title);
saveas(gcf,[store_folder '\Run' run.id '_yResidualsCorr.png'])

close all



%% plot  3D plot
close all
figure
scatterbar3(xCenterNominal,yCenterNominal,MMAmpSamplingArray,10);
hold on
grid on
%axis equal
colorbar
xlabel('x-axis, mm');
ylabel('y-axis, mm');
zlabel('Amplitude, V');
xlim([-10 110]);
ylim([-10 110]);
zlim([0 0.4]);
view(3);
str_title = ['MM Amp - Run' run.id '-' run.oscilloscope];
title(str_title);
saveas(gcf,[store_folder '\Run' run.id '_MMamp-3D.png'])
savefig(gcf,[store_folder '\Run' run.id '_MMamp-3D.fig'])


%% plot  3D plot
close all
figure
scatterbar3(xCenterNominal,yCenterNominal,MCPAmpSamplingArray,10);
hold on
grid on
%axis equal
colorbar
xlabel('x-axis, mm');
ylabel('y-axis, mm');
zlabel('Amplitude, V');
xlim([-10 110]);
ylim([-10 110]);
zlim([0 0.4]);
view(3);
str_title = ['MCP Amp - Run' run.id '-' run.oscilloscope];
title(str_title);
saveas(gcf,[store_folder '\Run' run.id '_MCPamp-3D.png'])
savefig(gcf,[store_folder '\Run' run.id '_MCPamp-3D.fig'])


%% plot  3D plot
close all
figure
scatterbar3(xCenterNominal,yCenterNominal,RMSSamplingArray,10);
hold on
grid on
%axis equal
colorbar
xlabel('x-axis, mm');
ylabel('y-axis, mm');
zlabel('RMS (ps)');
xlim([-10 110]);
ylim([-10 110]);
%zlim([10 50]);
view(3);
str_title = ['RMS - Run' run.id '-' run.oscilloscope];
title(str_title);
saveas(gcf,[store_folder '\Run' run.id '_RMS-3D.png'])
savefig(gcf,[store_folder '\Run' run.id '_RMS-3D.fig'])


%% plot  3D plot
close all
figure
scatterbar3(xCenterNominal,yCenterNominal,abs(SATSamplingArray),10);
hold on
%axis equal
grid on
colorbar
xlabel('x-axis, mm');
ylabel('y-axis, mm');
zlabel('SAT (ns)');
xlim([-10 110]);
ylim([-10 110]);
%zlim([22 30]);
view(3);
str_title = ['SAT - Run' run.id '-' run.oscilloscope];
title(str_title);
saveas(gcf,[store_folder '\Run' run.id '_SAT-3D.png'])
savefig(gcf,[store_folder '\Run' run.id '_SAT-3D.fig'])







%% plot area maps - cutting pad to right size from nominal position
close all
figure
for resultPos = 1:length(areaMapsArray)
    xMatrix = areaMapsArray(resultPos).xx+xOffset;
    yMatrix = areaMapsArray(resultPos).yy+yOffset;
    
    xCenterNominalPad = 100-xCenterNominal(resultPos);
    yCenterNominalPad = 100-yCenterNominal(resultPos);
    
    xMinPad = xCenterNominalPad-5;
    xMaxPad = xCenterNominalPad+5;
    
    yMinPad = yCenterNominalPad-5;
    yMaxPad = yCenterNominalPad+5;
    xMask = xMatrix>xMinPad & xMatrix<xMaxPad;
    yMask = yMatrix>yMinPad & yMatrix<yMaxPad;
    combinedMask = xMask & yMask;
    
    padData = areaMapsArray(resultPos).amp;
    
    padData(isnan(padData))=0;
    
    xMatrixMasked = xMatrix.*combinedMask;
    yMatrixMasked = yMatrix.*combinedMask;
    dataMatrixMasked = padData.*combinedMask;
    
    xMatrixMasked( all(~xMatrixMasked,2), : ) = [];
    xMatrixMasked( :, all(~xMatrixMasked,1) ) = [];
    dataMatrixMasked( all(~yMatrixMasked,2), : ) = [];
    dataMatrixMasked( :, all(~yMatrixMasked,1) ) = [];
    
    yMatrixMasked( all(~yMatrixMasked,2), : ) = [];
    yMatrixMasked( :, all(~yMatrixMasked,1) ) = [];
    padData(isnan(padData))=0;
    
    dataMatrixMasked(dataMatrixMasked==0) = NaN;
    
    if max(max(combinedMask)) == 1
        h=pcolor(xMatrixMasked,yMatrixMasked,dataMatrixMasked);
        set(h, 'EdgeColor', 'none');
    end
    hold on;
end
axis equal
colorbar
xlabel('x-axis, mm');
ylabel('y-axis, mm');
xlim([0 100])
ylim([0 100])
%caxis([0 0.4]);
grid on
h = colorbar;
h.Label.String = 'Signal amplitude (V)';
h.Label.Position(1) = 3;
str_title = ['MM Amplitude - Run' run.id '-' run.oscilloscope];
title(str_title);
saveas(gcf,[store_folder '\Run' run.id '_AmpMapCroppedPads.png'])







close all
figure
for resultPos = 1:length(areaMapsArray)
    xMatrix = areaMapsArray(resultPos).xx+xOffset;
    yMatrix = areaMapsArray(resultPos).yy+yOffset;
    
    xCenterNominalPad = 100-xCenterNominal(resultPos);
    yCenterNominalPad = 100-yCenterNominal(resultPos);
    
    xMinPad = xCenterNominalPad-5;
    xMaxPad = xCenterNominalPad+5;
    
    yMinPad = yCenterNominalPad-5;
    yMaxPad = yCenterNominalPad+5;
    xMask = xMatrix>xMinPad & xMatrix<xMaxPad;
    yMask = yMatrix>yMinPad & yMatrix<yMaxPad;
    combinedMask = xMask & yMask;
    
    padData = areaMapsArray(resultPos).rms*1000;
    
    padData(isnan(padData))=0;
    
    xMatrixMasked = xMatrix.*combinedMask;
    yMatrixMasked = yMatrix.*combinedMask;
    dataMatrixMasked = padData.*combinedMask;
    
    xMatrixMasked( all(~xMatrixMasked,2), : ) = [];
    xMatrixMasked( :, all(~xMatrixMasked,1) ) = [];
    dataMatrixMasked( all(~yMatrixMasked,2), : ) = [];
    dataMatrixMasked( :, all(~yMatrixMasked,1) ) = [];
    
    yMatrixMasked( all(~yMatrixMasked,2), : ) = [];
    yMatrixMasked( :, all(~yMatrixMasked,1) ) = [];
    padData(isnan(padData))=0;
    
    dataMatrixMasked(dataMatrixMasked==0) = NaN;
    
    if max(max(combinedMask)) == 1
        h=pcolor(xMatrixMasked,yMatrixMasked,dataMatrixMasked);
        set(h, 'EdgeColor', 'none');
    end
    hold on;
end
axis equal
colorbar
xlabel('x-axis, mm');
ylabel('y-axis, mm');
xlim([0 100])
ylim([0 100])
%caxis([10 50]);
grid on
h = colorbar;
h.Label.String = 'Timing resolution RMS (ps)';
h.Label.Position(1) = 3;
str_title = ['Timing Resolution - Run' run.id '-' run.oscilloscope];
title(str_title);
saveas(gcf,[store_folder '\Run' run.id '_RMSMapCroppedPads.png'])



close all
figure
for resultPos = 1:length(areaMapsArray)
    xMatrix = areaMapsArray(resultPos).xx+xOffset;
    yMatrix = areaMapsArray(resultPos).yy+yOffset;
    
    xCenterNominalPad = 100-xCenterNominal(resultPos);
    yCenterNominalPad = 100-yCenterNominal(resultPos);
    
    xMinPad = xCenterNominalPad-5;
    xMaxPad = xCenterNominalPad+5;
    
    yMinPad = yCenterNominalPad-5;
    yMaxPad = yCenterNominalPad+5;
    xMask = xMatrix>xMinPad & xMatrix<xMaxPad;
    yMask = yMatrix>yMinPad & yMatrix<yMaxPad;
    combinedMask = xMask & yMask;
    
    padData = abs(areaMapsArray(resultPos).sat);
    
    padData(isnan(padData))=0;
    
    xMatrixMasked = xMatrix.*combinedMask;
    yMatrixMasked = yMatrix.*combinedMask;
    dataMatrixMasked = padData.*combinedMask;
    
    xMatrixMasked( all(~xMatrixMasked,2), : ) = [];
    xMatrixMasked( :, all(~xMatrixMasked,1) ) = [];
    dataMatrixMasked( all(~yMatrixMasked,2), : ) = [];
    dataMatrixMasked( :, all(~yMatrixMasked,1) ) = [];
    
    yMatrixMasked( all(~yMatrixMasked,2), : ) = [];
    yMatrixMasked( :, all(~yMatrixMasked,1) ) = [];
    padData(isnan(padData))=0;
    
    dataMatrixMasked(dataMatrixMasked==0) = NaN;
    
    if max(max(combinedMask)) == 1
        h=pcolor(xMatrixMasked,yMatrixMasked,dataMatrixMasked);
        set(h, 'EdgeColor', 'none');
    end
    hold on;
end
axis equal
colorbar
xlabel('x-axis, mm');
ylabel('y-axis, mm');
xlim([0 100])
ylim([0 100])
%caxis([22 30]);
grid on
h = colorbar;
h.Label.String = 'SAT (ns)';
h.Label.Position(1) = 3;
str_title = ['SAT - Run' run.id '-' run.oscilloscope];
title(str_title);
saveas(gcf,[store_folder '\Run' run.id '_SATMapCroppedPads.png'])


%%
close all
figure
for resultPos = 1:length(areaMapsArray)
    xMatrix = areaMapsArray(resultPos).xx+xOffset;
    yMatrix = areaMapsArray(resultPos).yy+yOffset;
    
    xCenterNominalPad = 100-xCenterNominal(resultPos);
    yCenterNominalPad = 100-yCenterNominal(resultPos);
    
    xMinPad = xCenterNominalPad-5;
    xMaxPad = xCenterNominalPad+5;
    
    yMinPad = yCenterNominalPad-5;
    yMaxPad = yCenterNominalPad+5;
    xMask = xMatrix>xMinPad & xMatrix<xMaxPad;
    yMask = yMatrix>yMinPad & yMatrix<yMaxPad;
    combinedMask = xMask & yMask;
    
    padData = 1000*(areaMapsArray(resultPos).sat-twCorrArray(resultPos,1));
    
    padData(isnan(padData))=0;
    
    xMatrixMasked = xMatrix.*combinedMask;
    yMatrixMasked = yMatrix.*combinedMask;
    dataMatrixMasked = padData.*combinedMask;
    
    xMatrixMasked( all(~xMatrixMasked,2), : ) = [];
    xMatrixMasked( :, all(~xMatrixMasked,1) ) = [];
    dataMatrixMasked( all(~yMatrixMasked,2), : ) = [];
    dataMatrixMasked( :, all(~yMatrixMasked,1) ) = [];
    
    yMatrixMasked( all(~yMatrixMasked,2), : ) = [];
    yMatrixMasked( :, all(~yMatrixMasked,1) ) = [];
    padData(isnan(padData))=0;
    
    dataMatrixMasked(dataMatrixMasked==0) = NaN;
    
    if max(max(combinedMask)) == 1
        h=pcolor(xMatrixMasked,yMatrixMasked,dataMatrixMasked);
        set(h, 'EdgeColor', 'none');
    end
    hold on;
end
axis equal
colorbar
xlabel('x-axis, mm');
ylabel('y-axis, mm');
xlim([0 100])
ylim([0 100])
%caxis([-100 100]);
grid on
h = colorbar;
h.Label.String = 'Mean Corrected SAT (ps)';
h.Label.Position(1) = 3;
str_title = ['Mean Corrected SAT- Run' run.id '-' run.oscilloscope];
title(str_title);
saveas(gcf,[store_folder '\Run' run.id '_SATMap-CorrectedMean-CroppedPads.png'])





        
        %%
        
close all
figure
for resultPos = 1:length(areaMapsArray)
    xMatrix = areaMapsArray(resultPos).xx+xOffset;
    yMatrix = areaMapsArray(resultPos).yy+yOffset;
    
    xCenterNominalPad = 100-xCenterNominal(resultPos);
    yCenterNominalPad = 100-yCenterNominal(resultPos);
    
    xMinPad = xCenterNominalPad-5;
    xMaxPad = xCenterNominalPad+5;
    
    yMinPad = yCenterNominalPad-5;
    yMaxPad = yCenterNominalPad+5;
    xMask = xMatrix>xMinPad & xMatrix<xMaxPad;
    yMask = yMatrix>yMinPad & yMatrix<yMaxPad;
    combinedMask = xMask & yMask;
    
    padData = areaMapsArray(resultPos).amp_MCP;
    
    padData(isnan(padData))=0;
    
    xMatrixMasked = xMatrix.*combinedMask;
    yMatrixMasked = yMatrix.*combinedMask;
    dataMatrixMasked = padData.*combinedMask;
    
    xMatrixMasked( all(~xMatrixMasked,2), : ) = [];
    xMatrixMasked( :, all(~xMatrixMasked,1) ) = [];
    dataMatrixMasked( all(~yMatrixMasked,2), : ) = [];
    dataMatrixMasked( :, all(~yMatrixMasked,1) ) = [];
    
    yMatrixMasked( all(~yMatrixMasked,2), : ) = [];
    yMatrixMasked( :, all(~yMatrixMasked,1) ) = [];
    padData(isnan(padData))=0;
    
    dataMatrixMasked(dataMatrixMasked==0) = NaN;
    
    if max(max(combinedMask)) == 1
        h=pcolor(xMatrixMasked,yMatrixMasked,dataMatrixMasked);
        set(h, 'EdgeColor', 'none');
    end
    hold on;
end
axis equal
colorbar
xlabel('x-axis, mm');
ylabel('y-axis, mm');
xlim([0 100])
ylim([0 100])
%caxis([0 0.4]);
grid on
h = colorbar;
h.Label.String = 'Signal amplitude (V)';
h.Label.Position(1) = 3;
str_title = ['MCP Amplitude - Run' run.id '-' run.oscilloscope];
title(str_title);
saveas(gcf,[store_folder '\Run' run.id '_AmpMapMCPCroppedPads.png'])




pause(1);
%% plot area maps
close all

figure
for resultPos = length(padResults):-1:1
    acceptPad=0;
    
    if  padResults(resultPos).rmsSampling>0 && padResults(resultPos).numberEntriesSampling>20
        acceptPad=1;
    end
    if acceptPad
        xMatrix = padResults(resultPos).areaMaps.xx+xOffset;
        yMatrix = padResults(resultPos).areaMaps.yy+yOffset;
        
        h=pcolor(xMatrix,yMatrix,padResults(resultPos).areaMaps.amp);
        set(h, 'EdgeColor', 'none');
    end
    if resultPos==length(padResults)
        hold on
        
    end
end
set(h, 'EdgeColor', 'none');
axis equal
colorbar
xlabel('x-axis, mm');
ylabel('y-axis, mm');
xlim([0 100])
ylim([0 100])
%caxis([0 0.4]);
grid on
h = colorbar;
h.Label.String = 'Signal amplitude (V)';
h.Label.Position(1) = 3;
str_title = ['MM Amplitude - Run' run.id '-' run.oscilloscope];
title(str_title);
saveas(gcf,[store_folder '\Run' run.id '_AmpMap.png'])

pause(1);
%% plot area maps
close all

figure
for resultPos = length(padResults):-1:1
    acceptPad=0;
    
    if  padResults(resultPos).rmsSampling>0 && padResults(resultPos).numberEntriesSampling>20
        acceptPad=1;
    end
    if acceptPad
        xMatrix = padResults(resultPos).areaMaps.xx+xOffset;
        yMatrix = padResults(resultPos).areaMaps.yy+yOffset;
        
        h=pcolor(xMatrix,yMatrix,padResults(resultPos).areaMaps.amp);
        set(h, 'EdgeColor', 'none');
    end
    if resultPos==length(padResults)
        hold on
        
    end
end
set(h, 'EdgeColor', 'none');
axis equal
colorbar
xlabel('x-axis, mm');
ylabel('y-axis, mm');
xlim([0 100])
ylim([0 100])
grid on
h = colorbar;
h.Label.String = 'Signal amplitude (V)';
h.Label.Position(1) = 3;
str_title = ['MM Amplitude - Run' run.id '-' run.oscilloscope];
title(str_title);
saveas(gcf,[store_folder '\Run' run.id '_AmpMap.png'])




%% plot area maps
close all

figure
for resultPos = length(padResults):-1:1
    acceptPad=0;
    
    if  padResults(resultPos).rmsSampling>0 && padResults(resultPos).numberEntriesSampling>20
        acceptPad=1;
    end
    if acceptPad
        xMatrix = padResults(resultPos).areaMaps.xx+xOffset;
        yMatrix = padResults(resultPos).areaMaps.yy+yOffset;
        
        h=pcolor(xMatrix,yMatrix,padResults(resultPos).areaMaps.sat);
        set(h, 'EdgeColor', 'none');
    end
    if resultPos==length(padResults)
        hold on
        
    end
end
set(h, 'EdgeColor', 'none');
axis equal
colorbar
xlabel('x-axis, mm');
ylabel('y-axis, mm');
xlim([0 100])
ylim([0 100])
grid on
h = colorbar;
h.Label.String = 'SAT, ns';
h.Label.Position(1) = 3;

str_title = ['SAT - Run' run.id '-' run.oscilloscope];
title(str_title);
saveas(gcf,[store_folder '\Run' run.id '_SATMap.png'])


%% plot area maps
close all

figure
for resultPos = length(padResults):-1:1
    acceptPad=0;
    
    if  padResults(resultPos).rmsSampling>0 && padResults(resultPos).numberEntriesSampling>20
        acceptPad=1;
    end
    if acceptPad
        xMatrix = padResults(resultPos).areaMaps.xx+xOffset;
        yMatrix = padResults(resultPos).areaMaps.yy+yOffset;
        
        h=pcolor(xMatrix,yMatrix,padResults(resultPos).areaMaps.amp_MCP);
        set(h, 'EdgeColor', 'none');
    end
    if resultPos==length(padResults)
        hold on
        
    end
end
set(h, 'EdgeColor', 'none');
axis equal
colorbar
xlabel('x-axis, mm');
ylabel('y-axis, mm');
xlim([0 100])
ylim([0 100])
grid on
h = colorbar;
h.Label.String = 'MCP amplitude, V';
h.Label.Position(1) = 3;
str_title = ['MCP Amplitude - Run' run.id '-' run.oscilloscope];
title(str_title);

saveas(gcf,[store_folder '\Run' run.id '_MCPAmpMap.png'])

close all

figure
hist(MMAmpSamplingArray,100)
xlabel('MM amp in sampling area, V');
ylabel('Events');
grid on
xlim([0 0.4]);
title_str = sprintf('%s \n MM amp distribution \\mu = %3.3f V',runTitleString, mean(MMAmpSamplingArray));
title(title_str);
saveas(gcf,[store_folder '\Run' run.id '_MMAmpSampling-DistributionFull.png'])

figure
hist(MCPAmpSamplingArray,100)
xlabel('MCP amp in sampling area, V');
ylabel('Events');
grid on
xlim([0 0.4]);
title_str = sprintf('%s \n MCP amp distribution \\mu = %3.3f V',runTitleString, mean(MCPAmpSamplingArray));
title(title_str);
saveas(gcf,[store_folder '\Run' run.id '_MCPAmpSampling-DistributionFull.png'])

figure
hist(MCPAmpSamplingArray,100); hold on
hist(MMAmpSamplingArray,100)

xlabel('Amp in sampling area, V');
ylabel('Events');
grid on
xlim([0 0.4]);
title_str = sprintf('%s \n Amp distribution \\mu = %3.3f V',runTitleString, mean(MCPAmpSamplingArray));
title(title_str);
saveas(gcf,[store_folder '\Run' run.id '_MMMCPAmpSampling-DistributionFull.png'])




figure
hist(RMSSamplingArray,100)
xlabel('Timing resolution RMS in sampling area, ps');
ylabel('Events');
grid on
title_str = sprintf('%s \n Timing resolution distribution \\mu = %3.3f ps',runTitleString, mean(RMSSamplingArray));
title(title_str);
saveas(gcf,[store_folder '\Run' run.id '_RMSSampling-DistributionFull.png'])

figure
hist(RMSSamplingArray,100)
xlabel('Timing resolution RMS in sampling area, ps');
ylabel('Events');
xlim([0 50])
grid on
title_str = sprintf('%s \n Timing resolution distribution \\mu = %3.3f ps',runTitleString, mean(RMSSamplingArray));
title(str_title);
saveas(gcf,[store_folder '\Run' run.id '_RMSSampling-DistributionZoom.png'])

figure
hist(abs(SATSamplingArray),100)
xlabel('SAT in sampling area, ns');
ylabel('Events');
xlim([20 30])
grid on
title_str = sprintf('%s \n SAT  distribution \\mu = %3.3f ns',runTitleString, mean(SATSamplingArray));
title(title_str);
saveas(gcf,[store_folder '\Run' run.id '_SATSampling-Distribution.png'])


figure
scatter(MMAmpSamplingArray,RMSSamplingArray,'.')
xlabel('Mean amplitude in sampling area, V');
ylabel('Timing resolution in sampling area RMS, ps');
xlim([0 0.3])
ylim([15 100])
grid on
str_title = ['RMS vs. MeanAmp - Run' run.id '-' run.oscilloscope];
title(str_title);
saveas(gcf,[store_folder '\Run' run.id '_RMSvsAmp.png'])

figure
scatter(NumberEntriesSamplingArray,MMAmpSamplingArray,'.')
ylabel('Mean amplitude in sampling area, V');
xlabel('Events in sampling area');
grid on
str_title = ['Events in sampling area vs. MeanAmp - Run' run.id '-' run.oscilloscope];
title(str_title);
saveas(gcf,[store_folder '\Run' run.id '_EventCountvsAmp.png'])

bgAvgSamplingArray = bgAvgSamplingArray*1000;
bgRMSSamplingArray = bgRMSSamplingArray*1000;


figure
hist(bgAvgSamplingArray,20)
xlabel('Baseline level, mV');
ylabel('Events');
grid on
xlim([-10 10])
title_str = sprintf('%s \n Baseline level distribution \\mu = %3.3f mV',runTitleString, mean(bgAvgSamplingArray));
title(title_str);
saveas(gcf,[store_folder '\Run' run.id '_BaselineLevel.png'])

figure
hist(bgRMSSamplingArray,20)
xlabel('Baseline RMS, mV');
ylabel('Events');
xlim([0 30]);
grid on
title_str = sprintf('%s \n Baseline RMS distribution \\mu = %3.3f mV',runTitleString, mean(bgRMSSamplingArray));
title(title_str);
saveas(gcf,[store_folder '\Run' run.id '_BaselineRMS.png'])

figure
hist(riseTimeSamplingArray,100)
xlabel('Rise time, ns')
ylabel('Events');
xlim([0 1]);
grid on
title_str = sprintf('%s \n Rise time distribution \\mu = %3.3f mV',runTitleString, mean(riseTimeSamplingArray));
title(title_str);
saveas(gcf,[store_folder '\Run' run.id '_riseTimeDistribution.png'])


%with fixed limits
% VisualiseMultipad(padCHArray,NumberEntriesSamplingArray,['Run ' run.id ' - Hits in Sampling Area'],[store_folder '\Run' run.id '_NumberAnalysedSamplingArea.png'],'Hits', '%d',0,0);
% VisualiseMultipad(padCHArray,NumberEntriesFullPadArray,['Run ' run.id ' - Hits in Full Pad'],[store_folder '\Run' run.id '_NumberAnalysedFullPad.png'],'Hits', '%d',0,0);
% VisualiseMultipad(padCHArray,NumberEntriesFullPadArrayGlblCut,['Run ' run.id ' - Hits in Full Pad after global cut'],[store_folder '\Run' run.id '_NumberAnalysedFullPadGlblCut.png'],'Hits', '%d',0,0);
% VisualiseMultipad(padCHArray,RMSSamplingArray,['Run ' run.id ' - RMS in Sampling Area'],[store_folder '\Run' run.id '_RMSSamplingArea.png'],'RMS', '%0.1f ps',0,0);
% VisualiseMultipad(padCHArray,abs(SATSamplingArray),['Run ' run.id ' - SAT in Sampling Area'],[store_folder '\Run' run.id '_SATSamplingArea.png'],'SAT', '%0.1f ns',24,30);
% VisualiseMultipad(padCHArray,MCPAmpSamplingArray,['Run ' run.id ' - MCP Amp in Sampling Area'],[store_folder '\Run' run.id '_MCPAmpSamplingArea.png'],'Amp', '%0.2f V',0,0.4);
% VisualiseMultipad(padCHArray,MMAmpSamplingArray,['Run ' run.id ' - MM Amp in Sampling Area'],[store_folder '\Run' run.id '_MMAmpSamplingArea.png'],'Amp', '%0.2f V',0,0.4);
% VisualiseMultipad(padCHArray,padXArray,['Run ' run.id ' - X Position'],[store_folder '\Run' run.id '_XPos.png'],'X', '%0.1f mm',0,100);
% VisualiseMultipad(padCHArray,padYArray,['Run ' run.id ' - Y Position Area'],[store_folder '\Run' run.id '_YPos.png'] ,'Y', '%0.1f mm',0,100);
% VisualiseMultipad(padCHArray,abs(bgAvgSamplingArray),['Run ' run.id ' - Baseline Level'],[store_folder '\Run' run.id '_bgAvg.png'] ,'BG', '%0.1f mV',0,0);
% VisualiseMultipad(padCHArray,abs(bgRMSSamplingArray),['Run ' run.id ' - Baseline RMS'],[store_folder '\Run' run.id '_bgRMS.png'] ,'BG RMS', '%0.1f mV',0,0);
% VisualiseMultipad(padCHArray,riseTimeSamplingArray,['Run ' run.id ' - Rise time'],[store_folder '\Run' run.id '_riseTime.png'] ,'RT', '%0.2f ns',0.3,0.8);

%no fixed limits
VisualiseMultipad(padCHArray,NumberEntriesSamplingArray,['Run ' run.id ' - Hits in Sampling Area'],[store_folder '\Run' run.id '_NumberAnalysedSamplingArea.png'],'Hits', '%d',0,0);
VisualiseMultipad(padCHArray,RMSSamplingArray,['Run ' run.id ' - RMS in Sampling Area'],[store_folder '\Run' run.id '_RMSSamplingArea.png'],'RMS', '%0.1f ps',0,0);
VisualiseMultipad(padCHArray,abs(SATSamplingArray),['Run ' run.id ' - SAT in Sampling Area'],[store_folder '\Run' run.id '_SATSamplingArea.png'],'SAT', '%0.1f ns',0,0);
VisualiseMultipad(padCHArray,MCPAmpSamplingArray,['Run ' run.id ' - MCP Amp in Sampling Area'],[store_folder '\Run' run.id '_MCPAmpSamplingArea.png'],'Amp', '%0.2f V',0,0);
VisualiseMultipad(padCHArray,MMAmpSamplingArray,['Run ' run.id ' - MM Amp in Sampling Area'],[store_folder '\Run' run.id '_MMAmpSamplingArea.png'],'Amp', '%0.2f V',0,0);
VisualiseMultipad(padCHArray,padXArray,['Run ' run.id ' - X Position'],[store_folder '\Run' run.id '_XPos.png'],'X', '%0.1f mm',0,0);
VisualiseMultipad(padCHArray,padYArray,['Run ' run.id ' - Y Position Area'],[store_folder '\Run' run.id '_YPos.png'] ,'Y', '%0.1f mm',0,0);
VisualiseMultipad(padCHArray,abs(bgAvgSamplingArray),['Run ' run.id ' - Baseline Level'],[store_folder '\Run' run.id '_bgAvg.png'] ,'BG', '%0.1f mV',0,0);
VisualiseMultipad(padCHArray,abs(bgRMSSamplingArray),['Run ' run.id ' - Baseline RMS'],[store_folder '\Run' run.id '_bgRMS.png'] ,'BG RMS', '%0.1f mV',0,0);
VisualiseMultipad(padCHArray,riseTimeSamplingArray,['Run ' run.id ' - Rise time'],[store_folder '\Run' run.id '_riseTime.png'] ,'RT', '%0.2f ns',0,0);

figure
%%plot timewalk correction curves together
padIDArray = [];

for resultPos = 1:length(padResults)
    
    acceptPad = 0;
    
    if padResults(resultPos).rmsSampling>0 && padResults(resultPos).numberEntriesSampling>20
        acceptPad=1;
    end
    
    if acceptPad
        padIDArray = [padIDArray;padResults(resultPos).padID];
        plot(padResults(resultPos).twCorrEpeak,padResults(resultPos).twCorrSAT-padResults(resultPos).twCorr(1),'.'); hold on
    end
end
hold off
xlabel('Electron peak charge')
ylabel('SAT (offset subtracted), ns');
grid on
str_title = ['Rise time distribution - Run' run.id '-' run.oscilloscope];
title(str_title);
saveas(gcf,[store_folder '\Run' run.id '_TWcomparison.png'])

