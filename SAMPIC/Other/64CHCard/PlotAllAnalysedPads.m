%take padResults vector and visualise

padIDArray = [];
padCHArray = [];
padXArray = [];
padYArray = [];
RMSSamplingArray = [];
SATSamplingArray = [];
MCPAmpSamplingArray = [];
MMAmpSamplingArray = [];
NumberEntriesSamplingArray = [];

for resultPos = 1:length(padResults)
    
    acceptPad = 0;
    
    if padResults(resultPos).rmsSampling>0 && padResults(resultPos).numberEntriesSampling>20
        acceptPad=1;
    end
    
    if acceptPad
    padIDArray = [padIDArray;padResults(resultPos).padID];
    padCHArray = [padCHArray;getChannelForPadNumber(padResults(resultPos).padID)];
    padXArray = [padXArray;100-padResults(resultPos).xc-60];
    padYArray = [padYArray;100-padResults(resultPos).yc-53];
    RMSSamplingArray = [RMSSamplingArray;padResults(resultPos).rmsSampling];
    SATSamplingArray = [SATSamplingArray;padResults(resultPos).SATMeanSampling];
    MCPAmpSamplingArray = [MCPAmpSamplingArray;padResults(resultPos).mcpAmpSampling];
    MMAmpSamplingArray = [MMAmpSamplingArray;padResults(resultPos).mmAmpSampling];
    NumberEntriesSamplingArray = [NumberEntriesSamplingArray;padResults(resultPos).numberEntriesSampling];
    end
    end

store_folder = ['\\eosproject-smb\eos\project\p\picosec\testbeam\2022_October_h4\Results\Run' run.id '-' run.oscilloscope];

close all

%% plot area maps
figure
for resultPos = length(padResults):-1:1
    acceptPad=0;
    
    if  padResults(resultPos).rmsSampling>0 && padResults(resultPos).numberEntriesSampling>20
        acceptPad=1;
    end
    if acceptPad
        xMatrix = padResults(resultPos).areaMaps.xx+54;
        yMatrix = padResults(resultPos).areaMaps.yy+49;

        h=pcolor(xMatrix,yMatrix,padResults(resultPos).areaMaps.rms*1000);
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
caxis([0 50]);
grid on
h = colorbar;
h.Label.String = 'Time resolution, ps';
h.Label.Position(1) = 3;
str_title = ['Time resolution - Run' run.id '-' run.oscilloscope];
title(str_title);
saveas(gcf,[store_folder '\Run' run.id '_TimeResMap.png'])

%% plot area maps
figure
for resultPos = length(padResults):-1:1
    acceptPad=0;
    
    if  padResults(resultPos).rmsSampling>0 && padResults(resultPos).numberEntriesSampling>20
        acceptPad=1;
    end
    if acceptPad
        xMatrix = padResults(resultPos).areaMaps.xx+54;
        yMatrix = padResults(resultPos).areaMaps.yy+49;

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
figure
for resultPos = length(padResults):-1:1
    acceptPad=0;
    
    if  padResults(resultPos).rmsSampling>0 && padResults(resultPos).numberEntriesSampling>20
        acceptPad=1;
    end
    if acceptPad
        xMatrix = padResults(resultPos).areaMaps.xx+54;
        yMatrix = padResults(resultPos).areaMaps.yy+49;

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
figure
for resultPos = length(padResults):-1:1
    acceptPad=0;
    
    if  padResults(resultPos).rmsSampling>0 && padResults(resultPos).numberEntriesSampling>20
        acceptPad=1;
    end
    if acceptPad
        xMatrix = padResults(resultPos).areaMaps.xx+54;
        yMatrix = padResults(resultPos).areaMaps.yy+49;

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





VisualiseMultipad(padCHArray,NumberEntriesSamplingArray,'Hitmap Multipad Picosec',[store_folder '\Run' run.id '_NumberAnalysedSamplingArea.png'],'Hits: %d');

VisualiseMultipad(padCHArray,RMSSamplingArray,'RMS Sampling Area Multipad Picosec',[store_folder '\Run' run.id '_RMSSamplingArea.png'],'RMS: %0.1f ps');
VisualiseMultipad(padCHArray,abs(SATSamplingArray),'SAT Multipad Picosec',[store_folder '\Run' run.id '_SATSamplingArea.png'],'SAT: %0.1f ns');
VisualiseMultipad(padCHArray,MCPAmpSamplingArray,'MCP Amp Sampling Area Multipad Picosec',[store_folder '\Run' run.id '_MCPAmpSamplingArea.png'],'Amp: %0.2f V');
VisualiseMultipad(padCHArray,MMAmpSamplingArray,'MM Amp Sampling Area Multipad Picosec',[store_folder '\Run' run.id '_MMAmpSamplingArea.png'],'Amo: %0.2f V');

VisualiseMultipad(padCHArray,padXArray,'X Position Multipad Picosec',[store_folder '\Run' run.id '_XPos.png'],'X: %0.1f mm');
VisualiseMultipad(padCHArray,padYArray,'Y Position Area Multipad Picosec',[store_folder '\Run' run.id '_YPos.png'] ,'Y: %0.1f mm');

