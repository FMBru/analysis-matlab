
close all


addpath 'C:\Users\gdd.CERN\Documents\MATLAB\Picosec\SAMPIC\CrateSystem\Functions'
%addpath 'C:\Users\GDD\Documents\MATLAB\Picosec\SAMPIC\CrateSystem\Functions'

addpath '\\eosproject-smb\eos\project\p\picosec\testbeam\2025_July_h4\Matlab\SAMPIC\CrateSystem\MediumGranularity\Functions'
padSide = 0.5*3.5;



%% enter run IDs to analyse together
runIDsToAnalyse = ["206"];
loadVariables = true ; %should load or plot only
shouldDisplayMaps = true;
numberMapsToSave = 20;


%% geometry definition
%plotting
minC = 0;
maxC = 200;

plotInterval = 10;

%pad 19 is at -5.25 9.093 compared to pad 0
%pad 19 on tracker is at 36.0182 / 36.4170 (determined with matlab scrip[t
%loadTrackerFileAlignment
%-> pad 0 on tracker is

pad0Pos = [33.26 27.9]; %center of pad 0 determiend in other run aligned with pad0
pad0Pos = [33 31.4]; %from mean of all hits on central pad
pad0Pos = [30.71 59.66]; %July25 run

% padCenters= [
%     0   0;
%     0 6.062;
%     5.25 3.031;
%     5.25 -3.031;
%     0 -6.062;
%     -5.25 -3.031;
%     -5.25 3.031;
%     0 12.12;
%     5.25 9.093;
%     10.5 6.062;
%     10.5 0;
%     10.5 -6.062;
%     5.25 -9.093;
%     0 -12.12;
%     -5.25 -9.093;
%     -10.5 -6.062;
%     -10.5 0;
%     -10.5 6.062;
%     -5.25 9.093;
% ];
padCenters= [
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



saveMapCounter = 0;

for r=1:19 %loop through pads
    padNumber = r;
    padCenter = padCenters(r,:);
    padX = padCenter(1) + pad0Pos(1);
    padY = padCenter(2) + pad0Pos(2);

    padCenters(r,1) = padX;
    padCenters(r,2) = padY;
end





if loadVariables
    clear MM_data MCP_data time_diff time_diff_sigmoid MCP_maxy MM_maxy trackerX trackerY eventIDArray padIDArray;
    padIDArray = [];
    time_diff_sigmoid = [];
    MCP_maxy = [];
    MM_maxy = [];
    trackerX = [];
    trackerY = [];
    eventIDArray = [];
    MM_epeak = [];
    MCP_epeak = [];
    MM_riseTime = [];
    MM_bgAvg = [];
    MM_bgRMS = [];
    MM_DataArray = [];
    MCP_DataArray = [];
    acceptedSignalArray = [];

    for runPos = 1:length(runIDsToAnalyse)

        run.id = convertStringsToChars(runIDsToAnalyse(runPos));
        runTitleString = ['Medium Granularity - June 2024 - Run ' run.id];
        run.name = run.id;
        run.oscilloscope = 'Merged';

        %variablesFolder = ['C:\Users\GDD\Documents\Picosec\Run' run.id '-SAMPIC\variables'];
        variablesFolder = ['\\eosproject-smb\eos\project\p\picosec\testbeam\2025_July_h4\Results\Run' run.id '-SAMPIC\variables'];
            variablesFolder = ['F:\Processing\Picosec\Run' run.id '-SAMPIC\variables'];

        store_folderEvents = ['\\eosproject-smb\eos\project\p\picosec\testbeam\2025_July_h4\Results\Run' run.id '-SAMPIC\events'];
        mkdir(store_folderEvents);
        store_folder = ['\\eosproject-smb\eos\project\p\picosec\testbeam\2025_July_h4\Results\Run' run.id '-SAMPIC'];
        mkdir(store_folder);


        str_disp=sprintf('Loading variables for run %s', run.id);
        disp(str_disp);

        variablesFilesList = dir(variablesFolder);
        numberFiles = length(variablesFilesList);

        str_disp=sprintf('Found %d files to load', numberFiles);
        disp(str_disp);

        for i=1:length(variablesFilesList)
            fileInfo = variablesFilesList(i);

            str_disp=sprintf('Loading file %d / %d', i, numberFiles);
            disp(str_disp);

    if contains(fileInfo.name,'.mat') & contains(fileInfo.name,'processedBatch')
                filePath = [variablesFolder '\' fileInfo.name];
                segment = load(filePath);
                numberEventsLoadedSegment = length(segment.MM_data);

                padIDArray = [padIDArray segment.dutChannelArray];
                time_diff_sigmoid = [time_diff_sigmoid segment.time_diff_sigmoid];
                MCP_maxy = [MCP_maxy segment.MCP_maxy];
                MM_maxy = [MM_maxy segment.MM_maxy];
                trackerX = [trackerX segment.trackerX];
                trackerY = [trackerY segment.trackerY];
                eventIDArray = [eventIDArray segment.eventIDArray];
                MM_DataArray = [MM_DataArray segment.MM_data];
                MCP_DataArray = [MCP_DataArray segment.MCP_data];
                % acceptedSignalArray = [acceptedSignalArray segment.acceptedSignalArray];

                %         for pos = 1:numberEventsLoadedSegment
                %             %step through loaded file and assemble arrays
                %             MM_epeakTemp(pos) = segment.MM_data(pos).sig.charge.e_peak;
                %             MCP_epeakTemp(pos)= segment.MCP_data(pos).sig.charge.e_peak;% extract electrom peak charge
                %             MM_riseTimeTemp(pos)= segment.MM_data(pos).sigmoid.timepoint90- segment.MM_data(pos).sigmoid.timepoint10;% extract electrom peak charge
                %         end
                %         MM_epeak = [MM_epeak MM_epeakTemp];
                %         MCP_epeak = [MCP_epeak MCP_epeakTemp];
                %         MM_riseTime = [MM_riseTime MM_riseTimeTemp];
                %         clear MM_epeakTemp MCP_epeakTemp MM_riseTimeTemp

                clear segment
            end
        end

        str_disp=sprintf('Loaded %d hits', length(padIDArray));
        disp(str_disp);
    end

    %assemble events

    %determine unique eventIDs
    uniqueEventIDs = [];
    for i = 1:length(eventIDArray)
        if length(uniqueEventIDs)==0 || any(ismember(uniqueEventIDs,eventIDArray(i))) == 0
            uniqueEventIDs = [uniqueEventIDs;eventIDArray(i)];
        end
    end
end


%determine global filter
time_avg_raw = median(time_diff_sigmoid(MM_maxy>0.01));  % assume mu from the median vaule
time_min = time_avg_raw - 10.5;     % predicted resolution 100ps cut 3 sigma
time_max = time_avg_raw + 10.5;     % left and 3 sigma right from median

%MM amplitude cuts to accept hit - in V
min_accept = 0.01;
max_accept = 0.95;

glbl_cut = time_diff_sigmoid>time_min & time_diff_sigmoid<time_max & ...
    MCP_maxy>min_accept*max(MCP_maxy) & MM_maxy>min_accept*max(MM_maxy)&...
    MM_maxy<max_accept*max(MM_maxy)& MCP_maxy<max_accept*max(MCP_maxy) & ...
    trackerX~=0 & trackerY~=0 & ~isnan(trackerX) & ~isnan(trackerY); %& acceptedSignalArray==1;


% glbl_cut = time_diff_sigmoid>time_min & time_diff_sigmoid<time_max & ...
%     trackerX~=0 & trackerY~=0;
%


eventsArray = [];
residualsArrayX = [];
residualsArrayY = [];
eventTrackerXArray = [];
eventTrackerYArray = [];
multiplicityArray = [];
sumAmpArray = [];
check = 0;

%loop through event IDs to assemble events in structures
for eventPos = 1:length(uniqueEventIDs)
    eventTemp.eventID = uniqueEventIDs(eventPos);

    shouldPlot = false;
    if mod(eventPos,plotInterval)==0
        shouldPlot = true;
    end

    %search for all events with this ID
    eventIDFilter = eventIDArray==eventTemp.eventID & glbl_cut; %& eventIDArray>10000;
    str_disp=sprintf('Plotting event ID %d', eventTemp.eventID);
    %disp(str_disp);;

    chIDArrayEvent = padIDArray(eventIDFilter);

    if length(chIDArrayEvent)>0
        %valid event
        time_diff_sigmoidEvent = time_diff_sigmoid(eventIDFilter);
        MCP_maxyEvent = MCP_maxy(eventIDFilter);
        MM_maxyEvent = MM_maxy(eventIDFilter);



        trackerXEvent = trackerX(eventIDFilter);
        trackerYEvent = trackerY(eventIDFilter);
        MM_DataArrayEvent = MM_DataArray(eventIDFilter);

        %invert X relative to center
        trackerXEvent = -(trackerXEvent-pad0Pos(1))+pad0Pos(1);

        %        MCP_DataArrayEvent = MCP_DataArray(eventIDFilter);
        %acceptedSignalArrayEvent = acceptedSignalArray(eventIDFilter);
        numberAcceptedSignalsInEvent = length(MCP_maxyEvent)
        MM_maxySum = sum(MM_maxyEvent);

        %% loop through events - select specific events here
        % Constraint on the positions to consider -> avoid bad pads
        if max(MM_maxySum)>0.25 && numberAcceptedSignalsInEvent>3 && (trackerXEvent(1)<padCenters(6,1) || trackerXEvent(1)>padCenters(4,1)) && trackerYEvent(1)<(padCenters(2,2)) %show only events with strong signal
            size(MM_DataArrayEvent);

            waveformsArray = [];

            positionRecoX = 0;
            positionRecoY = 0;
            chIDArrayEvent = unique(chIDArrayEvent);
            %loop through hit pads assemble waveforms array
            numChPresent = length(chIDArrayEvent);
            %padIDArrayEvent

            chIDArrayEvent;
            padIDArrayEvent = chIDArrayEvent;

            for pos = 1:length(chIDArrayEvent)
                ampWeight = MM_maxyEvent(pos)/MM_maxySum;
                chNumber = chIDArrayEvent(pos); %+1;

                padIndex = getPadForChannelNumberMediumGranularity(chNumber);
                padIDArrayEvent (pos) = padIndex;

                if padIndex <=0
                    check = check + 1;
                    break
                end
                padIndex;
                posPadX = padCenters(padIndex,1)
                posPadY = padCenters(padIndex,2)

                %add up estimate from weights
                positionRecoX = positionRecoX+posPadX*ampWeight;
                positionRecoY = positionRecoY+posPadY*ampWeight;

                waveformTemp.padID = padIDArrayEvent(pos);
                waveformTemp.x = MM_DataArrayEvent(pos).waveform(:,1) - MM_DataArrayEvent(pos).sigmoid.timepoint;
                waveformTemp.riseTime=MM_DataArrayEvent(pos).sigmoid.timepoint90-MM_DataArrayEvent(pos).sigmoid.timepoint10;% extract electrom peak charge

                waveformTemp.y = MM_DataArrayEvent(pos).waveform(:,2);
                %waveformTemp.xRef = MCP_DataArrayEvent(pos).sig.xfit;
                %waveformTemp.yRef = MCP_DataArrayEvent(pos).sig.yfit;
                waveformsArray = [waveformsArray waveformTemp];
            end
            %                 positionRecoX
            %                 positionRecoY
            %                 chIDArrayEvent

            residualX = positionRecoX-trackerXEvent(1);
            residualY = positionRecoY-trackerYEvent(1);
            residualsArrayX = [residualsArrayX;residualX];
            residualsArrayY = [residualsArrayY;residualY];
            eventTrackerXArray = [eventTrackerXArray;trackerXEvent(1)];
            eventTrackerYArray = [eventTrackerYArray;trackerYEvent(1)];
            multiplicityArray = [multiplicityArray;length(chIDArrayEvent)];
            sumAmpArray = [sumAmpArray;MM_maxySum];

            if shouldDisplayMaps & saveMapCounter<numberMapsToSave
                VisualiseMediumGranularityEventWaveformsPosReco(padCenters,padIDArrayEvent,1000*MM_maxyEvent,waveformsArray,[runTitleString ' - Event ' num2str(eventTemp.eventID)],[store_folderEvents '\Event' int2str(eventTemp.eventID) '.png'],'','%0.1f mV',minC,maxC,trackerXEvent(1),trackerYEvent(1),positionRecoX,positionRecoY)
                hold off;
                pause(1);
                close all;

                saveMapCounter = saveMapCounter+1;
            end
        end
    end
end

%plot multiplicity
figure;
hist(multiplicityArray,[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19]);
grid on
xlabel('Multiplicity (pads hit)')
ylabel('Events');
title_str = sprintf('%s \n Pad Multiplicity    mean=%.1f    std=%.1f',runTitleString,mean(multiplicityArray), std(multiplicityArray));
title(title_str)
movegui(gcf,'south');

saveas(gcf,[store_folder '\Run' run.id '_multiplicity.png'])

figure;
hist(residualsArrayX,200);
grid on
xlabel('Residuals in X (mm)')
ylabel('Events');
%xlim([-10 10]);
title_str = sprintf('%s \n X Residuals    mean=%.2fmm    std=%.2fmm',runTitleString,mean(residualsArrayX),std(residualsArrayX));
title(title_str)
movegui(gcf,'north');

saveas(gcf,[store_folder '\Run' run.id '_residualsX.png'])

figure;
hist(residualsArrayY,200);
grid on
xlabel('Residuals in Y (mm)')
ylabel('Events');
%xlim([-10 10]);
title_str = sprintf('%s \n Y Residuals    mean=%.2fmm   std=%.2fmm',runTitleString,mean(residualsArrayY),std(residualsArrayY));
title(title_str)
movegui(gcf,'east');

saveas(gcf,[store_folder '\Run' run.id '_residualsY.png'])

figure;
hist(sumAmpArray,100);
grid on
xlabel('Sum of pad signal amplitudes (mV)')
ylabel('Events');
%xlim([-10 10]);
title_str = sprintf('%s \n Sum of pad signal amplitudes    mean=%.2fV   std=%.2fV',runTitleString,mean(sumAmpArray),std(sumAmpArray));
title(title_str)
movegui(gcf,'west');

saveas(gcf,[store_folder '\Run' run.id '_sumAmp.png'])


%% save results arrays
storeMatfilePath = [store_folder '\results.mat'];
m = matfile(storeMatfilePath,'Writable',true);

save(storeMatfilePath,'residualsArrayX');
save(storeMatfilePath,'residualsArrayY','-append');
save(storeMatfilePath,'eventTrackerXArray','-append');
save(storeMatfilePath,'eventTrackerYArray','-append');
save(storeMatfilePath,'multiplicityArray','-append');
save(storeMatfilePath,'sumAmpArray','-append');



%% 2D maps

%sampling area for 2D maps
area.step = 0.5;   % set grid resolution in mm
area.size = 10;      % set half-size of the observed square area, mm - for
area.radius = 1;    % set radius of averaging circular window cut, mm
area.x_vec =-area.size:area.step:area.size+area.step;    % define area of interest for 2D plots in X axis
area.y_vec =-area.size:area.step:area.size+area.step;    % define area of interest for 2D plots in X axis

% move rectangle to pad center
area.x_vec = area.x_vec + pad0Pos(1) - area.step/2;
area.y_vec = area.y_vec + pad0Pos(2) - area.step/2;

% ndgrid plot
[area.xx, area.yy] = ndgrid(area.x_vec, area.y_vec);
for i=1:length(area.x_vec)
    for j=1:length(area.y_vec)
        % make moving circular cut that respects the global cut (idx_cut)
        cut_circ = ((eventTrackerXArray - area.x_vec(i)).^2 + (eventTrackerYArray - area.y_vec(j)).^2) < area.radius^2; %& glbl_cut;
        area.resXMean(i,j) = mean(residualsArrayX(cut_circ));
        area.resXStd(i,j) = std(residualsArrayX(cut_circ));
        area.resYMean(i,j) = mean(residualsArrayY(cut_circ));
        area.resYStd(i,j) = std(residualsArrayY(cut_circ));
    end
end


%% plot x residuals as map
figure
h=pcolor(area.xx,area.yy,area.resXMean);
hold on
set(h, 'EdgeColor', 'none');
axis equal
colorbar
xlabel('x-axis, mm');
ylabel('y-axis, mm');
h = colorbar;
h.Label.String = 'Residual X Mean';
h.Label.Position(1) = 3;
str_title = sprintf('%s: ResidualX-Mean', runTitleString);
title(str_title);
for r=1:19 %loop through pads
    padNumber = r;
    padCenter = padCenters(r,:);
    padX = padCenter(1);
    padY = padCenter(2);
    hexagonOutlineRotated(padSide,padX,padY);
end
axis equal
saveas(gcf,[store_folder '\Run' run.id '_ResidualXMean-Map.png'])
movegui(gcf,'northeast');
pause(1);

%% plot x residuals Std as map
figure
h=pcolor(area.xx,area.yy,area.resXStd);
hold on
set(h, 'EdgeColor', 'none');
axis equal
colorbar
xlabel('x-axis, mm');
ylabel('y-axis, mm');
h = colorbar;
h.Label.String = 'Residual X Std';
h.Label.Position(1) = 3;
str_title = sprintf('%s: ResidualX Std', runTitleString);
title(str_title);
for r=1:19 %loop through pads
    padNumber = r;
    padCenter = padCenters(r,:);
    padX = padCenter(1);
    padY = padCenter(2);
    hexagonOutlineRotated(padSide,padX,padY);
end
axis equal
saveas(gcf,[store_folder '\Run' run.id '_ResidualXStd-Map.png'])
movegui(gcf,'northeast');
pause(1);


%% plot y residuals as map
figure
h=pcolor(area.xx,area.yy,area.resYMean);
hold on
set(h, 'EdgeColor', 'none');
axis equal
colorbar
xlabel('x-axis, mm');
ylabel('y-axis, mm');
h = colorbar;
h.Label.String = 'Residual Y Mean';
h.Label.Position(1) = 3;
str_title = sprintf('%s: ResidualY Mean', runTitleString);
title(str_title);
for r=1:19 %loop through pads
    padNumber = r;
    padCenter = padCenters(r,:);
    padX = padCenter(1);
    padY = padCenter(2);
    hexagonOutlineRotated(padSide,padX,padY);
end
axis equal
saveas(gcf,[store_folder '\Run' run.id '_ResidualYMean-Map.png'])
movegui(gcf,'northeast');
pause(1);

%% plot x residuals Std as map
figure
h=pcolor(area.xx,area.yy,area.resYStd);
hold on
set(h, 'EdgeColor', 'none');
axis equal
colorbar
xlabel('x-axis, mm');
ylabel('y-axis, mm');
h = colorbar;
h.Label.String = 'Residual Y Std';
h.Label.Position(1) = 3;
str_title = sprintf('%s: ResidualY Std', runTitleString);
title(str_title);
for r=1:19 %loop through pads
    padNumber = r;
    padCenter = padCenters(r,:);
    padX = padCenter(1);
    padY = padCenter(2);
    hexagonOutlineRotated(padSide,padX,padY);
end
axis equal
saveas(gcf,[store_folder '\Run' run.id '_ResidualYStd-Map.png'])
movegui(gcf,'northeast');
pause(1);

