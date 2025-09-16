
addpath 'C:\Users\GDD\Documents\MATLAB\Picosec\SAMPIC\CrateSystem\Functions'

%% enter run IDs to analyse together
runIDsToAnalyse = ["374"];
runTitleString = ['Picolarge Capacitive Sharing - Apr 2023 - Run ' run.id];
loadVariables = true; %should load or plot only 
shouldDisplayMaps = true;
numberMapsToSave = 100;


%% geometry definition 
%plotting
minC =0;
maxC =200;

plotInterval = 10;

pad0Pos = [27.5 26.5]; %center of pad 0 determiend in other run aligned with pad0

padCentersUnit= 0.5*[
    0   0;
    -8.66 15;
    -17.32 0;
    -8.66 -15;
    8.66 -15;
    17.32 0;
    8.66 15;
];

padCenters= 0.5*[
    0   0;
    -8.66 15;
    -17.32 0;
    -8.66 -15;
    8.66 -15;
    17.32 0;
    8.66 15;
];

saveMapCounter = 0;

for r=1:7 %loop through pads
        padNumber = r-1;
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

        run.name = run.id;
        run.oscilloscope = 'Merged';

        variablesFolder = ['\\eosproject-smb\eos\project\p\picosec\testbeam\2023_April_h4\Results\Run' run.id '-Picolarge\variables'];
        store_folderEvents = ['\\eosproject-smb\eos\project\p\picosec\testbeam\2023_April_h4\Results\Run' run.id '-Picolarge\events'];
        mkdir(store_folderEvents);
        store_folder = ['\\eosproject-smb\eos\project\p\picosec\testbeam\2023_April_h4\Results\Run' run.id '-Picolarge'];
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

            if contains(fileInfo.name,'.mat')
                filePath = [variablesFolder '\' fileInfo.name];
                segment = load(filePath);
                numberEventsLoadedSegment = length(segment.MM_data);

                padIDArray = [padIDArray segment.padIDArray];
                time_diff_sigmoid = [time_diff_sigmoid segment.time_diff_sigmoid];
                MCP_maxy = [MCP_maxy segment.MCP_maxy];
                MM_maxy = [MM_maxy segment.MM_maxy];
                trackerX = [trackerX segment.trackerX];
                trackerY = [trackerY segment.trackerY];
                eventIDArray = [eventIDArray segment.eventIDArray];
                MM_DataArray = [MM_DataArray segment.MM_data];
                MCP_DataArray = [MCP_DataArray segment.MCP_data];
                acceptedSignalArray = [acceptedSignalArray segment.acceptedSignalArray];
                
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

glbl_cut = time_diff_sigmoid>time_min & time_diff_sigmoid<time_max & ...
    MCP_maxy>0.01*max(MCP_maxy) & MM_maxy>0.01*max(MM_maxy)&...
    MM_maxy<0.95*max(MM_maxy)& MCP_maxy<0.95*max(MCP_maxy) & ...
    trackerX~=0 & trackerY~=0 & acceptedSignalArray==1;


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

%loop through event IDs to assemble events in structures
for eventPos = 1:length(uniqueEventIDs)
    eventTemp.eventID = uniqueEventIDs(eventPos);


    shouldPlot = false;
    if mod(eventPos,plotInterval)==0
        shouldPlot = true;
    end
    %search for all events with this ID
    eventIDFilter = eventIDArray==eventTemp.eventID & glbl_cut; %& eventIDArray>10000;

    padIDArrayEvent = padIDArray(eventIDFilter);

    if length(padIDArrayEvent)>0
        %valid event
        time_diff_sigmoidEvent = time_diff_sigmoid(eventIDFilter);
        MCP_maxyEvent = MCP_maxy(eventIDFilter);
        MM_maxyEvent = MM_maxy(eventIDFilter);
        trackerXEvent = trackerX(eventIDFilter);
        trackerYEvent = trackerY(eventIDFilter);
        MM_DataArrayEvent = MM_DataArray(eventIDFilter);


%        MCP_DataArrayEvent = MCP_DataArray(eventIDFilter);
        acceptedSignalArrayEvent = acceptedSignalArray(eventIDFilter);
        numberAcceptedSignalsInEvent = length(acceptedSignalArrayEvent);
        MM_maxySum = sum(MM_maxyEvent);

        if max(MM_maxyEvent)>0.000 %show only events with strong signal
                size(MM_DataArrayEvent);
                %VisualisePicolargeEvent(padIDArrayEvent,1000*MM_maxyEvent,['Picolarge event ' num2str(eventTemp.eventID)],[store_folder '\Event' int2str(eventTemp.eventID) '.png'],'Amp','%0.1f mV',minC,maxC,trackerXEvent(1),trackerYEvent(1))
                waveformsArray = [];

                positionRecoX =0;
                positionRecoY =0;

                %loop through hit pads assemble waveforms array
                numChPresent = length(padIDArrayEvent);

                for pos = 1:length(padIDArrayEvent)
                        ampWeight = MM_maxyEvent(pos)/MM_maxySum;
                        
                        padIndex = padIDArrayEvent(pos)+1;

                        posPadX = padCenters(padIndex,1);
                        posPadY = padCenters(padIndex,2);

                        %add up estimate from weights
                        positionRecoX = positionRecoX+posPadX*ampWeight;
                        positionRecoY = positionRecoY+posPadY*ampWeight;
                        
                       waveformTemp.padID = padIDArrayEvent(pos);
                       waveformTemp.x = MM_DataArrayEvent(pos).sig.xfit - MM_DataArrayEvent(pos).sigmoid.timepoint;
                       waveformTemp.riseTime=MM_DataArrayEvent(pos).sigmoid.timepoint90-MM_DataArrayEvent(pos).sigmoid.timepoint10;% extract electrom peak charge

                       waveformTemp.y = MM_DataArrayEvent(pos).sig.yfit;
                       %waveformTemp.xRef = MCP_DataArrayEvent(pos).sig.xfit;
                       %waveformTemp.yRef = MCP_DataArrayEvent(pos).sig.yfit;
                       waveformsArray = [waveformsArray waveformTemp];
                end
                residualX = positionRecoX-trackerXEvent(1);
                residualY = positionRecoY-trackerYEvent(1);
                residualsArrayX = [residualsArrayX;residualX];
                residualsArrayY = [residualsArrayY;residualY];
                eventTrackerXArray = [eventTrackerXArray;trackerXEvent(1)];
                eventTrackerYArray = [eventTrackerYArray;trackerYEvent(1)];
                multiplicityArray = [multiplicityArray;length(padIDArrayEvent)];
                sumAmpArray = [sumAmpArray;MM_maxySum];

                if shouldDisplayMaps & saveMapCounter<numberMapsToSave
                %VisualisePicolargeEventWaveforms(padIDArrayEvent,1000*MM_maxyEvent,waveformsArray,['Picolarge event ' num2str(eventTemp.eventID)],[store_folder '\Event' int2str(eventTemp.eventID) '.png'],'Amp','%0.1f mV',minC,maxC,trackerXEvent(1),trackerYEvent(1))
                VisualisePicolargeEventWaveformsPosReco(padIDArrayEvent,1000*MM_maxyEvent,waveformsArray,[runTitleString ' - Event ' num2str(eventTemp.eventID)],[store_folderEvents '\Event' int2str(eventTemp.eventID) '.png'],'','%0.1f mV',minC,maxC,trackerXEvent(1),trackerYEvent(1),positionRecoX,positionRecoY)
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
hist(multiplicityArray,[1 2 3 4 5 6 7]);
grid on
xlabel('Multiplicity (pads hit)')
ylabel('Events');
title_str = sprintf('%s \n Pad Multiplicity    mean=%.1f    std=%.1f',runTitleString,mean(multiplicityArray), std(multiplicityArray));
title(title_str)
saveas(gcf,[store_folder '\Run' run.id '_multiplicity.png'])

figure;
hist(residualsArrayX,200);
grid on
xlabel('Residuals in X (mm)')
ylabel('Events');
xlim([-10 10]);
title_str = sprintf('%s \n X Residuals    mean=%.2fmm    std=%.2fmm',runTitleString,mean(residualsArrayX),std(residualsArrayX));
title(title_str)
saveas(gcf,[store_folder '\Run' run.id '_residualsX.png'])

figure;
hist(residualsArrayY,200);
grid on
xlabel('Residuals in Y (mm)')
ylabel('Events');
xlim([-10 10]);
title_str = sprintf('%s \n Y Residuals    mean=%.2fmm   std=%.2fmm',runTitleString,mean(residualsArrayY),std(residualsArrayY));
title(title_str)
saveas(gcf,[store_folder '\Run' run.id '_residualsY.png'])

figure;
hist(sumAmpArray,100);
grid on
xlabel('Sum of pad signal amplitudes (mV)')
ylabel('Events');
%xlim([-10 10]);
title_str = sprintf('%s \n Sum of pad signal amplitudes    mean=%.2fV   std=%.2fV',runTitleString,mean(sumAmpArray),std(sumAmpArray));
title(title_str)
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
