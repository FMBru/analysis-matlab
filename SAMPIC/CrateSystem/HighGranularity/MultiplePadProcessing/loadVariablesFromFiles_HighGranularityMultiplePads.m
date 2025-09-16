
close all
clear all

%addpath 'C:\Users\gdd.CERN\Documents\MATLAB\Picosec\SAMPIC\CrateSystem\Functions'
addpath '\\eosproject-smb\eos\project\p\picosec\testbeam\2024_September_h4\Matlab\SAMPIC\CrateSystem\HighGranularity\Functions'
addpath '\\eosproject-smb\eos\project\p\picosec\testbeam\2024_September_h4\Matlab\SAMPIC\CrateSystem\Functions'
addpath '\\eosproject-smb\eos\project\p\picosec\testbeam\2024_September_h4\Matlab\CommonFunctions'
addpath '\\eosproject-smb\eos\project\p\picosec\testbeam\2024_September_h4\Matlab\SAMPIC\CrateSystem\HighGranularity\MultiplePadProcessing'
padSide = 0.5*3.5;

%% enter run IDs to analyse together
runIDsToAnalyse = ["296"];
loadVariables = true; %should load or plot only
shouldDisplayMaps = true;
numberMapsToSave = 50;

channelsEnabled = [];

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
pad0Pos = [35.42 33.88]; %from mean of all hits on central pad
pad0Pos = [32.44 32.03]; %high granularity


padCenters= [
    0   0;
    0 2.1651;
    1.875 1.0825;
    1.875 -1.0825;
    0 -2.1651;
    -1.875 -1.0825;
    -1.875 1.0825;
    0 4.3302;
    1.875 3.2477;
    3.75 2.1651;
    3.75 0;
    3.75 -2.1651;
    1.875 -3.2477;
    0 -4.3302;
    -1.875 -3.2477;
    -3.75 -2.1651;
    -3.75 0;
    -3.75 2.1651;
    -1.875 3.2477;
    0 6.4953;
    1.875 5.4127;
    3.75 4.3302;
    5.625 3.2477
    5.625 1.0826
    5.625 -1.0826
    5.625 -3.2477
    3.75 -4.3302;
    1.875 -5.4127;
    0 -6.4953;
    -1.875 -5.4127;
    -3.75 -4.3302;
    -5.625 -3.2477
    -5.625 -1.0826
    -5.625 1.0826
    -5.625 3.2477
    -3.75 4.3302;
    -1.875 5.4127;
    ];

saveMapCounter = 0;

for padNumber=1:37 %loop through pads
    padNumber = padNumber;
    padCenter = padCenters(padNumber,:);
    padX = padCenter(1) + pad0Pos(1);
    padY = padCenter(2) + pad0Pos(2);

    padCenters(padNumber,1) = padX;
    padCenters(padNumber,2) = padY;
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
        runTitleString = ['High Granularity - September 2024 - Run ' run.id];
        run.name = run.id;
        run.oscilloscope = 'Merged';

        variablesFolder = ['C:\Users\gdd.CERN\Documents\Picosec\Run' run.id '-SAMPIC\variables'];
        % variablesFolder = ['\\eosproject-smb\eos\project\p\picosec\testbeam\2024_September_h4\Results\Run' run.id '-SAMPIC\variables'];
        store_folderEvents = ['\\eosproject-smb\eos\project\p\picosec\testbeam\2024_September_h4\Results\Run' run.id '-SAMPIC\events'];
        mkdir(store_folderEvents);
        store_folder = ['\\eosproject-smb\eos\project\p\picosec\testbeam\2024_September_h4\Results\Run' run.id '-SAMPIC'];
        mkdir(store_folder);
        store_folderPads = ['\\eosproject-smb\eos\project\p\picosec\testbeam\2024_September_h4\Results\Run' run.id '-SAMPIC\Pads'];
        mkdir(store_folderPads);


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

            if contains(fileInfo.name,'.mat') && contains(fileInfo.name,'processedBatch')

                filePath = [variablesFolder '\' fileInfo.name];
                segment = load(filePath);
                numberEventsLoadedSegment = length(segment.MCP_maxy);

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

        if length(channelsEnabled)==0 || any(ismember(channelsEnabled,padIDArray(i))) == 0
            channelsEnabled = [channelsEnabled;padIDArray(i)];
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

glbl_cut = time_diff_sigmoid>time_min & time_diff_sigmoid<time_max & ...
    MCP_maxy>min_accept*max(MCP_maxy) & MM_maxy>0.005 &...
    MM_maxy<max_accept*max(MM_maxy)& MCP_maxy<max_accept*max(MCP_maxy) & ...
    trackerX~=0 & trackerY~=0 & ~isnan(trackerX) & ~isnan(trackerY); %& acceptedSignalArray==1;


% glbl_cut = time_diff_sigmoid>time_min & time_diff_sigmoid<time_max & ...
%     trackerX~=0 & trackerY~=0;
%

%% loop through pads and determine timewalk correction per pad

timewalkParameters = cell(1,37);


for padNumber=1:37 %loop through pads

    channelMask = padIDArray==getChannelForPadNumberHighGranularity(padNumber);

    time_diff_sigmoidDUT = time_diff_sigmoid(channelMask);
    %     padNumber
    %     length(time_diff_sigmoidDUT)
    %     pause(3)

    MCP_maxyDUT = MCP_maxy(channelMask);
    MM_maxyDUT = MM_maxy(channelMask);
    trackerXDUT = trackerX(channelMask);
    trackerYDUT = trackerY(channelMask);
    eventIDArrayDUT = eventIDArray(channelMask);
    %MM_riseTimeDUT = MM_riseTime(channelMask);
    %MM_bgAvgDUT = MM_bgAvg(channelMask);
    %MM_bgRMSDUT = MM_bgRMS(channelMask);
    glbl_cutDUT = glbl_cut(channelMask);

    clear twalk
    twalk.p = [0 0 0];
    twalk.pRes = [0 0 0];

    if length(time_diff_sigmoidDUT)>10
        % Time walk analysis (this needs to be improoved)
        twalk.n_epk = 70;                               % number of e-peak bins
        e_peak_srt = sort(MM_maxyDUT(glbl_cutDUT));         % make temporary sort of e-peaks
        % try to distrubute e_peak vector evenly over the amplitude range
        twalk.epeak_vec = e_peak_srt(1:round(length(e_peak_srt)/twalk.n_epk):end); % 50 bins

        time_avg_raw = median(time_diff_sigmoidDUT);
        time_min = time_avg_raw - 0.3;     % predicted resolution 100ps cut 3 sigma
        time_max = time_avg_raw + 0.3;
        timeDiffMask = time_diff_sigmoidDUT>time_min & time_diff_sigmoidDUT<time_max;

        meanSAT = median(time_diff_sigmoidDUT(timeDiffMask));

        if length(twalk.epeak_vec)>0% && abs(mean(twalk.mean_sat))>0
            for i=1:length(twalk.epeak_vec)-1
                temp_cut = MM_maxyDUT > twalk.epeak_vec(i) & MM_maxyDUT < twalk.epeak_vec(i+1) & glbl_cutDUT;
                twalk.mean_sat(i) = mean(time_diff_sigmoidDUT(temp_cut));
                twalk.e_peak(i) =  mean(MM_maxyDUT(temp_cut));
                twalk.npts(i) = sum(temp_cut);
                twalk.rms(i) = std(time_diff_sigmoidDUT(temp_cut));
                twalk.err(i) = std(time_diff_sigmoidDUT(temp_cut))./sqrt(twalk.npts(i)); % mean error
                twalk.e_peak_err_p(i) = twalk.e_peak(i)-(twalk.epeak_vec(i+1)); % e charge limits
                twalk.e_peak_err_n(i) = -twalk.e_peak(i)+(twalk.epeak_vec(i));
            end

            % fit correction function using minuit
            %             fit_data = [];
            %             fit_data(1,:) = twalk.e_peak;
            %             fit_data(2,:) = twalk.mean_sat;
            %             fit_data(3,:) = twalk.err;
            %             p0=[];
            %             p0(1) = min(twalk.mean_sat);
            %             p0(2) = 1;
            %             p0(3) = 0.5;
            %             cmd='min; ret';
            %             [p, err, chi] = fminuit('twalk_fn_minuit',p0,fit_data,'-b','-c',cmd);
            %             twalk.p = p;
            %             twalk.chi = chi;

            % fit correction function using minuit
            %fit only points where error is < lim
            % Define the fixed parameter
            twkMask = twalk.err < 1.08;

            if abs(mean(twalk.mean_sat))>0
                startIndex = ceil(length(twalk.mean_sat(twkMask)) * 0.7);


                debugParam = twalk.mean_sat(twkMask)
                fixedParam = median(debugParam(startIndex:end));

                % Prepare fitting data
                fit_data = [];
                fit_data(1,:) = twalk.e_peak(twkMask);
                fit_data(2,:) = twalk.mean_sat(twkMask);
                fit_data(3,:) = twalk.err(twkMask);
                fit_data(4,:) = fixedParam;

                % Initial guesses for the variable parameters
                p0=[];
                p0 = [1, 0.5]; % Initial guesses for the parameters to be optimized

                % Define the fitting function with the fixed parameter
                % fittingFunc = @(p) twalk_fn_minuit_fix_offset(p, fit_data, fixedParam);

                % Perform the fit using fminuit
                cmd = 'min; ret';
                [p, err, chi] = fminuit('twalk_fn_minuit_fix_offset', p0,fit_data, '-b', '-c', cmd);
                p;
                % Store the results
                twalk.p = [fixedParam p];
                twalk.chi = chi;

                % Plotting
                figure;
                errorbar(twalk.e_peak(twkMask), twalk.mean_sat(twkMask), twalk.err(twkMask), ...
                    'o black', 'MarkerSize', 2, 'MarkerFaceColor', 'black', 'CapSize', 0);
                hold on;
                plot(twalk.e_peak, twalk_fn_minuit(twalk.p, twalk.e_peak), ...
                    'LineWidth', 1.3, 'Color', 'Red');
                xlabel('Electron peak amplitude (V)');
                ylabel('SAT, ns');


                message =  sprintf('  Fit parameters: %f, %f, %f - MeanSAT: %f',mean(twalk.mean_sat),twalk.p(2),twalk.p(3),meanSAT);
                %             y_pos=get(gcf,'ylim');
                %             x_pos=get(gcf,'xlim');
                %             text(x_pos(1),0.75*y_pos(2),message)


                title_str = sprintf('%s \n SAT vs. amplitude - Pad %d \n %s', runTitleString, padNumber, message);
                title(title_str)
                grid
                movegui(gcf, 'northwest');

                % Save the figure as a .png file
                pngFilePath = [store_folderPads '/timewalk-Run' run.id '_pad' num2str(padNumber) '.png'];
                saveas(gcf, pngFilePath);

                % Define the output directory and filename for the .txt file
                outputDir = [store_folder '/tWalkData/'];
                if ~exist(outputDir, 'dir')
                    mkdir(outputDir);
                end

                % Construct the full file path for the .txt file
                fullPath = [outputDir 'tWalkData_Pad' num2str(padNumber) '.txt'];

                % Open the file for writing
                fileID = fopen(fullPath, 'w');
                if fileID == -1
                    error('Failed to open file for writing: %s', fullPath);
                end

                % Loop through the data and write to the file
                for i = 1:length(twalk.e_peak)
                    fprintf(fileID, '%f %f %f %f %f\n', twalk.e_peak(i), twalk.mean_sat(i), twalk.err(i), twalk.e_peak_err_n(i), twalk.e_peak_err_p(i));
                end

                % Close the file
                fclose(fileID);
                fprintf('Data successfully saved to %s\n', fullPath);


                % fit correction function using minuit
                fit_dataRes = [];
                fit_dataRes(1,:) = twalk.e_peak;
                fit_dataRes(2,:) = twalk.rms;
                fit_dataRes(3,:) = twalk.err;
                p0=[];
                p0(1) = mean(twalk.rms);
                p0(2) = 1;
                p0(3) = 0.5;
                cmd='min; ret';
                [pRes, err, chi] = fminuit('twalk_fn_minuit',p0,fit_dataRes,'-b','-c',cmd);
                twalk.pRes = pRes;
                %twalk.chi = chi;
                % plot resolution vs. e-charge
                figure
                hold on
                errorbar(twalk.e_peak,twalk.rms*1000,[],[],twalk.e_peak_err_n,twalk.e_peak_err_p,'o black', 'MarkerSize', 2,'MarkerFaceColor', 'black', 'CapSize',0);
                plot(twalk.e_peak,1000*twalk_fn_minuit(pRes,twalk.e_peak),'LineWidth',1.3, 'Color', 'Red');
                xlabel('Electron peak amplitude (V)')
                ylabel('Resolution, ps')
                title_str = sprintf('%s \n Resolution vs. amplitude  - Pad %d',runTitleString,padNumber);
                title(title_str)
                ylim([0 200]);
                grid
                saveas(gcf,[store_folderPads '\twResolution-Run' run.id '_pad' num2str(padNumber) '.png'])
                % Define the output directory and filename
                outputDir = [store_folder '/tWalkData/'];
                if ~exist(outputDir, 'dir')
                    mkdir(outputDir);
                end

                % Construct the full file path
                fullPath = [outputDir '/ResolutionData_Pad' num2str(padNumber) '.txt'];

                % Open the file for writing
                fileID = fopen(fullPath, 'w');
                if fileID == -1
                    error('Failed to open file for writing: %s', fullPath);
                end

                % Loop through the data and write to the file
                for i = 1:length(twalk.e_peak)
                    fprintf(fileID, '%f %f %f %f\n', twalk.e_peak(i), twalk.rms(i)*1000, twalk.e_peak_err_n(i), twalk.e_peak_err_p(i));
                end

                % Close the file
                fclose(fileID);
                fprintf('Data successfully saved to %s\n', fullPath);
                close all;
            end
            % make time walk correction
            %             if(twalk.en == 1)
            %                 time_diff_sigmoidDUT = time_diff_sigmoidDUT - twalk_fn_minuit(twalk.p, MM_maxyDUT);
            %             end
        end
    end
    %write fit parameters to file for this channel
    dlmwrite([store_folderPads '\twCorrection_pad' num2str(padNumber) '.txt'],[twalk.p;twalk.pRes])
    twalk_p(padNumber) = twalk.p(1);
    timewalkParameters{padNumber} = [twalk.p;twalk.pRes];

end


eventsArray = [];
residualsArrayX = [];
residualsArrayY = [];
eventTrackerXArray = [];
eventTrackerYArray = [];
multiplicityArray = [];
sumAmpArray = [];
eventTimeArray = [];
check = 0;
%*
sumApms_plot = [];
notsummedAmps_plot = [];
notsummedSecondaryAmps_plot = [];
nsThidAmp_plot = [];
time_res = cell(1,37);
padAmplitudes = cell(1,37);
%*
%% loop through event IDs to assemble events in structures
for eventPos = 1:length(uniqueEventIDs)
    eventTemp.eventID = uniqueEventIDs(eventPos);

    shouldPlot = false;
    if mod(eventPos,plotInterval)==0
        shouldPlot = true;
    end

    %search for all events with this ID
    eventIDFilter = eventIDArray==eventTemp.eventID & glbl_cut; %& eventIDArray>10000;
    str_disp=sprintf('Plotting event ID %d', eventTemp.eventID);
    %disp(str_disp);

    chIDArrayEvent = padIDArray(eventIDFilter);

    if length(chIDArrayEvent)>0
        %valid event
        time_diff_sigmoidEvent = time_diff_sigmoid(eventIDFilter);
        MCP_maxyEvent = MCP_maxy(eventIDFilter);
        MM_maxyEvent = MM_maxy(eventIDFilter);



        trackerXEvent = trackerX(eventIDFilter);
        trackerYEvent = trackerY(eventIDFilter);
        MM_DataArrayEvent = MM_DataArray(eventIDFilter);

        %% invert X relative to center % may be needed for some runs - e.g
        %269, not needed for run 232

        trackerXEvent = -(trackerXEvent-pad0Pos(1))+pad0Pos(1);

        %        MCP_DataArrayEvent = MCP_DataArray(eventIDFilter);
        %acceptedSignalArrayEvent = acceptedSignalArray(eventIDFilter);
        numberAcceptedSignalsInEvent = length(MCP_maxyEvent)
        MM_maxySum = sum(MM_maxyEvent);
        sortedAmps = sort(MM_maxyEvent, 'descend');

        %% loop through events - select specific events here

        if max(MM_maxySum)>0.10 %&& numberAcceptedSignalsInEvent>2% && trackerXEvent(1) > padCenters(6,1) & trackerYEvent(1) < padCenters(10,2)% && length(chIDArrayEvent)==3 %show only events with strong signal
            size(MM_DataArrayEvent);
            waveformsArray = [];
            %*
            sumApms_plot = [sumApms_plot, MM_maxySum];
            notsummedAmps_plot = [notsummedAmps_plot, sortedAmps(1)];
            if numberAcceptedSignalsInEvent>1
                notsummedSecondaryAmps_plot = [notsummedSecondaryAmps_plot; sortedAmps(2)];
            else
                notsummedSecondaryAmps_plot = [notsummedSecondaryAmps_plot; 0];
            end

            if numberAcceptedSignalsInEvent>2
                nsThidAmp_plot = [nsThidAmp_plot; sortedAmps(3)];
            else
                nsThidAmp_plot = [nsThidAmp_plot; 0];
            end

            %*
            positionRecoX = 0;
            positionRecoY = 0;
            chIDArrayEvent = unique(chIDArrayEvent);
            % Loop through hit pads assemble waveforms array
            numChPresent = length(chIDArrayEvent);
            %padIDArrayEvent

            chIDArrayEvent;
            padIDArrayEvent = chIDArrayEvent;

            sumTimeDiffWeighted = 0;
            sumWeightingFactor = 0;
            eventTime = 0;

            % Loop through active pads in current event
            for pos = 1:length(chIDArrayEvent)
                ampWeight = MM_maxyEvent(pos)/MM_maxySum;
                chNumber = chIDArrayEvent(pos); %+1;

                padIndex = getPadNumber(chNumber);
                padIDArrayEvent (pos) = padIndex;

                if padIndex <=0
                    check = check + 1;
                    break
                end
                padIndex;
                posPadX = padCenters(padIndex,1);
                posPadY = padCenters(padIndex,2);

                % Add up estimate from weights
                positionRecoX = positionRecoX+posPadX*ampWeight;
                positionRecoY = positionRecoY+posPadY*ampWeight;

                waveformTemp.padID = padIDArrayEvent(pos);
                waveformTemp.x = MM_DataArrayEvent(pos).waveform(:,1) - MM_DataArrayEvent(pos).sigmoid.timepoint;
                waveformTemp.riseTime=MM_DataArrayEvent(pos).sigmoid.timepoint90-MM_DataArrayEvent(pos).sigmoid.timepoint10;% extract electrom peak charge

                waveformTemp.y = MM_DataArrayEvent(pos).waveform(:,2);
                %waveformTemp.xRef = MCP_DataArrayEvent(pos).sig.xfit;
                %waveformTemp.yRef = MCP_DataArrayEvent(pos).sig.yfit;
                waveformsArray = [waveformsArray waveformTemp];

                padTimewalkParameters = timewalkParameters{padIndex};

                %time_diff_sigmoidEventTWCorrected =  time_diff_sigmoidEvent(pos) - twalk_p(padIndex);
                time_diff_sigmoidEventTWCorrected =  time_diff_sigmoidEvent(pos)  - twalk_fn_minuit(padTimewalkParameters(1,:), MM_maxyEvent(pos));
                % meanTimeDiffAfterCorrection = mean(time_diff_sigmoidEventTWCorrected);
                padResolution =  twalk_fn_minuit(padTimewalkParameters(2,:), MM_maxyEvent(pos));

                sumTimeDiffWeighted = sumTimeDiffWeighted + time_diff_sigmoidEventTWCorrected / (padResolution^2);
                sumWeightingFactor = sumWeightingFactor + 1/(padResolution^2);

                %*
                time_res{padIndex}=[time_res{padIndex}, time_diff_sigmoidEventTWCorrected];
                padAmplitudes{padIndex} = [padAmplitudes{padIndex}, MM_maxyEvent(pos)];

                eventTime = eventTime + ampWeight*time_diff_sigmoidEventTWCorrected;
                %*
            end

            %determine combined time for event with info from all pads
            %eventTime = sumTimeDiffWeighted/sumWeightingFactor;

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
            eventTimeArray = [eventTimeArray;eventTime];

            if shouldDisplayMaps & saveMapCounter<numberMapsToSave
                VisualiseHighGranularityEventWaveformsPosReco(padCenters,padIDArrayEvent,1000*MM_maxyEvent,waveformsArray,[runTitleString ' - Event ' num2str(eventTemp.eventID)],[store_folderEvents '\Event' int2str(eventTemp.eventID) '.png'],'','%0.1f mV',minC,maxC,trackerXEvent(1),trackerYEvent(1),positionRecoX,positionRecoY)
                hold off;
                pause(1);
                close all;

                saveMapCounter = saveMapCounter+1;
            end
        end
    end
end

%% Plot
%plot multiplicity
figure;
hist(multiplicityArray,[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37]);
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

% Export data to .txt for later analysis
outputDir = [store_folder '/PositionResExtractedData'];
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end
fullPath = [outputDir '/PositionResX.txt'];
fileID = fopen(fullPath, 'w');
fprintf(fileID, '%f\n', residualsArrayX);
fclose(fileID);

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

% Export data to .txt for later analysis
outputDir = [store_folder '/PositionResExtractedData'];
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end
fullPath = [outputDir '/PositionResY.txt'];
fileID = fopen(fullPath, 'w');
fprintf(fileID, '%f\n', residualsArrayY);
fclose(fileID);


saveas(gcf,[store_folder '\Run' run.id '_residualsY.png'])
%%
figure;
histogram(sumAmpArray,100, 'DisplayStyle', 'stairs');
grid on
xlabel('Sum of pad signal amplitudes (mV)')
ylabel('Events');
%xlim([-10 10]);
title_str = sprintf('%s \n Sum of pad signal amplitudes    mean=%.2fV   std=%.2fV',runTitleString,mean(sumAmpArray),std(sumAmpArray));
title(title_str)
movegui(gcf,'west');

saveas(gcf,[store_folder '\Run' run.id '_sumAmp.png'])






%% extract residual in central area only
samplingRadius = 3;
samplingXCenter = pad0Pos(1)+3;
samplingYCenter = pad0Pos(2)+3;
cut_circSampling = ((eventTrackerXArray - samplingXCenter).^2 + (eventTrackerYArray - samplingYCenter).^2) < samplingRadius^2; %& glbl_cut;

%plot multiplicity
figure;
hist(multiplicityArray(cut_circSampling),[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37]);
grid on
xlabel('Multiplicity (pads hit)')
ylabel('Events');
title_str = sprintf('%s \n Pad Multiplicity  r=%.1fmm  mean=%.1f    std=%.1f',runTitleString,samplingRadius,mean(multiplicityArray(cut_circSampling)), std(multiplicityArray(cut_circSampling)));
title(title_str)
movegui(gcf,'south');

saveas(gcf,[store_folder '\Run' run.id '_multiplicity-SamplingR' int2str(samplingRadius) 'mm.png'])

figure;
hist(residualsArrayX,200);
grid on
xlabel('Residuals in X (mm)')
ylabel('Events');
%xlim([-10 10]);
title_str = sprintf('%s \n X Residuals  r=%.1fmm   mean=%.2fmm    std=%.2fmm',runTitleString,samplingRadius,mean(residualsArrayX(cut_circSampling)),std(residualsArrayX(cut_circSampling)));
title(title_str)
movegui(gcf,'north');

% Export data to .txt for later analysis
outputDir = [store_folder '/PositionResExtractedData'];
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end
fullPath = [outputDir '/PositionResX.txt'];
fileID = fopen(fullPath, 'w');
fprintf(fileID, '%f\n', residualsArrayX);
fclose(fileID);

saveas(gcf,[store_folder '\Run' run.id '_residualsX-SamplingR' int2str(samplingRadius) 'mm.png'])


figure;
hist(residualsArrayY,200);
grid on
xlabel('Residuals in Y (mm)')
ylabel('Events');
%xlim([-10 10]);
title_str = sprintf('%s \n Y Residuals r=%.1fmm   mean=%.2fmm   std=%.2fmm',runTitleString,samplingRadius,mean(residualsArrayY(cut_circSampling)),std(residualsArrayY(cut_circSampling)));
title(title_str)
movegui(gcf,'east');

% Export data to .txt for later analysis
outputDir = [store_folder '/PositionResExtractedData'];
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end
fullPath = [outputDir '/PositionResY.txt'];
fileID = fopen(fullPath, 'w');
fprintf(fileID, '%f\n', residualsArrayY);
fclose(fileID);


saveas(gcf,[store_folder '\Run' run.id '_residualsY-SamplingR' int2str(samplingRadius) 'mm.png'])










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
area.step = 1;   % set grid resolution in mm
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
        area.eventTimeCombined(i,j) = rms(eventTimeArray(cut_circ));
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
for padNumber=1:37 %loop through pads
    padNumber = padNumber;
    padCenter = padCenters(padNumber,:);
    padX = padCenter(1);
    padY = padCenter(2);
    hexagonOutlineRotated(padSide,padX,padY);
end
caxis([0 3]);
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
for padNumber=1:37 %loop through pads
    padNumber = padNumber;
    padCenter = padCenters(padNumber,:);
    padX = padCenter(1);
    padY = padCenter(2);
    hexagonOutlineRotated(padSide,padX,padY);
end
caxis([0 2]);
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
for padNumber=1:37 %loop through pads
    padNumber = padNumber;
    padCenter = padCenters(padNumber,:);
    padX = padCenter(1);
    padY = padCenter(2);
    hexagonOutlineRotated(padSide,padX,padY);
end
axis equal
caxis([0 2]);
saveas(gcf,[store_folder '\Run' run.id '_ResidualYMean-Map.png'])
movegui(gcf,'northeast');
pause(1);



%% plot combined event Time res as map
figure
h=pcolor(area.xx,area.yy,1000*area.eventTimeCombined);
hold on
set(h, 'EdgeColor', 'none');
axis equal
colorbar
xlabel('x-axis, mm');
ylabel('y-axis, mm');
h = colorbar;
h.Label.String = 'RMS Combined Event Time (ps)';
h.Label.Position(1) = 3;
str_title = sprintf('%s: RMS Combined Event Time', runTitleString);
title(str_title);
for padNumber=1:37 %loop through pads
    padNumber = padNumber;
    padCenter = padCenters(padNumber,:);
    padX = padCenter(1);
    padY = padCenter(2);
    hexagonOutlineRotated(padSide,padX,padY);
end
axis equal
caxis([0 100]);
saveas(gcf,[store_folder '\Run' run.id '_RMSCombinedEventTime-Map.png'])
movegui(gcf,'northeast');
pause(1);

%% Plot time resolution per pad

mus = [];
sigmas = [];
rmss = [];
for i = 1:37
    padNumber = i;
    time_avg_raw = median(time_res{i});
    time_min = time_avg_raw - 0.3;     % predicted resolution 100ps cut 3 sigma
    time_max = time_avg_raw + 0.3;
    filtered_data = time_res{i}(time_res{i} >= time_min & time_res{i} <= time_max);

    % Export data to .txt for later analysis
    outputDir = [store_folder '/TimeResExtractedData'];
    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end
    fullPath = [outputDir '/TimeResPad' num2str(padNumber) '.txt'];
    fileID = fopen(fullPath, 'w');
    fprintf(fileID, '%f\n', filtered_data);
    fclose(fileID);

    if length(filtered_data) ~= 0
        pd = fitdist(filtered_data', 'Normal');
        rmss = [rmss, rms(filtered_data)];
        mus = [mus, pd.mu];
        sigmas = [sigmas, pd.sigma];
        pause (1);
        % Determine optimal number of bins (e.g., using Freedman-Diaconis rule)
        nbins = ceil((max(filtered_data) - min(filtered_data)) / (2 * iqr(filtered_data) * length(filtered_data)^(-1/3)));

        figure;
        histfit(filtered_data,nbins,"normal")

        message =  sprintf('  RMS_{tot} = %2.1f ps ',1000*std(filtered_data));

        timeMisalignment = mean(filtered_data);

        y_pos=get(gca,'ylim');
        x_pos=get(gca,'xlim');
        text(x_pos(1),0.75*y_pos(2),message)
        xlabel('Time difference (ns)')
        ylabel('Events')
        title_str = sprintf('%s \n Time diff hist  - Pad %d - mean: %f',runTitleString,padNumber,timeMisalignment);
        title(title_str)
        %ylim([0 200]);
        grid
        saveas(gcf,[store_folderPads '\timeDiffHist-Run' run.id '_pad' num2str(padNumber) '.png'])
        close all;
    else
        mus = [mus, 0];
        sigmas = [sigmas, 0];
        rmss = [rmss, 0];
    end

end
%channelsEnabled = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19]
%channelsEnabled = [1,2,3,4,5,6,7,8,9,10,11,12,13,101,102,105,107,108,109]
%channelsEnabled = [107,102,108,105,101,13,109,1,5,10,9,8,7,12,6,11,4,3,2];
%channelsEnabled = [65;66;67;68;69;70;71;72;73;74;75;76;77;78;79;80;81;82;83;84;85;86;87;88];

VisualiseHighGranularity(channelsEnabled,rmss,'Map Multipad Picosec Time Resolution',[store_folder '\Run' run.id '_TimeResMAP.png'],'Res (ns)','%0.5f %',0.03,0.06)

%% Time resolution combined event
time_avg_raw = median(eventTimeArray)
time_min = time_avg_raw - 0.3;     % predicted resolution 100ps cut 3 sigma
time_max = time_avg_raw + 0.3;
filteredEventTimeArrayData = eventTimeArray(eventTimeArray > time_min & eventTimeArray < time_max);
% Export data to .txt for later analysis
outputDir = [store_folder '/TimeResExtractedDataEvent'];
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end
fullPath = [outputDir '/TimeResEvents.txt'];
fileID = fopen(fullPath, 'w');
fprintf(fileID, '%f\n', filteredEventTimeArrayData);
fclose(fileID);

nbins = ceil((max(filteredEventTimeArrayData) - min(filteredEventTimeArrayData)) / (2 * iqr(filteredEventTimeArrayData) * length(filteredEventTimeArrayData)^(-1/3)));

figure;
histfit(filteredEventTimeArrayData,nbins,"normal")

message =  sprintf('  RMS_{tot} = %2.1f ps ',1000*std(filteredEventTimeArrayData));

timeMisalignment = mean(filteredEventTimeArrayData);

y_pos=get(gca,'ylim');
x_pos=get(gca,'xlim');
text(x_pos(1),0.75*y_pos(2),message)
xlabel('Time difference (ns)')
ylabel('Events')
title_str = sprintf('%s \n Time diff hist  - Events', runTitleString);
title(title_str)
%ylim([0 200]);
grid
saveas(gcf,[outputDir '\timeDiffHist-Run' run.id '_Event.png'])
close all

%% Plot amplitude map for individual pads

musAmplitudes = [];
sigmasAmplitudes = [];
rmssAmplitudes = [];

P = ["N","theta","nBar"];
fitfun = fittype(@(N,theta,nBar,x) (N./nBar).*((theta+1).^(theta+1)).*((x./nBar).^theta).*exp(-(theta+1).*x./nBar)./gamma(theta+1));
options = fitoptions(fitfun);
options.Upper = [1000 150 0.1];
options.Lower = [0 0 0];
options.StartPoint = [1, 0.01, 0.04]; % Adjust the starting points if needed

for i = 1:37
    filtered_data = padAmplitudes{i};
    % Export data to .txt for later analysis
    outputDir = [store_folder '/AmpExtractedData'];
    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end
    fullPath = [outputDir '/AmpPad' num2str(i) '.txt'];
    fileID = fopen(fullPath, 'w');
    fprintf(fileID, '%f\n', filtered_data);
    fclose(fileID);

    if ~isempty(filtered_data)
        pause(1);
        figure;
        h = histogram(filtered_data, 50);

        % Filter out empty bins for fitting
        bin_centers = h.BinEdges(1:end-1) + h.BinWidth / 2;
        non_zero_bins = h.Values > 0;

        try
            % Fit the curve using only non-zero bins
            [fitted_curve, gof] = fit(bin_centers(non_zero_bins)', h.Values(non_zero_bins)', fitfun, options);

            % Store the fitted parameters
            P = [P; fitted_curve.N, fitted_curve.theta, fitted_curve.nBar];

            % Plot the fitted curve
            if fitted_curve.theta > 0
                hold on;
                plot(bin_centers(non_zero_bins), fitted_curve(bin_centers(non_zero_bins)), 'LineWidth', 3, 'Color', [0, 0, 1]);

                % Calculate chi-squared and degrees of freedom
                ch2 = sum(((h.Values(non_zero_bins) - fitted_curve(bin_centers(non_zero_bins))').^2) ./ fitted_curve(bin_centers(non_zero_bins))');
                dof = sum(non_zero_bins) - 3;
                nch2 = ch2 / dof;
                np = 1 - chi2cdf(ch2, dof); % P(chi^2 > ch2)

                rmssAmplitudes = [rmssAmplitudes, rms(filtered_data)];
                musAmplitudes = [musAmplitudes, fitted_curve.nBar];

            end

            % Display results
            message = sprintf('amplitude = %2.1f mV', 1000*fitted_curve.nBar);
            y_pos = get(gca, 'ylim');
            x_pos = get(gca, 'xlim');
            text(x_pos(1), 0.75 * y_pos(2), message);
            xlabel('Amp (V)');
            ylabel('Events');
            title_str = sprintf('%s \n Amp hist  - Pad %d', runTitleString, i);
            title(title_str);
            grid on;

            % Save the figure
            saveas(gcf, [store_folderPads '\AmpHist-Run' run.id '_pad' num2str(i) '.png']);
            close all;
        catch exception
            disp(getReport(exception));
            musAmplitudes = [musAmplitudes, 0];
            sigmasAmplitudes = [sigmasAmplitudes, 0];
            rmssAmplitudes = [rmssAmplitudes, 0];

            continue;
        end
    else
        musAmplitudes = [musAmplitudes, 0];
        sigmasAmplitudes = [sigmasAmplitudes, 0];
        rmssAmplitudes = [rmssAmplitudes, 0];
    end
end

VisualiseHighGranularity(channelsEnabled, musAmplitudes, 'Map Multipad Picosec Mean Amplitude', [store_folder '\Run' run.id '_AmpMAP.png'], 'Amp (V)', '%0.5f %', 0.03, 0.06);

%%
figure;

% Subplot 1: Main Scatter Plot with Fitted Lines
subplot(2,1,1); % 2 rows, 1 column, 1st subplot
hold on;

% Plot the first scatter plot
scatter(sumApms_plot, notsummedAmps_plot, '.b'); % Blue dots for the first group

% Fit a line to the first data set
coeff1 = polyfit(sumApms_plot, notsummedAmps_plot, 1);
fittedLine1 = polyval(coeff1, sumApms_plot);
plot(sumApms_plot, fittedLine1, '-b', 'LineWidth', 1.5); % Plot the fitted line

% Calculate residuals for the first group
residuals1 = notsummedAmps_plot - fittedLine1;
RMS1 = rms(residuals1); % Root mean square of the residuals

% Plot the second scatter plot
scatter(sumApms_plot, notsummedSecondaryAmps_plot, '.r'); % Red dots for the second group

% Fit a line to the second data set
coeff2 = polyfit(sumApms_plot, notsummedSecondaryAmps_plot', 1);
fittedLine2 = polyval(coeff2, sumApms_plot);
plot(sumApms_plot, fittedLine2, '-r', 'LineWidth', 1.5); % Plot the fitted line

% Calculate residuals for the second group
residuals2 = notsummedSecondaryAmps_plot' - fittedLine2;
RMS2 = rms(residuals2'); % Root mean square of the residuals

% Plot the third scatter plot
scatter(sumApms_plot, nsThidAmp_plot, '.g'); % Green dots for the third group

% Fit a line to the third data set
coeff3 = polyfit(sumApms_plot, nsThidAmp_plot', 1);
fittedLine3 = polyval(coeff3, sumApms_plot);
plot(sumApms_plot, fittedLine3, '-g', 'LineWidth', 1.5); % Plot the fitted line

% Calculate residuals for the third group
residuals3 = nsThidAmp_plot' - fittedLine3;
RMS3 = rms(residuals3); % Root mean square of the residuals

% Add the legend with gradient (slope) values and RMS of residuals
legend('Leading pad', ...
    sprintf('Slope = %.4f, RMS = %.4f', coeff1(1), RMS1), ...
    'Secondary pad', ...
    sprintf('Slope = %.4f, RMS = %.4f', coeff2(1), RMS2), ...
    'Third pad', ...
    sprintf('Slope = %.4f, RMS = %.4f', coeff3(1), RMS3), ...
    'Location', 'best');

title('Scatter Plot with Fitted Lines');
hold off;

% Subplot 2: Residuals Plot
subplot(2,1,2); % 2 rows, 1 column, 2nd subplot
hold on;

% Plot residuals for all groups
scatter(sumApms_plot, residuals1, '.b'); % Blue residuals for the first group
scatter(sumApms_plot, residuals2', '.r'); % Red residuals for the second group
scatter(sumApms_plot, residuals3, '.g'); % Green residuals for the third group

% Add a horizontal line at y=0 for reference
yline(0, '-k', 'LineWidth', 0.5);

% Add labels and title to the residuals plot
xlabel('Sum Apms');
ylabel('Residuals');
title('Residuals of the Fitted Lines');

hold off;

% Link the x-axes of both subplots
h = findall(gcf, 'Type', 'axes');
linkaxes(h, 'x');

%%


function padID = getPadNumber(chID)
setMappingHighGranularity

padID = -1;

for k=1:length(mapping)
    if mapping(k,2) == chID
        padID = mapping(k,1);
        break;
    end
end

end
