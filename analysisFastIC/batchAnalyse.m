%List the following in batchAnalyseList.txt.txt for analysing in batch
%RunNumber,Oscilloscope,CH of DUT, CH of REF MCP, DetectorGeometry,
%FilesToAnalyse (0=all)
% DUT index
%DUT index: 1:MM1, 2:MM2, 3:MM3, 4:MM4, 5:MM5, 6:MM6
%e.g. 112,Pool2,C4,C1,2,0
% x
% x
% x
% 189,Pool3,4,1,1,0,FastIC-440C-275A-multi-channel-ch0
% x
% 189,Pool3,2,1,1,0,FastIC-440C-275A-multi-channel-ch1
% x
% 189,Pool2,4,1,1,0,FastIC-440C-275A-multi-channel-ch4
% x
% 189,Pool2,2,1,1,0,FastIC-440C-275A-multi-channel-ch5

% x

addpath '.\commonFunctions';
addpath '.\commonFunctions\scp';
addpath '.\functions';

fid = fopen('batchAnalyseList.txt');
% Read all lines & collect in cell array
analyseEntries = textscan(fid,'%s','delimiter','\n');
fclose(fid);

%loop through entries in batchAnalyseList
while  length(analyseEntries{1})>1

    %skipping first line
    close all;
    clear all;

    %% config
    shouldRemoveFilesAfterProcessing = false;
    shouldDownloadFiles = false;
    useEOSFolder = true;


    %% processing
    fid = fopen('batchAnalyseList.txt');
    % Read all lines & collect in cell array
    analyseEntries = textscan(fid,'%s','delimiter','\n');
    fclose(fid);
    %length(analyseEntries{1})

    firstStringCellArray =  analyseEntries{1}(2);
    firstString = char(firstStringCellArray);
    stringArray = split(firstString,',');

    batchProcess = 1;

    %parse line
    runIDString = stringArray{1}
    oscilloscopeString = stringArray{2}

    DUTnameString = stringArray{3};
    MCPnameString = stringArray{4};
    trackerDUTIndex = stringArray{5};
    noFileToAnalyze = str2num(stringArray{6});
    runInfoDesc = stringArray{7};

    DUTnameString = erase(DUTnameString,"C");
    MCPnameString = erase(MCPnameString,"C");

    %download files
    if shouldDownloadFiles
        if useEOSFolder
            shouldRemoveFilesAfterProcessing = false;
        else
            downloadRun
        end
    end

    pause(1);

    %process files, save analysed .mat
    if strcmp(oscilloscopeString,'SAMPIC')
        %is sampic, process fiels by reading in from SAMPIC ASCII, binary
        %process sampic events
        ProcessSAMPICData
    else
        %is oscilloscope, read waveforms directly
        processRawFiles
    end

    pause(1);

    %analyse data, plotting, upload results
    analyseRun

    pause(1);


    %delete data files
    dataFolder=convertStringsToChars(append('C:\Users\GDD\Documents\Picosec\May22\',oscilloscopeString,'\Run',runIDString));
    if shouldRemoveFilesAfterProcessing
        rmdir(dataFolder, 's')
    end

    pause(1);

    fid = fopen('batchAnalyseList.txt', 'r') ;              % Open source file.
    line1 = fgetl(fid) ;                                  % Read/discard line.
    lineDiscard = fgetl(fid) ;                                  % Read/discard line.
    buffer = fread(fid, Inf) ;                    % Read rest of the file.
    fclose(fid)

    fid = fopen('batchAnalyseList.txt', 'w')  ;   % Open destination file.
    fwrite(fid, line1) ;                         % Save to file.
    fwrite(fid, buffer) ;                         % Save to file.
    fclose(fid) ;

    %look if furher analysis entries are in batch file
    fid = fopen('batchAnalyseList.txt');
    analyseEntries = textscan(fid,'%s','delimiter','\n');
    fclose(fid);

end