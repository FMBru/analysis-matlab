%List the following in batchAnalyseList.txt for analysing in batch
%RunNumber,Oscilloscope,CH of DUT, CH of REF MCP, DetectorGeometry,
%FilesToAnalyse (0=all)
% DUT index 
%DUT index: 1:MM1, 2:MM2, 3:MM3, 4:MM4, 5:MM5, 6:MM6
%e.g. 112,Pool2,C4,C1,2,0
%pause(300)

addpath '\\eosproject-smb\eos\project\p\picosec\testbeam\2025_July_h4\Matlab\CommonFunctions';
addpath '\\eosproject-smb\eos\project\p\picosec\testbeam\2025_July_h4\Matlab\Analysis\Functions';

filename = 'batchAnalyseList.txt';
fid = fopen(filename);
% Read all lines & collect in cell array
analyseEntries = textscan(fid,'%s','delimiter','\n');
fclose(fid);
        
        %loop through entries in batchAnalyseList
while  length(analyseEntries{1})>1

    %skipping first line
    close all;
    clear all;
    filename = 'batchAnalyseList.txt';

    %% config
    shouldRemoveFilesAfterProcessing = false;
    shouldDownloadFiles = false;
    useEOSFolder = true;


    %% processing
    fid = fopen(filename);
    % Read all lines & collect in cell array
    analyseEntries = textscan(fid,'%s','delimiter','\n');
    fclose(fid);
    %length(analyseEntries{1})

    firstStringCellArray =  analyseEntries{1}(2);
    firstString = char(firstStringCellArray);
    stringArray = split(firstString,',');

    batchProcess = 1;

    %parse line
    runIDString = stringArray{1};
    oscilloscopeString = stringArray{2};

    DUTnameString = stringArray{3};
    MCPnameString = stringArray{4};
    trackerDUTIndex = stringArray{5};
    noFileToAnalyze = str2num(stringArray{6});
    runInfoDesc = stringArray{7};
    enableNoiseRej = str2num(stringArray{8});
    enableFilter = str2num(stringArray{9});
    enableIonTail = str2num(stringArray{10});

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
    ProcessRawFilesFrancesco
    end

    pause(1);

    %AnalyseRun %already in ProcessRawFiles

    pause(1);


    %delete data files
%     dataFolder=convertStringsToChars(append('C:\Users\GDD\Documents\Picosec\May22\',oscilloscopeString,'\Run',runIDString));
%     if shouldRemoveFilesAfterProcessing
%         rmdir(dataFolder, 's');
%     end

    pause(1);

    fid = fopen(filename, 'r') ;              % Open source file.
    line1 = fgetl(fid) ;                                  % Read/discard line.
    lineDiscard = fgetl(fid) ;                                  % Read/discard line.
    buffer = fread(fid, Inf) ;                    % Read rest of the file.
    fclose(fid);

    fid = fopen(filename, 'w')  ;   % Open destination file.
    fwrite(fid, line1) ;                         % Save to file.
    fwrite(fid, buffer) ;                         % Save to file.
    fclose(fid) ;

    %look if furher analysis entries are in batch file
    fid = fopen(filename);
    analyseEntries = textscan(fid,'%s','delimiter','\n');
    fclose(fid);
end
