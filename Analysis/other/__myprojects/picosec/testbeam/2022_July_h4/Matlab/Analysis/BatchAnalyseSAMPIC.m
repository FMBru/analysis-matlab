%List the following in batchAnalyseList.txt.txt for analysing in batch
%RunNumber,Oscilloscope,CH of DUT, CH of REF MCP, DetectorGeometry
% DUT index 
%1: miniCactus, 2: MM3 (VacChamber), 3: MM1(Multipad), 4: MM2 (on support plate), 5: MM4 (electron setup)
%e.g. 112,Pool2,C4,C1,2

addpath 'C:\Users\GDD\Documents\MATLAB\Picosec\CommonFunctions';
addpath 'C:\Users\GDD\Documents\MATLAB\Picosec\CommonFunctions\SCP';
addpath 'C:\Users\GDD\Documents\MATLAB\Picosec\SAMPIC';

fid = fopen('batchAnalyseList.txt');
% Read all lines & collect in cell array
analyseEntries = textscan(fid,'%s','delimiter','\n');
fclose(fid);

%loop through entries in batchAnalyseList
%while  length(analyseEntries{1})>0
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
    
    firstStringCellArray =  analyseEntries{1}(1);
    firstString = char(firstStringCellArray);
    stringArray = split(firstString,',');
    
    batchProcess = 1;
    
    %parse line
    runIDString = stringArray{1}
    oscilloscopeString = stringArray{2}
    
    DUTnameString = stringArray{3}
    MCPnameString = stringArray{4}
    trackerDUTIndex = stringArray{5}
    
    %download files
    if shouldDownloadFiles
        if useEOSFolder
            shouldRemoveFilesAfterProcessing = false;
        else
            downloadRun
        end
    end
    
    pause(1);
    
    %read in data
    disp('Read in SAMPIC data file');
    readSAMPICData
    
    %match events to each other based on time differences
        disp('Match hits to events');

    matchSAMPICSignals
    
    %process files, save analysed .mat
            disp('Processing files...');

    ProcessRawFilesSAMPIC
    
    pause(1);
    
    
            disp('Finished processing');
    %analyse data, plotting, upload results
    AnalyseRun
    
    pause(1);
    
    %delete data files
    dataFolder=convertStringsToChars(append('C:\Users\GDD\Documents\Picosec\May22\',oscilloscopeString,'\Run',runIDString));
    if shouldRemoveFilesAfterProcessing
        rmdir(dataFolder, 's')
    end
    
    pause(1);
    
    fid = fopen('batchAnalyseList.txt', 'r') ;              % Open source file.
    fgetl(fid) ;                                  % Read/discard line.
    buffer = fread(fid, Inf) ;                    % Read rest of the file.
    fclose(fid)
    
    fid = fopen('batchAnalyseList.txt', 'w')  ;   % Open destination file.
    fwrite(fid, buffer) ;                         % Save to file.
    fclose(fid) ;
    
    %look if furher analysis entries are in batch file
    fid = fopen('batchAnalyseList.txt');
    analyseEntries = textscan(fid,'%s','delimiter','\n');
    fclose(fid);
%end