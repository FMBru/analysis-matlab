%List the following in batchAnalyseList.txt.txt for analysing in batch
%RunNumber,Oscilloscope,CH of DUT, CH of REF MCP, DetectorGeometry,
%FilesToAnalyse (0=all)
% DUT index 
%1: miniCactus, 2: MM3 (VacChamber), 3: MM1(Multipad), 4: MM2 (on support plate), 5: MM4 (electron setup)
%e.g. 112,Pool2,C4,C1,2,0
clear all

addpath 'C:\Users\GDD\Documents\MATLAB\Picosec\CommonFunctions';
addpath 'C:\Users\GDD\Documents\MATLAB\Picosec\CommonFunctions\SCP';

fid = fopen('batchAnalyseList_PicolargeMultiplePads.txt');
% Read all lines & collect in cell array
analyseEntries = textscan(fid,'%s','delimiter','\n');
fclose(fid);

%loop through entries in batchAnalyseList
for  i=2:length(analyseEntries{1})
    i;
    %skipping first line
    close all;
    
    %% config
    shouldRemoveFilesAfterProcessing = false;
    shouldDownloadFiles = false;
    useEOSFolder = true;
    
    %% processing
%     fid = fopen('batchAnalyseList.txt');
%     % Read all lines & collect in cell array
%     analyseEntries = textscan(fid,'%s','delimiter','\n');
%     fclose(fid);
    %length(analyseEntries{1})
    
    firstStringCellArray =  analyseEntries{1}(i);
    firstString = char(firstStringCellArray);
    stringArray = split(firstString,',');
    
    batchProcess = 1;
    
    %parse line
    runIDString = stringArray{1};
    padIDString = stringArray{2};
    oscilloscopeString = stringArray{3};
    
    DUTnameString = stringArray{4};
    MCPnameString = stringArray{5};
    trackerDUTIndex = '2';
    noFileToAnalyze = str2num(stringArray{6});
    
    DUTnameString = erase(DUTnameString,"C");
    MCPnameString = erase(MCPnameString,"C");

    str_disp=sprintf('Starting to process Run %s - %s DUT:%d REF:%d for Pad: %d', runIDString, oscilloscopeString, str2num(DUTnameString), str2num(MCPnameString), str2num(padIDString));
        disp(str_disp);
    pause(1);
    
    
    ProcessRawFiles_PicolargeMultiplePads

    %analyse data, plotting, upload results
    %AnalyseRun
    
    pause(1);
    
end

%processed all files
%load variables

clear MM_data MCP_data time_diff time_diff_sigmoid MCP_maxy MM_maxy trackerX trackerY eventIDArray acceptedSignalArray;

loadVariablesFromFiles_PicolargeMultiplePads

