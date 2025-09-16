user = 'gdd';
password = 'Win_Admin'; %remove your password when done

oscilloscope = "Pool2";
runNumber = '123';

addpath 'C:\Users\GDD\Documents\MATLAB\Picosec\CommonFunctions';
addpath 'C:\Users\GDD\Documents\MATLAB\Picosec\CommonFunctions\SCP';

%addpath 'C:\Users\GDD\Documents\MATLAB\Picosec\SCP'
connObj  =  ssh2_config('lxplus.cern.ch',user,password,22);

    if exist('batchProcess','var') == 1
        runNumber = runIDString;
        oscilloscope = oscilloscopeString;
        shouldDownloadTrackerFile = true;
        shouldDownloadDataFiles = false;
    else 
        %not batch processing
        run.id = '131';
        run.oscilloscope = 'Pool3';
        shouldDownloadTrackerFile = true;
        shouldDownloadDataFiles = true;
    end

    



if shouldDownloadTrackerFile
    %download reconstructed tracker file
    trackerFileName=['asciiRun' runNumber '.dat'];
    dataFolderLocal=append('C:\Users\GDD\Documents\Picosec\May22\Tracker\reconstructed');
    dataFilePath = append('/eos/project/p/picosec/testbeam/2022_May_h4/tracker/reconstructed');
    ssh2_struct = scp_get(connObj, trackerFileName, convertStringsToChars(dataFolderLocal), convertStringsToChars(dataFilePath))
end


%download data file
%make connection
if shouldDownloadDataFiles
    %list files in folder
    dataFilePath = append('/eos/project/p/picosec/testbeam/2022_May_h4/',oscilloscope,'/Run',runNumber);
    dataFilePathCommand = append('ls ',dataFilePath);
    [connObj, command_result] = ssh2_command(connObj, convertStringsToChars(dataFilePathCommand), 1);
    
    %download all fles in folder
    dataFolderLocal=append('C:\Users\GDD\Documents\Picosec\May22\',oscilloscope,'\Run',runNumber);
    mkdir(dataFolderLocal);
    
    ssh2_struct = scp_get(connObj, command_result, convertStringsToChars(dataFolderLocal), convertStringsToChars(dataFilePath))
    
end
ssh2_close(connObj);
