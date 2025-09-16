user = 'fbrunbau';
password = 'XX!';

oscilloscope = "Pool1";
runNumber = '391';

%download reconstructed tracker file
trackerFileName=['asciiRun' runNumber '.dat'];
%ssh2_struct = scp_simple_get('lxplus.cern.ch', user, password, trackerFileName, 'C:\Users\GDD\Documents\Picosec\May22\Tracker\reconstructed', '/eos/project/p/picosec/testbeam/2022_May_h4/tracker/reconstructed/')

%download data file
%make connection
connObj  =  ssh2_config('lxplus.cern.ch',user,password,22);

%list files in folder
dataFilePath = append('/eos/project/p/picosec/testbeam/2022_May_h4/',oscilloscope,'/Run',runNumber);
dataFilePathCommand = append('ls ',dataFilePath);
[connObj, command_result] = ssh2_command(connObj, convertStringsToChars(dataFilePathCommand), 1);

command_result

%download all fles in folder
dataFolderLocal=append('C:\Users\GDD\Documents\Picosec\May22\',oscilloscope,'\Run',runNumber);

%ssh2_struct = scp_simple_get('lxplus.cern.ch', user, password, 'SealedModule_Collimatortest.gcode', convertStringsToChars(dataFolderLocal), convertStringsToChars(dataFilePath))
%ssh2_struct = scp_simple_get('lxplus.cern.ch', user, password, 'SealedModule_Collimatortest.gcode', 'C:\Users\GDD\Documents\Picosec\May22\Pool1\Run391', '/eos/project/p/picosec/testbeam/2022_May_h4/Pool1/Run391')


%ssh2_struct = scp_get(connObj, 'SealedModule_Collimatortest.gcode', convertStringsToChars(dataFolderLocal), convertStringsToChars(dataFilePath))
ssh2_struct = scp_get(connObj, command_result, convertStringsToChars(dataFolderLocal), convertStringsToChars(dataFilePath))

ssh2_close(connObj);
