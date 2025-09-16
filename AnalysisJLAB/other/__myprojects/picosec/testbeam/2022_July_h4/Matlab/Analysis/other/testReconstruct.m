 

user = 'gdd';
password = 'Win_Admin'; %remove your password when done

oscilloscope = "Pool2";
runNumber = '123';

addpath 'C:\Users\GDD\Documents\MATLAB\Picosec\CommonFunctions';
addpath 'C:\Users\GDD\Documents\MATLAB\Picosec\CommonFunctions\SCP';
connObj  =  ssh2_config('lxplus.cern.ch',user,password,22);

dataFilePathCommand = append('cd ');
   % [connObj, command_result] = ssh2_command(connObj, convertStringsToChars(dataFilePathCommand), 1);
    [connObj, command_result] = ssh2_command(connObj, convertStringsToChars("source /eos/project/p/picosec/analysis/anamicom/lmfit_lib.source"), 1);
    [connObj, command_result] = ssh2_command(connObj, convertStringsToChars("/eos/project/p/picosec/analysis/anamicom/anastrip 164"), 1);
        
        
ssh2_close(connObj);
