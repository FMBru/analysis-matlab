
runIDsToAnalyse = ["339"];

%read in and save all signals and signalsRef to variables
%ProcessSAMPICData_Multipad

%read in variables, process them and save matched signals
%processSignalsFromSavedFilesSegmented_Multipad

%load analysed data and plot
loadVariablesFromFiles_Multipad

%create plots for all pads
AnalyseAllChannels