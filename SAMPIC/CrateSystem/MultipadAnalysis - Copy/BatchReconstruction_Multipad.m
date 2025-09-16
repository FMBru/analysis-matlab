
runIDsToAnalyse = ["339"];
channelToAnalyse = [0 13 14 23 24 33 34];


%read in and save all signals and signalsRef to variables
ProcessSAMPICData_Multipad_Light

%read in variables, process them and save matched signals
processSignalsFromSavedFilesSegmented_Multipad_Light

%load analysed data and plot
loadVariablesFromFiles_Multipad

alignFromTracker

%create plots for all pads
AnalyseAllChannels