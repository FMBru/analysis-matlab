runIDsToAnalyse = ["206","218","276"];
runIDsToAnalyse = ["206"];


%read in and save all signals and signalsRef to variables
ProcessSAMPICData_MediumGranularityFull

%read in variables, process them and save matched signals
processSignalsFromSavedFilesSegmented

%load analysed data and plot
loadVariablesFromFiles_MediumGranularityMultiplePads_Light

%create plots for all pads
AnalyseAllChannels