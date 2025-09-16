
runIDsToAnalyse = ["311"];
channelToAnalyse = [0 12 13 14 22 23 24 22 33 34];

% runIDsToAnalyse = ["339"];
% channelToAnalyse = [0 23 24 33 34];


%read in and save all signals and signalsRef to variables
%ProcessSAMPICData_Multipad_Light

%read in variables, process them and save matched signals
%processSignalsFromSavedFilesSegmented_Multipad_Light

%load analysed data and plot
loadVariablesFromFiles_Multipad

alignFromTracker

%create plots for all pads
AnalyseAllChannels