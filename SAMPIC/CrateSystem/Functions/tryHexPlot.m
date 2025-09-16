
clear all
colorMap = colormap;

channelsEnabled = [1,2,3,4,5,6,7];
                            colorMapEntry = colorMap(50,:)
                            

close all
dataArrayInput=[3,2,3,4,5,1,6,7];
plotPath = "text.png";
plottingString = "%.1f mV"

VisualisePicolarge (channelsEnabled,dataArrayInput,"Uniformity map",plotPath,"Amp",plottingString,0,0)


