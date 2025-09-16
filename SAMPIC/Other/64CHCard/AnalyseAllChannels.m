padResultsArray = [];
dutChannelArrayToAnalyse = channelToAnalyse;
%dutChannelArray = [50;40;30;20]; %override which channels to analyse
shouldSaveMat=false;

for chPos = 1:length(dutChannelArrayToAnalyse)
    analysis.dutChannel =  dutChannelArrayToAnalyse(chPos);
    
    str_disp=sprintf('Analysing channel %d', chPos);
    disp(str_disp);
    
    pad.yc = 0;
    pad.xc = 0;
    padRMSSampling = 0;
    padSATMeanSampling = 0;
    padMeanMCPAmpSampling = 0;
    padMeanMMAmpSampling = 0;
    numberEntriesAnalysed = 0;
    
    AnalyseRunSAMPIC
    
    %padTemp.yc = 100-pad.yc-53;
    %padTemp.xc = 100-pad.xc-60;
    padTemp.yc = pad.yc;
    padTemp.xc = pad.xc;
    padTemp.padID = getPadForChannelNumber(chPos);
    padTemp.rmsSampling = padRMSSampling;
    padTemp.SATMeanSampling = padSATMeanSampling;
    padTemp.mcpAmpSampling = padMeanMCPAmpSampling;
    padTemp.mmAmpSampling = padMeanMMAmpSampling;
    padTemp.numberEntriesSampling = numberEntriesAnalysed
    padTemp.areaMaps = area;

    

    %if padTemp.padID>0 && numberEntriesAnalysed>0
        padResultsArray = [padResultsArray;padTemp];
   % end
    
end
padResults = padResultsArray;
save(['\\eosproject-smb\eos\project\p\picosec\testbeam\2022_October_h4\Results\Run' run.id '-' run.oscilloscope '\padResults.mat'], 'padResults'); % save interesting data for plots

PlotAllAnalysedPads

AnalyseCombinedChannels