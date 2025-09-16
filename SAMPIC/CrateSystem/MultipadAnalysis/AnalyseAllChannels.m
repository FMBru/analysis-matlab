clearvars -except runIDsToAnalyse channelsEnabled channelToAnalyse run MM_data MCP_data time_diff time_diff_sigmoid MCP_maxy MM_maxy trackerX trackerY eventIDArray padIDArray MM_epeak MCP_epeak MM_riseTime

for runPos=1:length(runIDsToAnalyse)
    
    run.id = convertStringsToChars(runIDsToAnalyse(runPos));
    
    run.id = convertStringsToChars(runIDsToAnalyse(runPos));
    runTitleString = ['Multipad - July 2025 - Run ' run.id];
    run.name = run.id;
    run.oscilloscope = 'SAMPIC';
    
    
    padResultsArray = [];
    channelToAnalyse=[12 13 14 22 23 24 32 33 34];
    %channelToAnalyse=32:48;
    tracker.dutIndex = 1;
    run.name = run.id;
    run.oscilloscope = 'SAMPIC';
    
    addpath '\\eosproject-smb\eos\project\p\picosec\testbeam\2025_July_h4\Matlab\CommonFunctions';
    addpath '\\eosproject-smb\eos\project\p\picosec\testbeam\2025_July_h4\Matlab\SAMPIC\CrateSystem\Functions';
    
    dutChannelArrayToAnalyse = channelToAnalyse;
    %dutChannelArray = [50;40;30;20]; %override which channels to analyse
    shouldSaveMat=false;
    %dutChannelArrayToAnalyse=27
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
        numberEntriesAnalysedFullPad = 0;
        numberEntriesAnalysedFullPadGlblCut = 0;
        
        bgAvg = 0;
        bgRMS = 0;
        area = 0;
        padRiseTime = 0;
        padTWCorrectionEpeak = 0;
        padTWCorrectionSAT = 0;
        padTWCorrection = 0;
        
        %AnalyseRunSAMPIC
        %analkyse single pad - estimate TW correctin parameters
        AnalyseRunSAMPICFromVariables %pre-loaded only required variables from loadVariablesFromFiles
        
        %padTemp.yc = 100-pad.yc-53;
        %padTemp.xc = 100-pad.xc-60;
        padTemp.yc = pad.yc;
        padTemp.xc = pad.xc;
        padTemp.padID = getPadForChannelNumber(chPos,run.id);
        padTemp.rmsSampling = padRMSSampling;
        padTemp.SATMeanSampling = padSATMeanSampling;
        padTemp.mcpAmpSampling = padMeanMCPAmpSampling;
        padTemp.mmAmpSampling = padMeanMMAmpSampling;
        padTemp.numberEntriesSampling = numberEntriesAnalysed
        padTemp.numberEntriesFullPad = numberEntriesAnalysedFullPad
        padTemp.numberEntriesFullPadGlblCut = numberEntriesAnalysedFullPadGlblCut
%         padTemp.bgAvg = bgAvg;
%         padTemp.bgRMS = bgRMS;
        padTemp.areaMaps = area;
        padTemp.riseTime = padRiseTime;
        padTemp.twCorrEpeak =  padTWCorrectionEpeak;
        padTemp.twCorrSAT =  padTWCorrectionSAT;
        padTemp.twCorr = padTWCorrection;
        
        
        %if padTemp.padID>0 && numberEntriesAnalysed>0
        padResultsArray = [padResultsArray;padTemp];
        % end
        
    end
    padResults = padResultsArray;
    save(['\\eosproject-smb\eos\project\p\picosec\testbeam\2025_July_h4\Results\Run' run.id '-' run.oscilloscope '\padResults.mat'], 'padResults'); % save interesting data for plots
    
    PlotAllAnalysedPads
end
