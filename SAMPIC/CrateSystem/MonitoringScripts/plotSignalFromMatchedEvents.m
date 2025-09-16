close all

doPlot = false
startTimeDiffMCPArray = [];
startTimeDiffScintArray = [];

counter = 0;

for eventIndex = 1:size(matchedEvents,1)
    matchedEvents(eventIndex).matchedScint = 0;

    
    close all
    
    signal = matchedEvents(eventIndex).signal;
    mcp = matchedEvents(eventIndex).mcp;
    scint = matchedEvents(eventIndex).scint;
    
    signalWaveform = signal.waveform;
    mcpWaveform = mcp.waveform;
    scintWaveform = scint.waveform;
    
    
    startTime = signalWaveform(1,1);
    startTimeMCP = mcpWaveform(1,1);
    startTimeScint = scintWaveform(1,1);
    
    startTimeDiffMCP = 1e9*(startTime-startTimeMCP);
    startTimeDiffScint = 1e9*(startTime-startTimeScint);
    if startTimeDiffMCP>-120 && startTimeDiffMCP<-100
        
        
        signalWaveform(:,1) = (signalWaveform(:,1)-startTime)*1e9;
        mcpWaveform(:,1) = (mcpWaveform(:,1)-startTimeMCP)*1e9;
        scintWaveform(:,1) = (scintWaveform(:,1)-startTime)*1e9;
        
        if startTimeDiffScint>-90 && startTimeDiffScint<-70
                %event selected for MCP and scint
               startTimeDiffMCPArray = [startTimeDiffMCPArray;startTimeDiffMCP];
                startTimeDiffScintArray = [startTimeDiffScintArray;startTimeDiffScint];
                matchedEvents(eventIndex).matchedScint = 1;
        counter = counter+1;
            if doPlot
                plot(signalWaveform(:,1),signalWaveform(:,2),'b.-'); hold on
                plot(mcpWaveform(:,1),mcpWaveform(:,2),'r.-');
                plot(scintWaveform(:,1),scintWaveform(:,2),'g.-');
                xlabel("Time (ns)");
                ylabel("Amplitude (V)");
                %legend("Detector signal","MCP","ScintillatorTrigger");
                pause (1);
                
            end
            
        end
        
    end
    
end