function padID = getPadForChannelNumber(chID, runNumber)
    run.id = runNumber;    
    setMapping;
    
    padID = -1;
    for k=1:length(mapping)
        if mapping(k,2) == chID
           padID = mapping(k,1);
           break;
        end
    end
    
end