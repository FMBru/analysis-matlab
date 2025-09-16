function padID = getPadForChannelNumberMediumGranularity(chID, runNumber)
    run.id = runNumber;    
    setMappingMediumGranularity;
    
    padID = -1;
    mapping;
    for k=1:length(mapping)
        if mapping(k,2) == chID
           padID = mapping(k,1);
           break;
        end
    end
    
end