function chID = getChannelForPadNumberHighGranularity(padID)
    setMappingHighGranularity
    
    chID = -1;
    
    for k=1:length(mapping)
        if mapping(k,1) == padID
           chID = mapping(k,2);
           break;
        end
    end
    
end