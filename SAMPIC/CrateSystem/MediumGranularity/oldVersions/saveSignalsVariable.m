
numberVars = floor(length(signals)/100000);


for k = 1:numberVars
    startIndex = 100000*(k-1)+1;
    endIndex = k*100000;
    signalsSegment = signals(startIndex:endIndex);
    chIDsArraySegment = chIDsArray(startIndex:endIndex);
    chTimesArraySegment = chTimesArray(startIndex:endIndex);
    chIndexArraySegment = chIndexArray(startIndex:endIndex);
    %save signals
    storeMatfilePath = [storeMatfileFolder '\signals_' int2str(k) '.mat'];
    mkdir(storeMatfileFolder)
    m = matfile(storeMatfilePath,'Writable',true);
    %
    save(storeMatfilePath,'signalsSegment');
    save(storeMatfilePath,'chIDsArraySegment','-append');
    save(storeMatfilePath,'chTimesArraySegment','-append');
    save(storeMatfilePath,'chIndexArraySegment','-append');

end