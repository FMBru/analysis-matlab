tracker.path = ['\\eosproject-smb\eos\project\p\picosec\testbeam\2024_June_h4\tracker\reconstructed\asciiRun192.dat'];

    trackerFile = fopen(tracker.path,'rt');
    triggerFilePath=[run.path '\' baseName '\' baseName '_trigger_data.bin'];
    D = textscan(trackerFile, '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f', 'Delimiter','\t', 'HeaderLines',2, 'CollectOutput',1);
    tracker.data = cell2mat(D);


    dataX = tracker.data(:,4);
        dataY = tracker.data(:,5);

        meanX = mean(dataX)        
        meanY = mean(dataY)