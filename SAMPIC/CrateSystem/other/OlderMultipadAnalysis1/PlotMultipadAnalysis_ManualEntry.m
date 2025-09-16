channelsEnabled = [12,19,82,89,56];
dataArray = [21.5,22.1,19.1,19.1,29.1]; %double DLC CsI Sept 24 - Time res single Gauss 275/470

VisualiseMultipad(channelsEnabled,dataArray,'Time resolution Double Picosec','\\eosproject-smb\eos\project\p\picosec\testbeam\2024_September_h4\Results\DoubleDLC-CsI-TimeRes.png','TimeR','%0.1f ps',0,0);


channelsEnabled = [12,19,82,89,56];
dataArray = [97.2,91.8,104.2,99.3,50.6]; %double DLC CsI Sept 24 - Mean amp 275/470

VisualiseMultipad(channelsEnabled,dataArray,'Signal amplitude Double Picosec','\\eosproject-smb\eos\project\p\picosec\testbeam\2024_September_h4\Results\DoubleDLC-CsI-MeanAmp.png','Amp','%0.1f mV',0,0);

