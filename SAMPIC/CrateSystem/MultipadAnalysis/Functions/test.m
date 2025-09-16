  
channelsEnabled = [101;102;104;105;107;108;109;1;2;3;4;5;6;7;8;9;10;11;12];
chCountsArrayNumbers = [1:19];

VisualiseMediumGranularity(channelsEnabled,chCountsArrayNumbers,'Hitmap Multipad Picosec',[store_folder '\Run' run.id '_Hitmap.png'],'Hits','%0.1f %',0,100);
