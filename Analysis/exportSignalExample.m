addpath '\\eosproject-smb\eos\project\p\picosec\testbeam\2024_September_h4\Matlab\CommonFunctions'       
ch_mcp_str='\\eosproject-smb\eos\project\p\picosec\testbeam\2024_September_h4\Pool3\Run339\C3--Trace--00022.trc'; 
ch_mcp = ReadLeCroyBinaryWaveform(ch_mcp_str);

x = ch_mcp.x(:,1);
y = ch_mcp.y(:,1);

plot(x,y)
waveform = [x y];
dlmwrite('exampleSignal.csv',waveform)