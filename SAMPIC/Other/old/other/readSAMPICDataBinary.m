filePath = '\\eosproject-smb\eos\project\p\picosec\testbeam\2022_May_h4\SAMPIC\Run420\Run420.bin';

fid = fopen(filePath);
%A = fread(fileID)

numberHeaderLines = 7;
for i=1:numberHeaderLines
tline = fgets(fid)
end
pause(1)

eventnumber = 10;

for i = 1:eventnumber


hitNumber = fread(fid,1,'int32',0)  
epochTime = fread(fid,1,'double',0) 
channel = fread(fid,1,'int32',0)  
TriggerCellTimeInstant = fread(fid,1,'double',0)  
CFDTimeInstant = fread(fid,1,'double',0)  
Baseline = fread(fid,1,'float64',0)  
PeakValue = fread(fid,1,'float64',0)  
Amplitude = fread(fid,1,'float32',0)  

dataSize = fread(fid,1,'int32',0) 

dataLength = dataSize;


dataVec = zeros(dataLength,1);
for pos=1:dataLength
    dataVec(pos) = fread(fid,1,'float32',0)  ;
end
close all

plot(dataVec,'.-')

pause(2);

end

fclose(fid);




