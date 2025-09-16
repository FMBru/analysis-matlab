        fid = fopen(dataFilePathForLatency);
        
        for i=1:100
fseek(fid,20,'bof');

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
            
            
fclose(fid);