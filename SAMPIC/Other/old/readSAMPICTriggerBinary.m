function [eventCounterEntries] = readSAMPICTriggerBinary(triggerFile)
%triggerFile = triggerFilePath;
fileID = fopen(triggerFile);
str_disp=sprintf('Reading in trigger file');
disp(str_disp);
%tic

%trigger file to read
%triggerFilePath = '\\eosproject-smb\eos\project\p\picosec\testbeam\2021_October_h4\SAMPIC\Run296\Run_SAMPIC_296_Data_10_30_2021_9h_21min_Ascii/Run_SAMPIC_Trigger_Data_10_30_2021_9h_21min.bin';

%determine number entries
fid = fileID;
fseek(fid, 0, 'eof');
filesize = ftell(fid);
numberEntries = filesize/8
frewind(fid)
eventnumber = numberEntries;

entries = [];
event_id_prev = -1;
event_id_ov = 0;

timestamp_prev = -1;
timestamp_ov = 0;

%loop through events
for i = 1:eventnumber
    
    % - Byte 0: TriggerID from FPGA ( only the LSB byte)
    
    %- Byte 1 (LSB) to Byte 2 (MSB):  TriggerIDFrom External Trigger (16 Bits)
    
    %- Bytes 3(LSB) to 7 (MSB) : TimeStamp from the FPGA (on 40 Bits).
    
    
    
    %extract byte by byte
    TriggerIDFPGA = fread(fid,1,'ubit8',0);
    TriggerIDSRSRaw = fread(fid,1,'ubit16',0);
    %byte1 = fread(fid,1,'ubit8',0);
    %byte2 = fread(fid,1,'ubit8',0);
    %byte3 = fread(fid,1,'ubit8',0);
    %byte4 = fread(fid,1,'ubit8',0);
    %byte5 = fread(fid,1,'ubit8',0);
    %timestamp = 256*256*256*256*byte5 + 256*256*256*byte4 + 256*256*byte3 + 256*byte2 + byte1;
    %timestamp = 7.5*fread(fid,1,'ubit40',0);
    periodFactor = (64*1000)/8512;
    
    timestampRaw = fread(fid,1,'ubit40',0);
    
    % count overflows
    if(timestampRaw < timestamp_prev)
        timestamp_ov = timestamp_ov + 1;
    end
    timestamp_prev = timestampRaw;
    
    timestamp = timestamp_ov * 1099511627776 + timestampRaw;
    timestamp = timestamp* periodFactor;
    
    % count overflows
    if(TriggerIDSRSRaw < event_id_prev)
        event_id_ov = event_id_ov + 1;
    end
    event_id_prev = TriggerIDSRSRaw;
    TriggerIDSRS = event_id_ov * 65536 + TriggerIDSRSRaw;
    %timestamp = fread(fid,1,'ubit40',0);
    %timestamp = bitshift(timestamp,2)
    %timestampShifted = bitshift(timestamp,8);
    
    entry = [TriggerIDSRS timestamp];
    %entry = [i TriggerIDFPGA TriggerIDSRS timestamp];

    entries(i,:) = entry;

    if mod(i,50000)==0
        str_disp=sprintf('Read trigger %d', i);
        disp(str_disp);
        %toc 
        %tic
%         if i>200000 % DEBUGGING to process only small number of triggers
%             break;
%         end
    end
    
    
    %pause(1);
end


fclose(fid);

%plot decoded trigger ID vs timestamp
shouldPlot = 0;

if shouldPlot ==1
    close all
    figure(1);
    plot(entries(:,4), entries(:,3),'.');
end
eventCounterEntries = entries;
end

