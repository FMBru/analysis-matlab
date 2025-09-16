run.id = '002';

run.path=['\\eosproject-smb\eos\project\p\picosec\testbeam\2022_May_h4\SAMPIC\TestRun' run.id '\'];
SAMPICFolderInfo = dir(run.path);

foundTriggerFile = false;
foundAsciiDataFile = false;
foundBinaryDataFile = false;

dataFilePath = '';
triggerFilePath = '';

for i=1:length(SAMPICFolderInfo)
    fileInfo = SAMPICFolderInfo(i);
    if contains(fileInfo.name,'SAMPIC_Trigger_Data')    %check if trigger file
        triggerFilePath=[run.path fileInfo.name];
           if exist(triggerFilePath,'file')==2 
                       foundTriggerFile= true;

           end
    elseif contains(fileInfo.name,'.dat')                %check if ascii data file
            dataFilePath=[run.path fileInfo.name];
           if exist(dataFilePath,'file')==2 
                       foundAsciiDataFile= true;
           end
    elseif contains(fileInfo.name,'.bin')               %check if binary data file
            dataFilePath=[run.path fileInfo.name];
           if exist(dataFilePath,'file')==2 
                       foundBinaryDataFile= true;
           end
    end
end

%read in trigger file
if 

fid = fopen('/Users/Florian/cernbox/SAMPIC/Run174/Run174.dat');

%init array for file reading
tline = fgetl(fid);
lineNumber = 1;
lineString = '';
samplingFrequency = 1;
timestep = 1;
hits = [];
firstCell0Time = 0;
cell0TimeArray = [];
unixTimeArray = [];

while ischar(tline) && lineNumber<100000000000
    
    if mod(lineNumber,1000)==0
        lineNumber
    end
       % disp(tline);
        lineString = tline;
        startString = lineString(1:3);
        if startString=='==='
           %is comment
           %disp('comment');
           if contains(lineString,'SAMPLING FREQUENCY')
                samplingFrequencyLineSplit = split(lineString);
                samplingFrequency = str2num(samplingFrequencyLineSplit{4})*1000000 ;   
                timestep = 1/samplingFrequency;
           end
        else
            %disp('dataLine');
            %disp(startString);
           %is data line
           if startString=='Hit'
               hitLine = lineString;
               chLine = fgetl(fid);
               samplesLine = fgetl(fid);
               
               hitLineSplit = split(hitLine);
               hitNumber = str2num(hitLineSplit{2});
               hitUnixTime = str2num(hitLineSplit{5});
               
                chLineSplit = split(chLine);
                chNumber = str2num(chLineSplit{2}) ;   
                cell0Time = str2num(chLineSplit{4}) ;  
                
                if firstCell0Time == 0
                   firstCell0Time = cell0Time; 
                end
                
                samplesLineSplit = split(samplesLine);
                numberDataSampels = length(samplesLineSplit)-2;
                %samplesVector = str2num(cell2mat(samplesLineSplit(2:(length(samplesLineSplit)-1))));
               
                sampleVector = [];
                for pos=1:(length(samplesLineSplit)-2)
                    pickIndex = pos+1;
                    sampleVector(pos) = str2num(samplesLineSplit{pickIndex});
                end
                
                reducedCell0Time = cell0Time-firstCell0Time;
                
                numberSamples = length(sampleVector);
                waveform = zeros(numberSamples,2);
                waveform(1:numberSamples,1) = (1:numberSamples)-1;
                waveform(1:numberSamples,1) = waveform(1:numberSamples,1)*timestep+cell0Time*0.000000001;
                waveform(1:numberSamples,2) = sampleVector;
                
                cell0TimeArray = [cell0TimeArray;reducedCell0Time];
                unixTimeArray = [unixTimeArray;hitUnixTime];
                
               hit.cell0Time = reducedCell0Time;
               hit.unixTimestamp = hitUnixTime; %or hitUnixTime
               hit.ch = chNumber;
               hit.hit = hitNumber;
               hit.waveform = waveform;
               
               hits = [hits;hit];
           end
        end

    
    tline = fgetl(fid);
    lineNumber = lineNumber+1;
end
fclose(fid);

numberHits = length(hits)

%CH0: MPC
%CH1: Detector
%CH2: SmallAreaTrigger
%{
       timeDifferences = [];

for i = 1:numberHits
   hit = hits(i);
   if hit.ch == 1

       
       detectorHitTime = hit.cell0Time;
      %is detector signal, go through others to find match
      for j = 1:numberHits
          testHit = hits(j);
          if testHit.ch==0 && i ~= j 
              testHitTimestamp = testHit.cell0Time;
              timeDiff = testHitTimestamp-detectorHitTime;
              timeDifferences = [timeDifferences;timeDiff];
              if timeDiff > -0.001 && timeDiff < 0.001
                 %select as matching events
                 
              end
          end
      end
   end
end

%}